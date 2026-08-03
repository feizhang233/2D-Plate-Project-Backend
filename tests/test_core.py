import unittest
from math import sqrt

import numpy as np

from mindlin_plate import (
    MindlinMaterial,
    assemble_system,
    hard_simply_supported_dofs,
    kinematic_matrices,
    mitc4_shear_B,
    q4_consistent_load,
    q4_element_matrices,
    rectangular_mesh,
    simply_supported_sinusoidal_solution,
    sinusoidal_load,
    solve_dirichlet,
)


COORDS = np.array(
    [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]], dtype=float
)


def nodal_field(field):
    result = np.zeros(12)
    for index, (x, y) in enumerate(COORDS):
        result[3 * index : 3 * index + 3] = field(x, y)
    return result


class TestContinuum(unittest.TestCase):
    def test_reference_numerical_example(self):
        material = MindlinMaterial(210e9, 0.3, 0.1)
        answer = simply_supported_sinusoidal_solution(1.0, material, 10e3)
        self.assertAlmostEqual(material.shear_modulus, 80.769230769e9, delta=100.0)
        self.assertAlmostEqual(material.bending_rigidity, 1.923076923e7, delta=1.0)
        self.assertAlmostEqual(answer.bending_deflection, 1.33458e-6, delta=1e-11)
        self.assertAlmostEqual(answer.shear_deflection, 7.52672e-8, delta=1e-12)
        self.assertAlmostEqual(answer.center_deflection, 1.40984e-6, delta=1e-11)

    def test_surface_stress_is_antisymmetric(self):
        material = MindlinMaterial(70e9, 0.25, 0.02)
        curvature = np.array([0.1, -0.05, 0.02])
        top = material.bending_stress(curvature, 0.01)
        bottom = material.bending_stress(curvature, -0.01)
        np.testing.assert_allclose(top, -bottom)

    def test_parabolic_shear_stress_recovers_resultant(self):
        material = MindlinMaterial(70e9, 0.25, 0.02)
        shear_force = np.array([125.0, -80.0])
        half_t = material.thickness / 2.0
        np.testing.assert_allclose(
            material.parabolic_shear_stress(shear_force, half_t), 0.0, atol=1e-12
        )
        points, weights = np.polynomial.legendre.leggauss(3)
        integrated = sum(
            weight
            * material.parabolic_shear_stress(shear_force, half_t * point)
            * half_t
            for point, weight in zip(points, weights)
        )
        np.testing.assert_allclose(integrated, shear_force, rtol=1e-14, atol=1e-12)


class TestQ4(unittest.TestCase):
    def test_partition_and_mapping(self):
        point, _, _ = kinematic_matrices(COORDS, 0.2, -0.4)
        self.assertAlmostEqual(point.shape.sum(), 1.0)
        np.testing.assert_allclose(point.physical, [0.1, -0.2])
        self.assertAlmostEqual(point.det_jacobian, 0.25)

    def test_skew_affine_mapping_reproduces_physical_gradient(self):
        coordinates = np.array([[0.0, 0.0], [2.0, 0.0], [3.0, 1.0], [1.0, 1.0]])
        nodal_values = 2.0 * coordinates[:, 0] - 3.0 * coordinates[:, 1] + 4.0
        for xi, eta in [(-0.3, 0.2), (0.5, -0.6)]:
            point, _, _ = kinematic_matrices(coordinates, xi, eta)
            np.testing.assert_allclose(
                point.physical_gradients.T @ nodal_values, [2.0, -3.0], atol=1e-12
            )

    def test_rigid_and_patch_modes(self):
        modes = [
            nodal_field(lambda x, y: (1.0, 0.0, 0.0)),
            nodal_field(lambda x, y: (2.0 * x, 2.0, 0.0)),
            nodal_field(lambda x, y: (-3.0 * y, 0.0, -3.0)),
        ]
        twist = nodal_field(lambda x, y: (x * y, y, x))
        for xi in (-1 / sqrt(3), 1 / sqrt(3)):
            for eta in (-1 / sqrt(3), 1 / sqrt(3)):
                _, bending, shear = kinematic_matrices(COORDS, xi, eta)
                for mode in modes:
                    np.testing.assert_allclose(bending @ mode, 0.0, atol=1e-12)
                    np.testing.assert_allclose(shear @ mode, 0.0, atol=1e-12)
                np.testing.assert_allclose(bending @ twist, [0.0, 0.0, 2.0], atol=1e-12)
                np.testing.assert_allclose(shear @ twist, 0.0, atol=1e-12)

    def test_mitc4_removes_reference_pure_bending_shear(self):
        curvature = 1.7
        mode = nodal_field(lambda x, y: (0.5 * curvature * x * x, curvature * x, 0.0))
        for xi in (-1 / sqrt(3), 0.2):
            for eta in (-1 / sqrt(3), 0.3):
                np.testing.assert_allclose(mitc4_shear_B(COORDS, xi, eta) @ mode, 0.0, atol=1e-12)


class TestElement(unittest.TestCase):
    def test_constant_load_and_symmetry(self):
        material = MindlinMaterial(210e9, 0.3, 0.1)
        force = q4_consistent_load(COORDS, 12.0)
        np.testing.assert_allclose(force[0::3], 3.0)
        np.testing.assert_allclose(force[1::3], 0.0)
        np.testing.assert_allclose(force[2::3], 0.0)
        for scheme in ("full", "reduced", "mitc4"):
            matrix = q4_element_matrices(COORDS, material, scheme).total
            np.testing.assert_allclose(matrix, matrix.T, rtol=1e-13, atol=1e-7)

    def test_full_pure_bending_spurious_energy_formula(self):
        material = MindlinMaterial(210e9, 0.3, 0.01)
        curvature = 2.3
        mode = nodal_field(lambda x, y: (0.5 * curvature * x * x, curvature * x, 0.0))
        element = q4_element_matrices(COORDS, material, "full")
        bending_energy = 0.5 * mode @ element.bending @ mode
        shear_energy = 0.5 * mode @ element.shear @ mode
        expected_ratio = (
            material.shear_correction
            * (1.0 - material.poisson)
            / 2.0
            * (1.0 / material.thickness) ** 2
        )
        self.assertAlmostEqual(shear_energy / bending_energy, expected_ratio, places=8)
        for scheme in ("reduced", "mitc4"):
            shear = q4_element_matrices(COORDS, material, scheme).shear
            self.assertLess(abs(0.5 * mode @ shear @ mode), shear_energy * 1e-12)

    def test_free_full_element_has_three_rigid_modes(self):
        material = MindlinMaterial(210e9, 0.3, 0.1)
        matrix = q4_element_matrices(COORDS, material, "full").total
        eigenvalues = np.linalg.eigvalsh(matrix)
        tolerance = max(abs(eigenvalues[-1]) * 1e-10, 1e-8)
        self.assertEqual(np.count_nonzero(np.abs(eigenvalues) < tolerance), 3)


class TestAssembly(unittest.TestCase):
    def test_sinusoidal_plate_solves_and_balances(self):
        mesh = rectangular_mesh(1.0, 1.0, 4, 4)
        material = MindlinMaterial(210e9, 0.3, 0.1)
        stiffness, force = assemble_system(
            mesh, material, sinusoidal_load(10e3, 1.0, 1.0), "mitc4"
        )
        solution = solve_dirichlet(stiffness, force, hard_simply_supported_dofs(mesh))
        self.assertLess(np.linalg.norm(solution.reactions[solution.free_dofs], ord=np.inf), 1e-4)
        center = np.argmin(np.linalg.norm(mesh.nodes - [0.5, 0.5], axis=1))
        self.assertGreater(solution.displacement[3 * center], 0.0)
        self.assertAlmostEqual(
            solution.reactions[0::3].sum(), -force[0::3].sum(), delta=1e-5
        )


if __name__ == "__main__":
    unittest.main()
