import unittest
from math import cos, pi, sin, sqrt

import numpy as np

from mindlin_plate import (
    MindlinMaterial,
    assemble_plate_system,
    assemble_system,
    dkq_element_matrices,
    dkq_kinematic_matrices,
    element_response,
    hard_simply_supported_dofs,
    plate_characteristic_length,
    q8_shape_functions,
    rectangular_mesh,
    select_plate_method,
    simply_supported_sinusoidal_solution,
    sinusoidal_load,
    solve_dirichlet,
)


COORDS = np.array(
    [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]], dtype=float
)


def nodal_quadratic(coordinates, kappa_x, kappa_y, kappa_xy):
    result = np.zeros(12, dtype=float)
    for node, (x, y) in enumerate(coordinates):
        result[3 * node : 3 * node + 3] = (
            0.5 * kappa_x * x * x
            + 0.5 * kappa_y * y * y
            + kappa_xy * x * y,
            kappa_x * x + kappa_xy * y,
            kappa_y * y + kappa_xy * x,
        )
    return result


class TestDiscreteKirchhoffCore(unittest.TestCase):
    def test_q8_partition_and_gradient_partition(self):
        values, gradients = q8_shape_functions(0.23, -0.37)
        self.assertAlmostEqual(values.sum(), 1.0)
        np.testing.assert_allclose(gradients.sum(axis=0), 0.0, atol=1e-14)

    def test_constant_curvature_patch_on_skew_quadrilateral(self):
        coordinates = np.array(
            [[0.0, 0.0], [2.0, 0.1], [2.6, 1.2], [0.2, 1.0]], dtype=float
        )
        dofs = nodal_quadratic(coordinates, 1.7, -0.8, 0.45)
        expected = np.array([1.7, -0.8, 0.9])
        value = 1.0 / sqrt(3.0)
        for xi, eta in ((-value, -value), (value, -value), (0.2, 0.35)):
            _, _, bending = dkq_kinematic_matrices(coordinates, xi, eta)
            np.testing.assert_allclose(bending @ dofs, expected, atol=1e-12)

    def test_free_element_has_only_three_rigid_modes(self):
        material = MindlinMaterial(210e9, 0.3, 0.01)
        element = dkq_element_matrices(COORDS, material)
        np.testing.assert_allclose(element.shear, 0.0)
        eigenvalues = np.linalg.eigvalsh(element.total)
        tolerance = max(float(np.max(np.abs(eigenvalues))) * 1e-10, 1e-8)
        self.assertEqual(np.count_nonzero(np.abs(eigenvalues) <= tolerance), 3)
        self.assertGreaterEqual(eigenvalues[0], -tolerance)

    def test_stiffness_is_invariant_under_planar_rotation(self):
        material = MindlinMaterial(70e9, 0.23, 0.01)
        angle = 0.43
        rotation = np.array(
            [[cos(angle), -sin(angle)], [sin(angle), cos(angle)]], dtype=float
        )
        rotated_coordinates = COORDS @ rotation.T
        original = dkq_element_matrices(COORDS, material).total
        rotated = dkq_element_matrices(rotated_coordinates, material).total
        transform = np.zeros((12, 12), dtype=float)
        for node in range(4):
            transform[3 * node, 3 * node] = 1.0
            transform[
                3 * node + 1 : 3 * node + 3,
                3 * node + 1 : 3 * node + 3,
            ] = rotation
        np.testing.assert_allclose(
            rotated, transform @ original @ transform.T, rtol=1e-13, atol=1e-7
        )

    def test_kirchhoff_response_has_no_transverse_shear(self):
        material = MindlinMaterial(70e9, 0.25, 0.01)
        dofs = nodal_quadratic(COORDS, 0.8, -0.3, 0.2)
        response = element_response(
            COORDS, dofs, material, plate_method="K"
        )
        self.assertEqual(response.plate_method, "K")
        np.testing.assert_allclose(response.shear_strain, 0.0)
        np.testing.assert_allclose(response.shear_force, 0.0)
        np.testing.assert_allclose(response.curvature, [0.8, -0.3, 0.4])


class TestPlateMethodSelection(unittest.TestCase):
    def test_rotated_planform_uses_short_span(self):
        rectangle = np.array(
            [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]], dtype=float
        )
        angle = 0.37
        rotation = np.array(
            [[cos(angle), -sin(angle)], [sin(angle), cos(angle)]], dtype=float
        )
        rotated = rectangle @ rotation.T
        self.assertAlmostEqual(plate_characteristic_length(rotated), 1.0)
        self.assertEqual(select_plate_method(rotated, 0.05).method, "K")
        self.assertEqual(select_plate_method(rotated, 0.051).method, "M")

    def test_auto_assembly_matches_selected_explicit_method(self):
        mesh = rectangular_mesh(1.0, 1.0, 2, 2)
        thin = MindlinMaterial(210e9, 0.3, 0.01)
        thin_auto = assemble_plate_system(mesh, thin, 10.0)
        thin_k, thin_force = assemble_system(
            mesh, thin, 10.0, plate_method="K"
        )
        self.assertEqual(thin_auto.selection.method, "K")
        np.testing.assert_allclose(thin_auto.stiffness, thin_k)
        np.testing.assert_allclose(thin_auto.force, thin_force)

        thick = MindlinMaterial(210e9, 0.3, 0.1)
        thick_auto = assemble_plate_system(mesh, thick, 10.0)
        thick_m, thick_force = assemble_system(
            mesh, thick, 10.0, "mitc4", plate_method="M"
        )
        self.assertEqual(thick_auto.selection.method, "M")
        np.testing.assert_allclose(thick_auto.stiffness, thick_m)
        np.testing.assert_allclose(thick_auto.force, thick_force)

    def test_thin_square_converges_to_kirchhoff_benchmark(self):
        mesh = rectangular_mesh(1.0, 1.0, 4, 4)
        material = MindlinMaterial(210e9, 0.3, 0.01)
        system = assemble_plate_system(
            mesh, material, sinusoidal_load(10e3, 1.0, 1.0)
        )
        solution = solve_dirichlet(
            system.stiffness, system.force, hard_simply_supported_dofs(mesh)
        )
        center = int(np.argmin(np.linalg.norm(mesh.nodes - [0.5, 0.5], axis=1)))
        computed = solution.displacement[3 * center]
        exact = simply_supported_sinusoidal_solution(
            1.0, material, 10e3
        ).bending_deflection
        self.assertGreater(computed / exact, 0.99)
        self.assertLess(computed / exact, 1.01)


if __name__ == "__main__":
    unittest.main()
