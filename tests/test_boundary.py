import unittest
from math import pi

import numpy as np

from mindlin_plate import (
    MindlinMaterial,
    assemble_boundary_load,
    edge_angle,
    exterior_edges,
    global_dof_transform,
    local_edge_constraints,
    nodal_rotation_transform,
    q4_edge_consistent_load,
    q4_element_matrices,
    rectangular_mesh,
    select_exterior_edges,
    solve_with_local_constraints,
    symmetry_dofs,
    transform_system_to_local,
)


COORDS = np.array(
    [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]], dtype=float
)


class TestNaturalBoundaryLoads(unittest.TestCase):
    def test_constant_edge_resultants(self):
        expected = np.array([10.0, 2.0, -3.0])
        for edge in range(4):
            force = q4_edge_consistent_load(
                COORDS, edge, transverse_shear=10.0, moment=(2.0, -3.0)
            )
            np.testing.assert_allclose(force.reshape(4, 3).sum(axis=0), expected)

    def test_exterior_edge_selection_and_assembly(self):
        mesh = rectangular_mesh(2.0, 1.0, 2, 1)
        self.assertEqual(len(exterior_edges(mesh)), 6)
        top = select_exterior_edges(mesh, lambda midpoint: np.isclose(midpoint[1], 1.0))
        self.assertEqual(len(top), 2)
        force = assemble_boundary_load(
            mesh, top, transverse_shear=4.0, moment=(1.5, -0.5)
        )
        np.testing.assert_allclose(
            force.reshape(-1, 3).sum(axis=0), [8.0, 3.0, -1.0]
        )


class TestBoundaryCoordinates(unittest.TestCase):
    def test_rotation_and_energy_invariance(self):
        angle = edge_angle([0.0, 0.0], [1.0, 1.0])
        self.assertAlmostEqual(angle, pi / 4.0)
        nodal = nodal_rotation_transform(angle)
        np.testing.assert_allclose(nodal.T @ nodal, np.eye(3), atol=1e-14)

        transform = global_dof_transform(4, {node: angle for node in range(4)})
        material = MindlinMaterial(210e9, 0.3, 0.1)
        matrix = q4_element_matrices(COORDS, material, "mitc4").total
        local_matrix, _ = transform_system_to_local(matrix, np.zeros(12), transform)
        local_dofs = np.arange(12, dtype=float) / 10.0
        global_dofs = transform @ local_dofs
        self.assertAlmostEqual(
            global_dofs @ matrix @ global_dofs,
            local_dofs @ local_matrix @ local_dofs,
            delta=1e-5,
        )

    def test_local_constraint_kinds_and_solution(self):
        self.assertEqual(set(local_edge_constraints([0, 2], "clamped")), {0, 1, 2, 6, 7, 8})
        self.assertEqual(set(local_edge_constraints([0, 2], "hard_simply_supported")), {0, 1, 6, 7})
        self.assertEqual(set(local_edge_constraints([0, 2], "symmetry")), {2, 8})

        angle = pi / 6.0
        solution = solve_with_local_constraints(
            np.eye(3), np.array([1.0, 2.0, 3.0]), 1, {0: angle}, {2: 0.0}
        )
        normal = np.array([-np.sin(angle), np.cos(angle)])
        self.assertAlmostEqual(normal @ solution.displacement[1:3], 0.0, places=13)

    def test_axis_aligned_symmetry_dofs(self):
        mesh = rectangular_mesh(1.0, 1.0, 2, 2)
        self.assertEqual(set(symmetry_dofs(mesh, "x_min")), {1, 10, 19})
        self.assertEqual(set(symmetry_dofs(mesh, "y_max")), {20, 23, 26})


if __name__ == "__main__":
    unittest.main()

