import unittest

import numpy as np

from mindlin_plate import (
    MindlinMaterial,
    constant_shear_patch_test,
    element_spectrum,
    recover_element_gauss_points,
    rectangular_mesh,
    split_nodal_dofs,
    square_distortion_comparison,
    square_thickness_scan,
)


COORDS = np.array(
    [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]], dtype=float
)


class TestVerificationGates(unittest.TestCase):
    def test_zero_modes_distinguish_stable_mitc_from_reduced_baseline(self):
        material = MindlinMaterial(210e9, 0.3, 0.1)
        spectra = {
            scheme: element_spectrum(COORDS, material, scheme)
            for scheme in ("full", "reduced", "mitc4")
        }
        self.assertTrue(spectra["full"].has_only_rigid_zero_modes)
        self.assertTrue(spectra["mitc4"].has_only_rigid_zero_modes)
        self.assertEqual(spectra["reduced"].zero_modes, 5)

    def test_constant_shear_patch_on_skew_element(self):
        skew = np.array([[0.0, 0.0], [2.0, 0.0], [2.7, 1.0], [0.3, 1.0]])
        for scheme in ("full", "reduced", "mitc4"):
            result = constant_shear_patch_test(skew, (1.2, -0.7), scheme)
            self.assertTrue(result.passed)

    def test_complete_thickness_scan_reaches_1e_minus_4(self):
        rows = square_thickness_scan(
            1.0, 210e9, 0.3, 10e3, ratios=(1e-1, 1e-4), divisions=4
        )
        self.assertEqual(rows[-1].thickness_ratio, 1e-4)
        normalized = rows[-1].normalized_deflections
        self.assertLess(normalized["full"], 1e-4)
        self.assertGreater(normalized["mitc4"], 0.95)

    def test_distorted_mesh_comparison(self):
        material = MindlinMaterial(210e9, 0.3, 0.01)
        result = square_distortion_comparison(
            1.0, material, 10e3, divisions=4, shear_scheme="mitc4"
        )
        self.assertGreater(result.minimum_gauss_det_jacobian, 0.0)
        self.assertLess(result.relative_change, 0.05)


class TestPostprocessing(unittest.TestCase):
    def test_dof_split_and_gauss_recovery(self):
        mesh = rectangular_mesh(1.0, 1.0, 1, 1)
        displacement = np.arange(mesh.ndof, dtype=float) * 1e-6
        w, theta = split_nodal_dofs(displacement)
        np.testing.assert_allclose(w, displacement[0::3])
        np.testing.assert_allclose(theta, displacement.reshape(-1, 3)[:, 1:])
        results = recover_element_gauss_points(
            mesh, displacement, MindlinMaterial(70e9, 0.25, 0.02), "mitc4"
        )
        self.assertEqual(len(results), 4)
        for result in results:
            np.testing.assert_allclose(
                result.response.stress_top, -result.response.stress_bottom
            )


if __name__ == "__main__":
    unittest.main()

