import json
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

from mindlin_plate.api import create_app
from mindlin_plate.service import (
    MAX_SERVICE_DOFS,
    PlateInputError,
    serialize_analysis,
    solve_plate_case,
)


TEMPLATE = (
    Path(__file__).resolve().parents[1]
    / "mindlin_plate"
    / "templates"
    / "rectangular_plate.json"
)


def small_case():
    case = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    case["mesh"]["elements_x"] = 2
    case["mesh"]["elements_y"] = 2
    return case


class TestPlateService(unittest.TestCase):
    def test_template_solves_and_balances(self):
        result = serialize_analysis(solve_plate_case(small_case()))
        summary = result["summary"]
        self.assertEqual(summary["plate_method"], "K")
        self.assertEqual(summary["node_count"], 9)
        self.assertEqual(len(result["element_results"]), 4)
        self.assertLess(abs(summary["transverse_equilibrium_error"]), 1e-8)
        self.assertGreater(summary["peak_deflection"], 0.0)

    def test_unconstrained_custom_case_is_rejected(self):
        case = small_case()
        case["boundary"] = {"type": "custom", "prescribed_dofs": {}}
        with self.assertRaisesRegex(PlateInputError, "requires prescribed_dofs"):
            solve_plate_case(case)

    def test_dense_problem_limit_is_checked_before_assembly(self):
        case = small_case()
        case["mesh"]["elements_x"] = 29
        case["mesh"]["elements_y"] = 29
        with self.assertRaisesRegex(
            PlateInputError, f"dense service limit is {MAX_SERVICE_DOFS}"
        ):
            solve_plate_case(case)


class TestPlateApi(unittest.TestCase):
    def test_template_solve_results_and_png_artifacts(self):
        case = small_case()
        case["postprocess"]["plots"] = ["deflection", "stress_top"]
        with tempfile.TemporaryDirectory() as directory:
            client = TestClient(create_app(directory))
            self.assertEqual(client.get("/health").json(), {"status": "ok"})
            template_response = client.get("/api/v1/templates/rectangular-plate")
            self.assertEqual(template_response.status_code, 200)
            self.assertEqual(template_response.json()["mesh"]["type"], "rectangular")

            solve_response = client.post("/api/v1/analyses", json=case)
            self.assertEqual(solve_response.status_code, 201, solve_response.text)
            payload = solve_response.json()
            self.assertEqual(payload["status"], "completed")
            self.assertEqual(len(payload["images"]), 2)

            for artifact in payload["images"]:
                image = client.get(artifact["url"])
                self.assertEqual(image.status_code, 200)
                self.assertEqual(image.headers["content-type"], "image/png")
                self.assertTrue(image.content.startswith(b"\x89PNG\r\n\x1a\n"))

            result = client.get(f"/api/v1/analyses/{payload['analysis_id']}")
            self.assertEqual(result.status_code, 200)
            self.assertEqual(result.json()["analysis_id"], payload["analysis_id"])

    def test_invalid_custom_boundary_returns_422(self):
        case = small_case()
        case["boundary"] = {"type": "custom", "prescribed_dofs": {}}
        case["postprocess"]["plots"] = []
        with tempfile.TemporaryDirectory() as directory:
            response = TestClient(create_app(directory)).post(
                "/api/v1/analyses", json=case
            )
        self.assertEqual(response.status_code, 422)
        self.assertIn("prescribed_dofs", response.json()["detail"])


if __name__ == "__main__":
    unittest.main()
