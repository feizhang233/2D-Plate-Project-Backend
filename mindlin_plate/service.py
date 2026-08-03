"""Application service that turns JSON-like plate data into core results.

This module deliberately has no FastAPI or Pydantic dependency.  It is the
small adapter between transport-layer JSON and the NumPy mathematics, which
makes the same workflow usable from the HTTP API, tests, or a Python script.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sin
from typing import Any, Mapping

import numpy as np

from .assembly import (
    AssembledPlateSystem,
    LinearSolution,
    assemble_plate_system,
    clamped_dofs,
    hard_simply_supported_dofs,
    soft_simply_supported_dofs,
    solve_dirichlet,
)
from .material import MindlinMaterial
from .mesh import Mesh, rectangular_mesh
from .postprocess import RecoveredElement, recover_element_centers


MAX_SERVICE_DOFS = 2500


class PlateInputError(ValueError):
    """Raised when otherwise valid JSON cannot describe a solvable case."""


@dataclass(frozen=True)
class PlateAnalysis:
    """The complete in-memory result of one plate analysis."""

    name: str
    mesh: Mesh
    material: MindlinMaterial
    system: AssembledPlateSystem
    solution: LinearSolution
    recovered: list[RecoveredElement]
    shear_scheme: str
    boundary_type: str


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise PlateInputError(f"{name} must be a JSON object")
    return value


def _make_mesh(data: Mapping[str, Any]) -> Mesh:
    kind = str(data.get("type", "rectangular")).strip().lower()
    try:
        if kind == "rectangular":
            return rectangular_mesh(
                float(data["length_x"]),
                float(data["length_y"]),
                int(data["elements_x"]),
                int(data["elements_y"]),
            )
        if kind == "custom":
            if data.get("nodes") is None or data.get("elements") is None:
                raise PlateInputError("custom mesh requires nodes and elements")
            return Mesh(
                np.asarray(data["nodes"], dtype=float),
                np.asarray(data["elements"], dtype=np.int64),
            )
    except KeyError as exc:
        raise PlateInputError(f"mesh is missing required field: {exc.args[0]}") from exc
    except (TypeError, ValueError) as exc:
        raise PlateInputError(f"invalid mesh: {exc}") from exc
    raise PlateInputError("mesh.type must be 'rectangular' or 'custom'")


def _make_material(data: Mapping[str, Any]) -> MindlinMaterial:
    try:
        return MindlinMaterial(
            young=float(data["young"]),
            poisson=float(data["poisson"]),
            thickness=float(data["thickness"]),
            shear_correction=float(data.get("shear_correction", 5.0 / 6.0)),
        )
    except KeyError as exc:
        raise PlateInputError(f"material is missing required field: {exc.args[0]}") from exc
    except (TypeError, ValueError) as exc:
        raise PlateInputError(f"invalid material: {exc}") from exc


def _make_load(data: Mapping[str, Any], mesh: Mesh):
    kind = str(data.get("type", "uniform")).strip().lower()
    try:
        magnitude = float(data["magnitude"])
    except KeyError as exc:
        raise PlateInputError("load is missing required field: magnitude") from exc
    except (TypeError, ValueError) as exc:
        raise PlateInputError(f"invalid load magnitude: {exc}") from exc
    if kind == "uniform":
        return magnitude
    if kind == "sinusoidal":
        x_min, y_min = mesh.nodes.min(axis=0)
        x_max, y_max = mesh.nodes.max(axis=0)
        length_x = float(x_max - x_min)
        length_y = float(y_max - y_min)
        if length_x <= 0.0 or length_y <= 0.0:
            raise PlateInputError("sinusoidal load requires non-zero x and y spans")

        def load_function(x: float, y: float) -> float:
            return magnitude * sin(pi * (x - x_min) / length_x) * sin(
                pi * (y - y_min) / length_y
            )

        return load_function
    raise PlateInputError("load.type must be 'uniform' or 'sinusoidal'")


def _make_constraints(data: Mapping[str, Any], mesh: Mesh) -> tuple[str, dict[int, float]]:
    kind = str(data.get("type", "hard_simply_supported")).strip().lower()
    factories = {
        "clamped": clamped_dofs,
        "soft_simply_supported": soft_simply_supported_dofs,
        "hard_simply_supported": hard_simply_supported_dofs,
    }
    if kind in factories:
        constraints = factories[kind](mesh)
    elif kind == "custom":
        raw = data.get("prescribed_dofs")
        if not isinstance(raw, Mapping) or not raw:
            raise PlateInputError("custom boundary requires prescribed_dofs")
        try:
            constraints = {int(index): float(value) for index, value in raw.items()}
        except (TypeError, ValueError) as exc:
            raise PlateInputError(f"invalid prescribed_dofs: {exc}") from exc
    else:
        raise PlateInputError(
            "boundary.type must be clamped, soft_simply_supported, "
            "hard_simply_supported, or custom"
        )
    if not constraints:
        raise PlateInputError("at least one displacement constraint is required")
    return kind, constraints


def solve_plate_case(case: Mapping[str, Any]) -> PlateAnalysis:
    """Validate JSON-like input, solve it, and recover element-center fields."""

    root = _mapping(case, "request")
    mesh = _make_mesh(_mapping(root.get("mesh"), "mesh"))
    if mesh.ndof > MAX_SERVICE_DOFS:
        raise PlateInputError(
            f"mesh has {mesh.ndof} dofs; the dense service limit is "
            f"{MAX_SERVICE_DOFS} dofs"
        )
    material = _make_material(_mapping(root.get("material"), "material"))
    load = _make_load(_mapping(root.get("load"), "load"), mesh)
    boundary_type, prescribed = _make_constraints(
        _mapping(root.get("boundary"), "boundary"), mesh
    )
    options = _mapping(root.get("analysis", {}), "analysis")
    plate_method = str(options.get("plate_method", "auto"))
    shear_scheme = str(options.get("shear_scheme", "mitc4")).strip().lower()
    if shear_scheme not in {"full", "reduced", "mitc4"}:
        raise PlateInputError("analysis.shear_scheme must be full, reduced, or mitc4")
    try:
        thinness_threshold = float(options.get("thinness_threshold", 1.0 / 20.0))
        system = assemble_plate_system(
            mesh=mesh,
            material=material,
            load=load,
            plate_method=plate_method,
            shear_scheme=shear_scheme,
            thinness_threshold=thinness_threshold,
        )
        solution = solve_dirichlet(system.stiffness, system.force, prescribed)
        recovered = recover_element_centers(
            mesh,
            solution.displacement,
            material,
            shear_scheme=shear_scheme,
            plate_method=system.selection.method,
            thinness_threshold=thinness_threshold,
        )
    except np.linalg.LinAlgError as exc:
        raise PlateInputError(
            "the constrained stiffness matrix is singular; check boundary conditions"
        ) from exc
    except (TypeError, ValueError) as exc:
        raise PlateInputError(str(exc)) from exc
    return PlateAnalysis(
        name=str(root.get("name", "plate-analysis")),
        mesh=mesh,
        material=material,
        system=system,
        solution=solution,
        recovered=recovered,
        shear_scheme=shear_scheme,
        boundary_type=boundary_type,
    )


def _equivalent_plane_stress(stress: np.ndarray) -> float:
    sx, sy, txy = stress
    return float(np.sqrt(max(sx * sx - sx * sy + sy * sy + 3.0 * txy * txy, 0.0)))


def serialize_analysis(analysis: PlateAnalysis) -> dict[str, Any]:
    """Convert NumPy-heavy analysis data into a JSON-compatible response."""

    mesh = analysis.mesh
    displacement = analysis.solution.displacement.reshape(-1, 3)
    reactions = analysis.solution.reactions.reshape(-1, 3)
    force = analysis.system.force.reshape(-1, 3)
    w = displacement[:, 0]
    center_point = mesh.nodes.mean(axis=0)
    center_node = int(np.argmin(np.linalg.norm(mesh.nodes - center_point, axis=1)))
    peak_node = int(np.argmax(np.abs(w)))
    free = analysis.solution.free_dofs
    free_residual = (
        float(np.linalg.norm(analysis.solution.reactions[free], ord=np.inf))
        if len(free)
        else 0.0
    )
    selection = analysis.system.selection
    summary = {
        "plate_method": selection.method,
        "theory": selection.theory,
        "characteristic_length": selection.characteristic_length,
        "thickness_ratio": selection.thickness_ratio,
        "thinness_threshold": selection.thinness_threshold,
        "node_count": int(len(mesh.nodes)),
        "element_count": int(len(mesh.elements)),
        "dof_count": int(mesh.ndof),
        "constrained_dof_count": int(len(analysis.solution.constrained_dofs)),
        "center_node": center_node,
        "center_deflection": float(w[center_node]),
        "peak_deflection_node": peak_node,
        "peak_deflection": float(w[peak_node]),
        "total_transverse_load": float(force[:, 0].sum()),
        "total_transverse_reaction": float(reactions[:, 0].sum()),
        "transverse_equilibrium_error": float(
            reactions[:, 0].sum() + force[:, 0].sum()
        ),
        "free_residual_inf_norm": free_residual,
    }
    nodal_results = [
        {
            "node": int(index),
            "x": float(point[0]),
            "y": float(point[1]),
            "w": float(values[0]),
            "theta_x": float(values[1]),
            "theta_y": float(values[2]),
            "reaction_w": float(reactions[index, 0]),
            "reaction_theta_x": float(reactions[index, 1]),
            "reaction_theta_y": float(reactions[index, 2]),
        }
        for index, (point, values) in enumerate(zip(mesh.nodes, displacement))
    ]
    element_results: list[dict[str, Any]] = []
    for item in analysis.recovered:
        response = item.response
        element_results.append(
            {
                "element": int(item.element_index),
                "center": response.physical_point.tolist(),
                "curvature": response.curvature.tolist(),
                "shear_strain": response.shear_strain.tolist(),
                "bending_moment": response.bending_moment.tolist(),
                "shear_force": response.shear_force.tolist(),
                "stress_top": response.stress_top.tolist(),
                "stress_bottom": response.stress_bottom.tolist(),
                "equivalent_stress_top": _equivalent_plane_stress(
                    response.stress_top
                ),
            }
        )
    return {
        "name": analysis.name,
        "summary": summary,
        "mesh": {
            "nodes": mesh.nodes.tolist(),
            "elements": mesh.elements.tolist(),
        },
        "nodal_results": nodal_results,
        "element_results": element_results,
    }
