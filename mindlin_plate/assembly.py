"""Dense educational assembly and Dirichlet solution routines."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .element import (
    Load,
    MomentLoad,
    ShearScheme,
    plate_element_matrices,
    q4_consistent_load,
    q4_edge_consistent_load,
)
from .material import MindlinMaterial
from .mesh import BoundaryEdge, Mesh
from .theory import (
    DEFAULT_THINNESS_RATIO,
    PlateMethodRequest,
    PlateMethodSelection,
    normalize_plate_method,
    select_plate_method,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class LinearSolution:
    displacement: FloatArray
    reactions: FloatArray
    free_dofs: IntArray
    constrained_dofs: IntArray


@dataclass(frozen=True)
class AssembledPlateSystem:
    """An assembled system together with the K/M decision used to build it."""

    stiffness: FloatArray
    force: FloatArray
    selection: PlateMethodSelection


def element_dof_indices(connectivity: ArrayLike) -> IntArray:
    nodes = np.asarray(connectivity, dtype=np.int64)
    return np.array([3 * node + local for node in nodes for local in range(3)], dtype=np.int64)


def assemble_system(
    mesh: Mesh,
    material: MindlinMaterial,
    load: Load,
    shear_scheme: ShearScheme = "full",
    plate_method: PlateMethodRequest = "M",
    thinness_threshold: float = DEFAULT_THINNESS_RATIO,
) -> tuple[FloatArray, FloatArray]:
    """Assemble an explicit or automatically selected plate formulation.

    The default remains the original M formulation for backward compatibility.
    Pass ``plate_method="auto"`` to select K/M from the complete plate shape.
    """

    if str(plate_method).strip().lower() == "auto":
        method = select_plate_method(
            mesh, material.thickness, thinness_threshold
        ).method
    else:
        method = normalize_plate_method(plate_method)
    stiffness = np.zeros((mesh.ndof, mesh.ndof), dtype=float)
    force = np.zeros(mesh.ndof, dtype=float)
    for connectivity in mesh.elements:
        coordinates = mesh.nodes[connectivity]
        dofs = element_dof_indices(connectivity)
        element = plate_element_matrices(
            coordinates, material, method, shear_scheme
        )
        element_force = q4_consistent_load(coordinates, load)
        stiffness[np.ix_(dofs, dofs)] += element.total
        force[dofs] += element_force
    return stiffness, force


def assemble_plate_system(
    mesh: Mesh,
    material: MindlinMaterial,
    load: Load,
    plate_method: PlateMethodRequest = "auto",
    shear_scheme: ShearScheme = "mitc4",
    thinness_threshold: float = DEFAULT_THINNESS_RATIO,
) -> AssembledPlateSystem:
    """High-level assembly with automatic K/M selection enabled by default."""

    classification = select_plate_method(
        mesh, material.thickness, thinness_threshold
    )
    if str(plate_method).strip().lower() == "auto":
        selection = classification
    else:
        selection = replace(
            classification, method=normalize_plate_method(plate_method)
        )
    stiffness, force = assemble_system(
        mesh,
        material,
        load,
        shear_scheme,
        selection.method,
        thinness_threshold,
    )
    return AssembledPlateSystem(stiffness, force, selection)


def assemble_boundary_load(
    mesh: Mesh,
    edges: list[BoundaryEdge],
    transverse_shear: Load = 0.0,
    moment: MomentLoad = (0.0, 0.0),
) -> FloatArray:
    """Assemble natural shear and moment data on selected exterior edges."""

    force = np.zeros(mesh.ndof, dtype=float)
    for edge in edges:
        if edge.element_index < 0 or edge.element_index >= len(mesh.elements):
            raise ValueError("boundary edge has an invalid element index")
        connectivity = mesh.elements[edge.element_index]
        dofs = element_dof_indices(connectivity)
        force[dofs] += q4_edge_consistent_load(
            mesh.nodes[connectivity], edge.local_edge, transverse_shear, moment
        )
    return force


def solve_dirichlet(
    stiffness: ArrayLike, force: ArrayLike, prescribed: Mapping[int, float]
) -> LinearSolution:
    matrix = np.asarray(stiffness, dtype=float)
    rhs = np.asarray(force, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("stiffness must be square")
    if rhs.shape != (matrix.shape[0],):
        raise ValueError("force size does not match stiffness")
    constrained = np.array(sorted(prescribed), dtype=np.int64)
    if constrained.size and (constrained.min() < 0 or constrained.max() >= len(rhs)):
        raise ValueError("prescribed dof index is out of range")
    all_dofs = np.arange(len(rhs), dtype=np.int64)
    free = np.setdiff1d(all_dofs, constrained, assume_unique=True)
    displacement = np.zeros(len(rhs), dtype=float)
    if constrained.size:
        displacement[constrained] = [prescribed[int(i)] for i in constrained]
    reduced_rhs = rhs[free]
    if constrained.size:
        reduced_rhs = reduced_rhs - matrix[np.ix_(free, constrained)] @ displacement[constrained]
    displacement[free] = np.linalg.solve(matrix[np.ix_(free, free)], reduced_rhs)
    reactions = matrix @ displacement - rhs
    return LinearSolution(displacement, reactions, free, constrained)


def _bounds(mesh: Mesh) -> tuple[float, float, float, float, float]:
    x_min, y_min = mesh.nodes.min(axis=0)
    x_max, y_max = mesh.nodes.max(axis=0)
    tol = 1e-10 * max(x_max - x_min, y_max - y_min, 1.0)
    return x_min, x_max, y_min, y_max, tol


def _on(value: float, target: float, tolerance: float) -> bool:
    return abs(value - target) <= tolerance


def clamped_dofs(mesh: Mesh) -> dict[int, float]:
    x_min, x_max, y_min, y_max, tol = _bounds(mesh)
    result: dict[int, float] = {}
    for node, (x, y) in enumerate(mesh.nodes):
        if _on(x, x_min, tol) or _on(x, x_max, tol) or _on(y, y_min, tol) or _on(y, y_max, tol):
            for local in range(3):
                result[3 * node + local] = 0.0
    return result


def soft_simply_supported_dofs(mesh: Mesh) -> dict[int, float]:
    x_min, x_max, y_min, y_max, tol = _bounds(mesh)
    result: dict[int, float] = {}
    for node, (x, y) in enumerate(mesh.nodes):
        if _on(x, x_min, tol) or _on(x, x_max, tol) or _on(y, y_min, tol) or _on(y, y_max, tol):
            result[3 * node] = 0.0
    return result


def hard_simply_supported_dofs(mesh: Mesh) -> dict[int, float]:
    """Set w=0 and tangential rotation=0 on axis-aligned boundaries."""

    x_min, x_max, y_min, y_max, tol = _bounds(mesh)
    result: dict[int, float] = {}
    for node, (x, y) in enumerate(mesh.nodes):
        on_x = _on(x, x_min, tol) or _on(x, x_max, tol)
        on_y = _on(y, y_min, tol) or _on(y, y_max, tol)
        if on_x or on_y:
            result[3 * node] = 0.0
        if on_x:  # edge tangent is y
            result[3 * node + 2] = 0.0
        if on_y:  # edge tangent is x
            result[3 * node + 1] = 0.0
    return result


def symmetry_dofs(mesh: Mesh, boundary: str) -> dict[int, float]:
    """Constrain theta_n=0 on one axis-aligned boundary."""

    x_min, x_max, y_min, y_max, tol = _bounds(mesh)
    definitions = {
        "x_min": (0, x_min, 1),
        "x_max": (0, x_max, 1),
        "y_min": (1, y_min, 2),
        "y_max": (1, y_max, 2),
    }
    if boundary not in definitions:
        raise ValueError("boundary must be x_min, x_max, y_min, or y_max")
    coordinate, target, rotation_dof = definitions[boundary]
    return {
        3 * node + rotation_dof: 0.0
        for node, point in enumerate(mesh.nodes)
        if _on(float(point[coordinate]), target, tol)
    }
