"""Natural loads and local-coordinate constraints for plate boundaries."""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, cos, sin
from typing import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .assembly import solve_dirichlet


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class LocalConstraintSolution:
    """Solution of constraints expressed in node-local coordinates.

    ``displacement`` and ``reactions`` are returned in global coordinates;
    the corresponding local arrays and the transformation are retained for
    checking or applying additional local results.
    """

    displacement: FloatArray
    reactions: FloatArray
    local_displacement: FloatArray
    local_reactions: FloatArray
    transformation: FloatArray
    free_dofs: NDArray[np.int64]
    constrained_dofs: NDArray[np.int64]


def edge_angle(point_a: ArrayLike, point_b: ArrayLike) -> float:
    """Return the global angle of the local edge tangent x' axis."""

    a = np.asarray(point_a, dtype=float)
    b = np.asarray(point_b, dtype=float)
    if a.shape != (2,) or b.shape != (2,):
        raise ValueError("edge points must have shape (2,)")
    delta = b - a
    if np.linalg.norm(delta) <= np.finfo(float).eps:
        raise ValueError("edge points must be distinct")
    return atan2(float(delta[1]), float(delta[0]))


def nodal_rotation_transform(angle: float) -> FloatArray:
    """Return L such that [w,theta_x,theta_y] = L [w,theta_x',theta_y']."""

    c = cos(angle)
    s = sin(angle)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)


def global_dof_transform(
    node_count: int, local_angles: Mapping[int, float]
) -> FloatArray:
    """Build the block transformation for nodes using local boundary axes."""

    if node_count < 1:
        raise ValueError("node_count must be positive")
    transform = np.eye(3 * node_count, dtype=float)
    for node, angle in local_angles.items():
        if node < 0 or node >= node_count:
            raise ValueError(f"local angle node {node} is out of range")
        dofs = slice(3 * node, 3 * node + 3)
        transform[dofs, dofs] = nodal_rotation_transform(float(angle))
    return transform


def transform_system_to_local(
    stiffness: ArrayLike, force: ArrayLike, transformation: ArrayLike
) -> tuple[FloatArray, FloatArray]:
    """Apply K'=L^T K L and f'=L^T f."""

    matrix = np.asarray(stiffness, dtype=float)
    rhs = np.asarray(force, dtype=float)
    transform = np.asarray(transformation, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("stiffness must be square")
    if transform.shape != matrix.shape or rhs.shape != (matrix.shape[0],):
        raise ValueError("system and transformation sizes do not match")
    return transform.T @ matrix @ transform, transform.T @ rhs


def local_edge_constraints(
    node_ids: ArrayLike, kind: str
) -> dict[int, float]:
    """Return constraints in an edge basis whose x' axis is tangential.

    Supported kinds:
    - ``clamped``: w=theta_x'=theta_y'=0
    - ``hard_simply_supported``: w=theta_x'=0
    - ``soft_simply_supported``: w=0
    - ``symmetry``: theta_y'=theta_n=0
    """

    nodes = np.asarray(node_ids, dtype=np.int64)
    local_dofs = {
        "clamped": (0, 1, 2),
        "hard_simply_supported": (0, 1),
        "soft_simply_supported": (0,),
        "symmetry": (2,),
    }
    if kind not in local_dofs:
        raise ValueError(f"unknown local edge constraint kind: {kind}")
    result: dict[int, float] = {}
    for node in nodes:
        if node < 0:
            raise ValueError("node indices must be non-negative")
        for local in local_dofs[kind]:
            result[3 * int(node) + local] = 0.0
    return result


def solve_with_local_constraints(
    stiffness: ArrayLike,
    force: ArrayLike,
    node_count: int,
    local_angles: Mapping[int, float],
    prescribed_local: Mapping[int, float],
) -> LocalConstraintSolution:
    """Transform, constrain, solve, and return both global and local results."""

    matrix = np.asarray(stiffness, dtype=float)
    rhs = np.asarray(force, dtype=float)
    transform = global_dof_transform(node_count, local_angles)
    local_matrix, local_force = transform_system_to_local(matrix, rhs, transform)
    local_solution = solve_dirichlet(local_matrix, local_force, prescribed_local)
    displacement = transform @ local_solution.displacement
    reactions = matrix @ displacement - rhs
    local_reactions = transform.T @ reactions
    return LocalConstraintSolution(
        displacement=displacement,
        reactions=reactions,
        local_displacement=local_solution.displacement,
        local_reactions=local_reactions,
        transformation=transform,
        free_dofs=local_solution.free_dofs,
        constrained_dofs=local_solution.constrained_dofs,
    )

