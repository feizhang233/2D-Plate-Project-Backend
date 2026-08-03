"""Element-center response recovery."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
from numpy.typing import ArrayLike

from .assembly import element_dof_indices
from .element import ElementResponse, ShearScheme, element_response
from .material import MindlinMaterial
from .mesh import Mesh
from .theory import (
    DEFAULT_THINNESS_RATIO,
    PlateMethodRequest,
    normalize_plate_method,
    select_plate_method,
)


@dataclass(frozen=True)
class RecoveredElement:
    element_index: int
    response: ElementResponse


@dataclass(frozen=True)
class RecoveredPoint:
    element_index: int
    natural_point: tuple[float, float]
    response: ElementResponse


def split_nodal_dofs(displacement: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """Split [w,theta_x,theta_y] interleaving into w and theta arrays."""

    values = np.asarray(displacement, dtype=float)
    if values.ndim != 1 or values.size % 3:
        raise ValueError("displacement must be a flat array with 3 dofs per node")
    return values[0::3].copy(), values.reshape(-1, 3)[:, 1:3].copy()


def recover_element_points(
    mesh: Mesh,
    displacement: ArrayLike,
    material: MindlinMaterial,
    natural_points: list[tuple[float, float]],
    shear_scheme: ShearScheme = "full",
    plate_method: PlateMethodRequest = "M",
    thinness_threshold: float = DEFAULT_THINNESS_RATIO,
) -> list[RecoveredPoint]:
    """Recover curvature, resultants, and surface stresses at requested points."""

    values = np.asarray(displacement, dtype=float)
    if values.shape != (mesh.ndof,):
        raise ValueError("displacement size does not match mesh")
    method = (
        select_plate_method(mesh, material.thickness, thinness_threshold).method
        if str(plate_method).strip().lower() == "auto"
        else normalize_plate_method(plate_method)
    )
    recovered: list[RecoveredPoint] = []
    for index, connectivity in enumerate(mesh.elements):
        dofs = values[element_dof_indices(connectivity)]
        for xi, eta in natural_points:
            response = element_response(
                mesh.nodes[connectivity],
                dofs,
                material,
                xi,
                eta,
                shear_scheme,
                method,
            )
            recovered.append(RecoveredPoint(index, (xi, eta), response))
    return recovered


def recover_element_gauss_points(
    mesh: Mesh,
    displacement: ArrayLike,
    material: MindlinMaterial,
    shear_scheme: ShearScheme = "full",
    plate_method: PlateMethodRequest = "M",
    thinness_threshold: float = DEFAULT_THINNESS_RATIO,
) -> list[RecoveredPoint]:
    """Recover all fields at the four 2x2 Gauss points of every element."""

    value = 1.0 / sqrt(3.0)
    points = [(-value, -value), (value, -value), (value, value), (-value, value)]
    return recover_element_points(
        mesh,
        displacement,
        material,
        points,
        shear_scheme,
        plate_method,
        thinness_threshold,
    )


def recover_element_centers(
    mesh: Mesh,
    displacement: ArrayLike,
    material: MindlinMaterial,
    shear_scheme: ShearScheme = "full",
    plate_method: PlateMethodRequest = "M",
    thinness_threshold: float = DEFAULT_THINNESS_RATIO,
) -> list[RecoveredElement]:
    points = recover_element_points(
        mesh,
        displacement,
        material,
        [(0.0, 0.0)],
        shear_scheme,
        plate_method,
        thinness_threshold,
    )
    return [RecoveredElement(item.element_index, item.response) for item in points]
