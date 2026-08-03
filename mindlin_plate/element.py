"""Q4 element stiffness, loading, energies, and point response."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .kirchhoff import dkq_kinematic_matrices
from .material import MindlinMaterial
from .q4 import kinematic_matrices, map_q4, mitc4_shear_B
from .theory import PlateMethodRequest, normalize_plate_method


FloatArray = NDArray[np.float64]
ShearScheme = Literal["full", "reduced", "mitc4"]
Load = float | Callable[[float, float], float]
MomentLoad = ArrayLike | Callable[[float, float], ArrayLike]


@dataclass(frozen=True)
class ElementMatrices:
    bending: FloatArray
    shear: FloatArray

    @property
    def total(self) -> FloatArray:
        return self.bending + self.shear


@dataclass(frozen=True)
class ElementResponse:
    physical_point: FloatArray
    curvature: FloatArray
    shear_strain: FloatArray
    bending_moment: FloatArray
    shear_force: FloatArray
    stress_top: FloatArray
    stress_bottom: FloatArray
    plate_method: str = "M"


def gauss_rule(order: int) -> tuple[FloatArray, FloatArray]:
    if order == 1:
        return np.array([[0.0, 0.0]], dtype=float), np.array([4.0], dtype=float)
    if order == 2:
        value = 1.0 / np.sqrt(3.0)
        points = np.array(
            [[-value, -value], [value, -value], [value, value], [-value, value]],
            dtype=float,
        )
        return points, np.ones(4, dtype=float)
    raise ValueError("only one- and two-point tensor-product rules are supported")


def q4_element_matrices(
    coordinates: ArrayLike,
    material: MindlinMaterial,
    shear_scheme: ShearScheme = "full",
) -> ElementMatrices:
    """Integrate Q4 bending and selected shear stiffness."""

    if shear_scheme not in {"full", "reduced", "mitc4"}:
        raise ValueError("shear_scheme must be 'full', 'reduced', or 'mitc4'")
    coords = np.asarray(coordinates, dtype=float)
    bending = np.zeros((12, 12), dtype=float)
    shear = np.zeros((12, 12), dtype=float)

    points, weights = gauss_rule(2)
    for (xi, eta), weight in zip(points, weights):
        point, b_bending, _ = kinematic_matrices(coords, float(xi), float(eta))
        bending += (
            b_bending.T @ material.bending_matrix @ b_bending
        ) * point.det_jacobian * weight

    shear_order = 1 if shear_scheme == "reduced" else 2
    points, weights = gauss_rule(shear_order)
    for (xi, eta), weight in zip(points, weights):
        point, _, raw_shear = kinematic_matrices(coords, float(xi), float(eta))
        b_shear = (
            mitc4_shear_B(coords, float(xi), float(eta))
            if shear_scheme == "mitc4"
            else raw_shear
        )
        shear += (
            b_shear.T @ material.shear_matrix @ b_shear
        ) * point.det_jacobian * weight

    # Suppress only round-off asymmetry, not a modeling error.
    bending = 0.5 * (bending + bending.T)
    shear = 0.5 * (shear + shear.T)
    return ElementMatrices(bending=bending, shear=shear)


def dkq_element_matrices(
    coordinates: ArrayLike,
    material: MindlinMaterial,
) -> ElementMatrices:
    """Integrate the 12-DOF discrete Kirchhoff quadrilateral stiffness.

    Kirchhoff thin-plate theory contains bending energy only; the returned
    shear matrix is therefore exactly zero so that the result remains
    compatible with :class:`ElementMatrices` and the existing assembler.
    """

    coords = np.asarray(coordinates, dtype=float)
    bending = np.zeros((12, 12), dtype=float)
    points, weights = gauss_rule(2)
    for (xi, eta), weight in zip(points, weights):
        point, _, b_bending = dkq_kinematic_matrices(
            coords, float(xi), float(eta)
        )
        bending += (
            b_bending.T @ material.bending_matrix @ b_bending
        ) * point.det_jacobian * weight
    bending = 0.5 * (bending + bending.T)
    return ElementMatrices(bending=bending, shear=np.zeros((12, 12), dtype=float))


def plate_element_matrices(
    coordinates: ArrayLike,
    material: MindlinMaterial,
    plate_method: PlateMethodRequest = "M",
    shear_scheme: ShearScheme = "mitc4",
) -> ElementMatrices:
    """Build one explicit K- or M-method plate element."""

    if str(plate_method).strip().lower() == "auto":
        raise ValueError("auto plate selection requires the full mesh planform")
    method = normalize_plate_method(plate_method)
    if method == "K":
        return dkq_element_matrices(coordinates, material)
    return q4_element_matrices(coordinates, material, shear_scheme)


def q4_consistent_load(coordinates: ArrayLike, load: Load) -> FloatArray:
    """Integrate a transverse distributed load into 12 nodal forces."""

    coords = np.asarray(coordinates, dtype=float)
    result = np.zeros(12, dtype=float)
    points, weights = gauss_rule(2)
    for (xi, eta), weight in zip(points, weights):
        point = map_q4(coords, float(xi), float(eta))
        q = float(load(*point.physical)) if callable(load) else float(load)
        for i, value in enumerate(point.shape):
            result[3 * i] += value * q * point.det_jacobian * weight
    return result


def _edge_natural_point(edge: int, coordinate: float) -> tuple[float, float, int]:
    """Return xi, eta and the Jacobian tangent-column index for one Q4 edge."""

    if edge == 0:
        return coordinate, -1.0, 0
    if edge == 1:
        return 1.0, coordinate, 1
    if edge == 2:
        return coordinate, 1.0, 0
    if edge == 3:
        return -1.0, coordinate, 1
    raise ValueError("Q4 local edge must be 0, 1, 2, or 3")


def q4_edge_consistent_load(
    coordinates: ArrayLike,
    edge: int,
    transverse_shear: Load = 0.0,
    moment: MomentLoad = (0.0, 0.0),
) -> FloatArray:
    """Integrate natural edge data [V_bar, m_x_bar, m_y_bar].

    Edge numbering follows the CCW node pairs (1-2, 2-3, 3-4, 4-1),
    represented by zero-based local edge indices 0-3.
    """

    coords = np.asarray(coordinates, dtype=float)
    result = np.zeros(12, dtype=float)
    value = 1.0 / np.sqrt(3.0)
    for coordinate in (-value, value):
        xi, eta, tangent_column = _edge_natural_point(edge, coordinate)
        point = map_q4(coords, xi, eta)
        line_jacobian = float(np.linalg.norm(point.jacobian[:, tangent_column]))
        shear_value = (
            float(transverse_shear(*point.physical))
            if callable(transverse_shear)
            else float(transverse_shear)
        )
        moment_value = (
            np.asarray(moment(*point.physical), dtype=float)
            if callable(moment)
            else np.asarray(moment, dtype=float)
        )
        if moment_value.shape != (2,):
            raise ValueError("edge moment must have shape (2,)")
        traction = np.array([shear_value, moment_value[0], moment_value[1]])
        for i, shape_value in enumerate(point.shape):
            result[3 * i : 3 * i + 3] += shape_value * traction * line_jacobian
    return result


def element_response(
    coordinates: ArrayLike,
    element_dofs: ArrayLike,
    material: MindlinMaterial,
    xi: float = 0.0,
    eta: float = 0.0,
    shear_scheme: ShearScheme = "full",
    plate_method: PlateMethodRequest = "M",
) -> ElementResponse:
    coords = np.asarray(coordinates, dtype=float)
    dofs = np.asarray(element_dofs, dtype=float)
    if dofs.shape != (12,):
        raise ValueError("element_dofs must have shape (12,)")
    if str(plate_method).strip().lower() == "auto":
        raise ValueError("auto plate selection requires the full mesh planform")
    method = normalize_plate_method(plate_method)
    if method == "K":
        point, _, b_bending = dkq_kinematic_matrices(coords, xi, eta)
        b_shear = np.zeros((2, 12), dtype=float)
    else:
        point, b_bending, raw_shear = kinematic_matrices(coords, xi, eta)
        b_shear = (
            mitc4_shear_B(coords, xi, eta)
            if shear_scheme == "mitc4"
            else raw_shear
        )
    curvature = b_bending @ dofs
    shear_strain = b_shear @ dofs
    moment = material.bending_matrix @ curvature
    shear_force = material.shear_matrix @ shear_strain
    half_t = material.thickness / 2.0
    return ElementResponse(
        physical_point=point.physical,
        curvature=curvature,
        shear_strain=shear_strain,
        bending_moment=moment,
        shear_force=shear_force,
        stress_top=material.bending_stress(curvature, half_t),
        stress_bottom=material.bending_stress(curvature, -half_t),
        plate_method=method,
    )
