"""Q4 interpolation, mapping, raw strain matrices, and MITC4 tying."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Q4Point:
    natural: FloatArray
    physical: FloatArray
    shape: FloatArray
    natural_gradients: FloatArray
    physical_gradients: FloatArray
    jacobian: FloatArray
    det_jacobian: float


def _coordinates(coordinates: ArrayLike) -> FloatArray:
    coords = np.asarray(coordinates, dtype=float)
    if coords.shape != (4, 2):
        raise ValueError("Q4 coordinates must have shape (4,2)")
    return coords


def shape_functions(xi: float, eta: float) -> tuple[FloatArray, FloatArray]:
    """Return Q4 shape values and gradients with columns [d/dxi, d/deta]."""

    signs = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    sx = signs[:, 0]
    se = signs[:, 1]
    values = 0.25 * (1.0 + sx * xi) * (1.0 + se * eta)
    gradients = np.column_stack(
        (0.25 * sx * (1.0 + se * eta), 0.25 * se * (1.0 + sx * xi))
    )
    return values, gradients


def map_q4(coordinates: ArrayLike, xi: float, eta: float) -> Q4Point:
    """Map one natural-coordinate point and validate its Jacobian."""

    coords = _coordinates(coordinates)
    values, natural_gradients = shape_functions(xi, eta)
    # J = [[x,xi, x,eta], [y,xi, y,eta]], as used in the reference.
    jacobian = coords.T @ natural_gradients
    det_jacobian = float(np.linalg.det(jacobian))
    scale = max(float(np.linalg.norm(jacobian, ord=np.inf)), 1.0)
    if det_jacobian <= np.finfo(float).eps * scale * scale:
        raise ValueError(f"Q4 element has non-positive/degenerate det(J)={det_jacobian:.6e}")
    # [N,xi, N,eta]^T = J^T [N,x, N,y]^T, hence grad_x N = J^-T grad_xi N.
    physical_gradients = np.linalg.solve(jacobian.T, natural_gradients.T).T
    return Q4Point(
        natural=np.array([xi, eta], dtype=float),
        physical=values @ coords,
        shape=values,
        natural_gradients=natural_gradients,
        physical_gradients=physical_gradients,
        jacobian=jacobian,
        det_jacobian=det_jacobian,
    )


def bending_B(physical_gradients: ArrayLike) -> FloatArray:
    gradients = np.asarray(physical_gradients, dtype=float)
    if gradients.shape != (4, 2):
        raise ValueError("physical_gradients must have shape (4,2)")
    matrix = np.zeros((3, 12), dtype=float)
    for i, (dn_dx, dn_dy) in enumerate(gradients):
        col = 3 * i
        matrix[0, col + 1] = dn_dx
        matrix[1, col + 2] = dn_dy
        matrix[2, col + 1] = dn_dy
        matrix[2, col + 2] = dn_dx
    return matrix


def shear_B(shape: ArrayLike, physical_gradients: ArrayLike) -> FloatArray:
    values = np.asarray(shape, dtype=float)
    gradients = np.asarray(physical_gradients, dtype=float)
    if values.shape != (4,) or gradients.shape != (4, 2):
        raise ValueError("expected shape (4,) and physical_gradients (4,2)")
    matrix = np.zeros((2, 12), dtype=float)
    for i, ((dn_dx, dn_dy), value) in enumerate(zip(gradients, values)):
        col = 3 * i
        matrix[0, col : col + 3] = [dn_dx, -value, 0.0]
        matrix[1, col : col + 3] = [dn_dy, 0.0, -value]
    return matrix


def kinematic_matrices(
    coordinates: ArrayLike, xi: float, eta: float
) -> tuple[Q4Point, FloatArray, FloatArray]:
    point = map_q4(coordinates, xi, eta)
    return point, bending_B(point.physical_gradients), shear_B(
        point.shape, point.physical_gradients
    )


def _covariant_raw_shear_B(
    coordinates: FloatArray, xi: float, eta: float
) -> FloatArray:
    point, _, raw = kinematic_matrices(coordinates, xi, eta)
    return point.jacobian.T @ raw


def mitc4_shear_B(coordinates: ArrayLike, xi: float, eta: float) -> FloatArray:
    """Return the MITC4 assumed physical shear-strain matrix.

    A=(0,1), C=(0,-1) tie gamma_xi. B=(-1,0), D=(1,0) tie
    gamma_eta. The tied covariant field is transformed back with the actual
    Jacobian at the requested point.
    """

    coords = _coordinates(coordinates)
    cov_a = _covariant_raw_shear_B(coords, 0.0, 1.0)
    cov_c = _covariant_raw_shear_B(coords, 0.0, -1.0)
    cov_b = _covariant_raw_shear_B(coords, -1.0, 0.0)
    cov_d = _covariant_raw_shear_B(coords, 1.0, 0.0)

    reduced_covariant = np.empty((2, 12), dtype=float)
    reduced_covariant[0] = 0.5 * (1.0 - eta) * cov_c[0] + 0.5 * (
        1.0 + eta
    ) * cov_a[0]
    reduced_covariant[1] = 0.5 * (1.0 - xi) * cov_b[1] + 0.5 * (
        1.0 + xi
    ) * cov_d[1]

    point = map_q4(coords, xi, eta)
    return np.linalg.solve(point.jacobian.T, reduced_covariant)
