"""Discrete Kirchhoff quadrilateral (DKQ) thin-plate kinematics.

The implementation keeps the package-wide nodal convention
``[w, theta_x, theta_y]`` with ``theta_x = w,x`` and
``theta_y = w,y``.  Kirchhoff constraints are imposed on every element edge:

* the tangent slope at the edge midpoint follows the cubic Hermite edge field;
* the normal slope at the edge midpoint is the average of the end values.

The four corner and four constrained midpoint slopes are interpolated with Q8
serendipity functions.  Differentiating that field produces the DKQ curvature
matrix without adding midside degrees of freedom.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .q4 import Q4Point, map_q4


FloatArray = NDArray[np.float64]


def q8_shape_functions(xi: float, eta: float) -> tuple[FloatArray, FloatArray]:
    """Return Q8 values and natural gradients in DKQ edge-node order.

    Nodes 1--4 are the Q4 corners and nodes 5--8 are the midpoints of edges
    1-2, 2-3, 3-4, and 4-1, respectively.
    """

    values = np.array(
        [
            -0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 + xi + eta),
            -0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - xi + eta),
            -0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - xi - eta),
            -0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 + xi - eta),
            0.5 * (1.0 - xi * xi) * (1.0 - eta),
            0.5 * (1.0 + xi) * (1.0 - eta * eta),
            0.5 * (1.0 - xi * xi) * (1.0 + eta),
            0.5 * (1.0 - xi) * (1.0 - eta * eta),
        ],
        dtype=float,
    )
    derivative_xi = np.array(
        [
            0.25 * (1.0 - eta) * (2.0 * xi + eta),
            0.25 * (1.0 - eta) * (2.0 * xi - eta),
            0.25 * (1.0 + eta) * (2.0 * xi + eta),
            0.25 * (1.0 + eta) * (2.0 * xi - eta),
            -xi * (1.0 - eta),
            0.5 * (1.0 - eta * eta),
            -xi * (1.0 + eta),
            -0.5 * (1.0 - eta * eta),
        ],
        dtype=float,
    )
    derivative_eta = np.array(
        [
            0.25 * (1.0 - xi) * (xi + 2.0 * eta),
            0.25 * (1.0 + xi) * (-xi + 2.0 * eta),
            0.25 * (1.0 + xi) * (xi + 2.0 * eta),
            0.25 * (1.0 - xi) * (-xi + 2.0 * eta),
            -0.5 * (1.0 - xi * xi),
            -(1.0 + xi) * eta,
            0.5 * (1.0 - xi * xi),
            -(1.0 - xi) * eta,
        ],
        dtype=float,
    )
    return values, np.column_stack((derivative_xi, derivative_eta))


def _edge_midpoint_gradient_operators(
    coordinates: FloatArray,
) -> tuple[FloatArray, ...]:
    """Express the four constrained midpoint gradients in corner DOFs."""

    operators: list[FloatArray] = []
    scale = max(float(np.ptp(coordinates, axis=0).max()), 1.0)
    tolerance = np.finfo(float).eps * scale
    for start, end in ((0, 1), (1, 2), (2, 3), (3, 0)):
        edge = coordinates[end] - coordinates[start]
        length = float(np.linalg.norm(edge))
        if length <= tolerance:
            raise ValueError("DKQ element contains a zero-length edge")
        tangent = edge / length
        normal = np.array([-tangent[1], tangent[0]], dtype=float)

        # g_s(mid) = 3(w_j-w_i)/(2L) - (g_s(i)+g_s(j))/4
        # g_n(mid) = (g_n(i)+g_n(j))/2
        # and grad(w) = tangent*g_s + normal*g_n.
        operator = np.zeros((2, 12), dtype=float)
        operator[:, 3 * start] -= 1.5 * tangent / length
        operator[:, 3 * end] += 1.5 * tangent / length
        for node in (start, end):
            operator[:, 3 * node + 1] += (
                -0.25 * tangent * tangent[0] + 0.5 * normal * normal[0]
            )
            operator[:, 3 * node + 2] += (
                -0.25 * tangent * tangent[1] + 0.5 * normal * normal[1]
            )
        operators.append(operator)
    return tuple(operators)


def dkq_kinematic_matrices(
    coordinates: ArrayLike, xi: float, eta: float
) -> tuple[Q4Point, FloatArray, FloatArray]:
    """Return mapped point, slope interpolation, and DKQ curvature matrix.

    The slope matrix maps the 12 corner DOFs to ``[w,x, w,y]``.  The curvature
    matrix maps them to ``[w,xx, w,yy, 2*w,xy]``.
    """

    coords = np.asarray(coordinates, dtype=float)
    if coords.shape != (4, 2):
        raise ValueError("DKQ coordinates must have shape (4,2)")
    point = map_q4(coords, xi, eta)
    values, natural_gradients = q8_shape_functions(xi, eta)

    slopes = np.zeros((2, 12), dtype=float)
    natural_slope_gradients = np.zeros((2, 12, 2), dtype=float)
    for node in range(4):
        slopes[0, 3 * node + 1] += values[node]
        slopes[1, 3 * node + 2] += values[node]
        natural_slope_gradients[0, 3 * node + 1] += natural_gradients[node]
        natural_slope_gradients[1, 3 * node + 2] += natural_gradients[node]

    operators = _edge_midpoint_gradient_operators(coords)
    for midpoint, operator in enumerate(operators, start=4):
        slopes += values[midpoint] * operator
        natural_slope_gradients += (
            operator[:, :, np.newaxis] * natural_gradients[midpoint]
        )

    physical_slope_gradients = np.empty_like(natural_slope_gradients)
    for component in range(2):
        physical_slope_gradients[component] = np.linalg.solve(
            point.jacobian.T, natural_slope_gradients[component].T
        ).T

    bending = np.vstack(
        (
            physical_slope_gradients[0, :, 0],
            physical_slope_gradients[1, :, 1],
            physical_slope_gradients[0, :, 1]
            + physical_slope_gradients[1, :, 0],
        )
    )
    return point, slopes, bending


def dkq_bending_B(
    coordinates: ArrayLike, xi: float, eta: float
) -> FloatArray:
    """Return the 3-by-12 DKQ curvature-displacement matrix."""

    return dkq_kinematic_matrices(coordinates, xi, eta)[2]
