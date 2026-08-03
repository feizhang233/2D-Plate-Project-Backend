"""Geometry-based selection between Kirchhoff (K) and Mindlin (M) plates."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]
PlateMethod = Literal["K", "M"]
PlateMethodRequest = Literal[
    "K", "M", "k", "m", "kirchhoff", "mindlin", "dkq", "auto"
]
DEFAULT_THINNESS_RATIO = 1.0 / 20.0


@dataclass(frozen=True)
class PlateMethodSelection:
    """Result of classifying a plate from its planform and thickness."""

    method: PlateMethod
    characteristic_length: float
    thickness: float
    thickness_ratio: float
    thinness_threshold: float

    @property
    def theory(self) -> str:
        return "kirchhoff" if self.method == "K" else "mindlin"

    @property
    def is_thin(self) -> bool:
        return self.thickness_ratio <= self.thinness_threshold


def _points(shape: ArrayLike | object) -> FloatArray:
    nodes = getattr(shape, "nodes", shape)
    points = np.asarray(nodes, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("plate shape must provide planar nodes with shape (n,2)")
    points = np.unique(points, axis=0)
    if len(points) < 3:
        raise ValueError("plate shape must contain at least three distinct points")
    return points


def _convex_hull(points: FloatArray) -> FloatArray:
    """Return the CCW convex hull using Andrew's monotone-chain algorithm."""

    ordered = sorted(map(tuple, points.tolist()))

    def cross(
        origin: tuple[float, float],
        first: tuple[float, float],
        second: tuple[float, float],
    ) -> float:
        return (first[0] - origin[0]) * (second[1] - origin[1]) - (
            first[1] - origin[1]
        ) * (second[0] - origin[0])

    lower: list[tuple[float, float]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0.0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0.0:
            upper.pop()
        upper.append(point)
    hull = np.asarray(lower[:-1] + upper[:-1], dtype=float)
    if len(hull) < 3:
        raise ValueError("plate shape is collinear or has zero planform area")
    return hull


def plate_characteristic_length(shape: ArrayLike | object) -> float:
    """Return the rotation-invariant minimum planform span.

    The minimum width of the convex hull is the governing short span for a
    rectangle and remains well-defined for rotated, skew, and polygonal plates.
    It describes the whole plate, so mesh refinement does not change the K/M
    classification.
    """

    hull = _convex_hull(_points(shape))
    widths: list[float] = []
    for start, end in zip(hull, np.roll(hull, -1, axis=0)):
        edge = end - start
        length = float(np.linalg.norm(edge))
        if length == 0.0:
            continue
        normal = np.array([-edge[1], edge[0]], dtype=float) / length
        projection = hull @ normal
        widths.append(float(projection.max() - projection.min()))
    characteristic_length = min(widths, default=0.0)
    scale = max(float(np.ptp(hull, axis=0).max()), 1.0)
    if characteristic_length <= np.finfo(float).eps * scale:
        raise ValueError("plate shape has zero characteristic span")
    return characteristic_length


def normalize_plate_method(method: str) -> PlateMethod:
    """Normalize user-facing K/M names to the compact method code."""

    normalized = method.strip().lower()
    if normalized in {"k", "kirchhoff", "dkq"}:
        return "K"
    if normalized in {"m", "mindlin"}:
        return "M"
    raise ValueError("plate method must be K/kirchhoff/dkq or M/mindlin")


def select_plate_method(
    shape: ArrayLike | object,
    thickness: float,
    thinness_threshold: float = DEFAULT_THINNESS_RATIO,
) -> PlateMethodSelection:
    """Choose K for a thin plate and M when transverse shear is significant.

    ``thickness / characteristic_length <= thinness_threshold`` selects the
    discrete Kirchhoff method.  The default threshold is 1/20.
    """

    if thickness <= 0.0:
        raise ValueError("thickness must be positive")
    if not (0.0 < thinness_threshold < 1.0):
        raise ValueError("thinness_threshold must lie between zero and one")
    characteristic_length = plate_characteristic_length(shape)
    ratio = float(thickness / characteristic_length)
    method: PlateMethod = "K" if ratio <= thinness_threshold else "M"
    return PlateMethodSelection(
        method=method,
        characteristic_length=characteristic_length,
        thickness=float(thickness),
        thickness_ratio=ratio,
        thinness_threshold=float(thinness_threshold),
    )


def choose_plate_method(
    shape: ArrayLike | object,
    thickness: float,
    thinness_threshold: float = DEFAULT_THINNESS_RATIO,
) -> PlateMethod:
    """Return only the selected ``"K"`` or ``"M"`` method code."""

    return select_plate_method(shape, thickness, thinness_threshold).method
