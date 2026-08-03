"""Analytical benchmark from chapter 11 of the reference."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sin
from typing import Callable

from .material import MindlinMaterial


@dataclass(frozen=True)
class SinusoidalSolution:
    center_deflection: float
    rotation_amplitude: float
    bending_deflection: float
    shear_deflection: float


def sinusoidal_load(q0: float, length_x: float, length_y: float) -> Callable[[float, float], float]:
    if length_x <= 0.0 or length_y <= 0.0:
        raise ValueError("plate lengths must be positive")

    def load(x: float, y: float) -> float:
        return q0 * sin(pi * x / length_x) * sin(pi * y / length_y)

    return load


def simply_supported_sinusoidal_solution(
    side: float, material: MindlinMaterial, q0: float
) -> SinusoidalSolution:
    """Single-mode solution for the square plate described on pages 36-42."""

    if side <= 0.0:
        raise ValueError("side must be positive")
    p = pi / side
    bending = q0 / (4.0 * material.bending_rigidity * p**4)
    shear = q0 / (2.0 * material.shear_rigidity * p**2)
    total = bending + shear
    rotation = material.shear_rigidity * p / (
        material.shear_rigidity + 2.0 * material.bending_rigidity * p**2
    ) * total
    return SinusoidalSolution(total, rotation, bending, shear)

