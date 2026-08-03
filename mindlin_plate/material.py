"""Constitutive matrices from chapters 3-4 of the reference."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class MindlinMaterial:
    """Homogeneous isotropic Reissner-Mindlin plate material.

    Parameters use any consistent unit system. ``shear_correction`` defaults to
    the rectangular-section energy-equivalent value 5/6.
    """

    young: float
    poisson: float
    thickness: float
    shear_correction: float = 5.0 / 6.0

    def __post_init__(self) -> None:
        if self.young <= 0.0:
            raise ValueError("young must be positive")
        if not (-1.0 < self.poisson < 0.5):
            raise ValueError("poisson must satisfy -1 < nu < 0.5")
        if self.thickness <= 0.0:
            raise ValueError("thickness must be positive")
        if self.shear_correction <= 0.0:
            raise ValueError("shear_correction must be positive")

    @property
    def shear_modulus(self) -> float:
        return self.young / (2.0 * (1.0 + self.poisson))

    @property
    def plane_stress_matrix(self) -> FloatArray:
        nu = self.poisson
        return self.young / (1.0 - nu * nu) * np.array(
            [[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, (1.0 - nu) / 2.0]],
            dtype=float,
        )

    @property
    def bending_rigidity(self) -> float:
        t = self.thickness
        return self.young * t**3 / (12.0 * (1.0 - self.poisson**2))

    @property
    def bending_matrix(self) -> FloatArray:
        return (self.thickness**3 / 12.0) * self.plane_stress_matrix

    @property
    def shear_rigidity(self) -> float:
        return self.shear_correction * self.shear_modulus * self.thickness

    @property
    def shear_matrix(self) -> FloatArray:
        return self.shear_rigidity * np.eye(2, dtype=float)

    def bending_stress(self, curvature: ArrayLike, z: float) -> FloatArray:
        """Return [sigma_x, sigma_y, tau_xy] at thickness coordinate ``z``."""

        if abs(z) > self.thickness / 2.0 + 1e-14:
            raise ValueError("z must lie inside the plate thickness")
        kappa = np.asarray(curvature, dtype=float)
        if kappa.shape != (3,):
            raise ValueError("curvature must have shape (3,)")
        return -z * (self.plane_stress_matrix @ kappa)

    def parabolic_shear_stress(self, shear_force: ArrayLike, z: float) -> FloatArray:
        """Recover the reference's parabolic shear stress from [Qx, Qy]."""

        if abs(z) > self.thickness / 2.0 + 1e-14:
            raise ValueError("z must lie inside the plate thickness")
        force = np.asarray(shear_force, dtype=float)
        if force.shape != (2,):
            raise ValueError("shear_force must have shape (2,)")
        t = self.thickness
        return 3.0 / (2.0 * t) * (1.0 - 4.0 * z * z / (t * t)) * force

