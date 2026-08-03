"""Continuum Reissner-Mindlin kinematics from chapter 2."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def displacement_at_z(w: float, theta: ArrayLike, z: float) -> FloatArray:
    """Return the 3-D displacement [-z*theta_x, -z*theta_y, w]."""

    rotation = np.asarray(theta, dtype=float)
    if rotation.shape != (2,):
        raise ValueError("theta must have shape (2,)")
    return np.array([-z * rotation[0], -z * rotation[1], w], dtype=float)


def generalized_strains(
    grad_w: ArrayLike, theta: ArrayLike, grad_theta: ArrayLike
) -> tuple[FloatArray, FloatArray]:
    """Return curvature and transverse engineering shear strain.

    ``grad_theta`` is arranged as::

        [[theta_x,x, theta_x,y],
         [theta_y,x, theta_y,y]]
    """

    dw = np.asarray(grad_w, dtype=float)
    rotation = np.asarray(theta, dtype=float)
    drotation = np.asarray(grad_theta, dtype=float)
    if dw.shape != (2,) or rotation.shape != (2,) or drotation.shape != (2, 2):
        raise ValueError("expected grad_w (2,), theta (2,), grad_theta (2,2)")
    curvature = np.array(
        [drotation[0, 0], drotation[1, 1], drotation[0, 1] + drotation[1, 0]],
        dtype=float,
    )
    shear = dw - rotation
    return curvature, shear


def small_strains_at_z(
    grad_w: ArrayLike, theta: ArrayLike, grad_theta: ArrayLike, z: float
) -> FloatArray:
    """Return [eps_x, eps_y, gamma_xy, gamma_xz, gamma_yz]."""

    curvature, shear = generalized_strains(grad_w, theta, grad_theta)
    return np.concatenate((-z * curvature, shear))

