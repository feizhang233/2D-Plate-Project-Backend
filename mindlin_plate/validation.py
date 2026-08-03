"""Verification gates from chapters 12-15 of the reference."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .assembly import assemble_system, hard_simply_supported_dofs, solve_dirichlet
from .benchmarks import simply_supported_sinusoidal_solution, sinusoidal_load
from .element import ShearScheme, q4_element_matrices
from .material import MindlinMaterial
from .mesh import Mesh, rectangular_mesh
from .q4 import kinematic_matrices, mitc4_shear_B


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ElementSpectrum:
    scheme: str
    eigenvalues: FloatArray
    zero_modes: int
    symmetry_error: float
    tolerance: float
    positive_semidefinite: bool

    @property
    def has_only_rigid_zero_modes(self) -> bool:
        return self.zero_modes == 3 and self.positive_semidefinite


@dataclass(frozen=True)
class PatchTestResult:
    scheme: str
    expected_shear: FloatArray
    maximum_curvature_error: float
    maximum_shear_error: float

    @property
    def passed(self) -> bool:
        return max(self.maximum_curvature_error, self.maximum_shear_error) < 1e-10


@dataclass(frozen=True)
class ThicknessScanRow:
    thickness_ratio: float
    analytical_deflection: float
    deflections: dict[str, float]
    normalized_deflections: dict[str, float]


@dataclass(frozen=True)
class DistortionComparison:
    regular_normalized_deflection: float
    distorted_normalized_deflection: float
    relative_change: float
    minimum_gauss_det_jacobian: float


@dataclass(frozen=True)
class CoreValidationReport:
    spectra: dict[str, ElementSpectrum]
    constant_shear: dict[str, PatchTestResult]
    thickness_scan: list[ThicknessScanRow]
    distortion: DistortionComparison


def element_spectrum(
    coordinates: ArrayLike,
    material: MindlinMaterial,
    shear_scheme: ShearScheme,
    relative_tolerance: float = 1e-10,
) -> ElementSpectrum:
    """Check symmetry, semidefiniteness, and physical rigid zero modes."""

    matrix = q4_element_matrices(coordinates, material, shear_scheme).total
    eigenvalues = np.linalg.eigvalsh(matrix)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    tolerance = relative_tolerance * scale
    symmetry_error = float(np.linalg.norm(matrix - matrix.T, ord=np.inf) / scale)
    return ElementSpectrum(
        scheme=shear_scheme,
        eigenvalues=eigenvalues,
        zero_modes=int(np.count_nonzero(np.abs(eigenvalues) <= tolerance)),
        symmetry_error=symmetry_error,
        tolerance=tolerance,
        positive_semidefinite=bool(eigenvalues[0] >= -tolerance),
    )


def constant_shear_patch_test(
    coordinates: ArrayLike,
    shear: ArrayLike = (1.0, -0.5),
    shear_scheme: ShearScheme = "full",
) -> PatchTestResult:
    """Verify w=gamma_x*x+gamma_y*y, theta=0 at multiple element points."""

    coords = np.asarray(coordinates, dtype=float)
    expected = np.asarray(shear, dtype=float)
    if coords.shape != (4, 2) or expected.shape != (2,):
        raise ValueError("expected Q4 coordinates and a two-component shear")
    dofs = np.zeros(12, dtype=float)
    dofs[0::3] = coords @ expected
    value = 1.0 / sqrt(3.0)
    sample_points = [(-value, -value), (value, -value), (value, value), (-value, value)]
    curvature_error = 0.0
    shear_error = 0.0
    for xi, eta in sample_points:
        _, bending, raw_shear = kinematic_matrices(coords, xi, eta)
        shear_matrix = (
            mitc4_shear_B(coords, xi, eta)
            if shear_scheme == "mitc4"
            else raw_shear
        )
        curvature_error = max(curvature_error, float(np.linalg.norm(bending @ dofs)))
        shear_error = max(
            shear_error, float(np.linalg.norm(shear_matrix @ dofs - expected))
        )
    return PatchTestResult(
        scheme=shear_scheme,
        expected_shear=expected,
        maximum_curvature_error=curvature_error,
        maximum_shear_error=shear_error,
    )


def _center_node(mesh: Mesh) -> int:
    center = 0.5 * (mesh.nodes.min(axis=0) + mesh.nodes.max(axis=0))
    return int(np.argmin(np.linalg.norm(mesh.nodes - center, axis=1)))


def _solve_square(
    mesh: Mesh,
    material: MindlinMaterial,
    side: float,
    q0: float,
    scheme: ShearScheme,
) -> float:
    stiffness, force = assemble_system(
        mesh, material, sinusoidal_load(q0, side, side), scheme
    )
    solution = solve_dirichlet(stiffness, force, hard_simply_supported_dofs(mesh))
    return float(solution.displacement[3 * _center_node(mesh)])


def square_thickness_scan(
    side: float,
    young: float,
    poisson: float,
    q0: float,
    ratios: Iterable[float] = (1e-1, 1e-2, 1e-3, 1e-4),
    divisions: int = 6,
    schemes: tuple[ShearScheme, ...] = ("full", "reduced", "mitc4"),
) -> list[ThicknessScanRow]:
    """Run the fixed-mesh thin-limit scan required by section 15.3."""

    mesh = rectangular_mesh(side, side, divisions, divisions)
    rows: list[ThicknessScanRow] = []
    for ratio in ratios:
        if ratio <= 0.0:
            raise ValueError("thickness ratios must be positive")
        material = MindlinMaterial(young, poisson, side * ratio)
        exact = simply_supported_sinusoidal_solution(side, material, q0).center_deflection
        deflections = {
            scheme: _solve_square(mesh, material, side, q0, scheme) for scheme in schemes
        }
        rows.append(
            ThicknessScanRow(
                thickness_ratio=float(ratio),
                analytical_deflection=exact,
                deflections=deflections,
                normalized_deflections={
                    scheme: value / exact for scheme, value in deflections.items()
                },
            )
        )
    return rows


def _minimum_gauss_det_jacobian(mesh: Mesh) -> float:
    value = 1.0 / sqrt(3.0)
    points = [(-value, -value), (value, -value), (value, value), (-value, value)]
    return min(
        kinematic_matrices(mesh.nodes[element], xi, eta)[0].det_jacobian
        for element in mesh.elements
        for xi, eta in points
    )


def square_distortion_comparison(
    side: float,
    material: MindlinMaterial,
    q0: float,
    divisions: int = 6,
    shear_scheme: ShearScheme = "mitc4",
    distortion_x: float = 0.06,
    distortion_y: float = 0.025,
) -> DistortionComparison:
    """Compare regular and smoothly distorted meshes with fixed boundaries."""

    def transform(x: float, y: float) -> tuple[float, float]:
        xi = x / side
        eta = y / side
        envelope = 16.0 * xi * (1.0 - xi) * eta * (1.0 - eta)
        return (
            x + distortion_x * side * envelope,
            y + distortion_y * side * envelope,
        )

    regular = rectangular_mesh(side, side, divisions, divisions)
    distorted = rectangular_mesh(side, side, divisions, divisions, transform=transform)
    exact = simply_supported_sinusoidal_solution(side, material, q0).center_deflection
    regular_w = _solve_square(regular, material, side, q0, shear_scheme)
    distorted_w = _solve_square(distorted, material, side, q0, shear_scheme)
    return DistortionComparison(
        regular_normalized_deflection=regular_w / exact,
        distorted_normalized_deflection=distorted_w / exact,
        relative_change=abs(distorted_w - regular_w) / abs(regular_w),
        minimum_gauss_det_jacobian=_minimum_gauss_det_jacobian(distorted),
    )


def run_core_validation(
    side: float = 1.0,
    young: float = 210e9,
    poisson: float = 0.3,
    q0: float = 10e3,
    divisions: int = 6,
) -> CoreValidationReport:
    """Run the complete reusable validation set for the mathematical core."""

    coordinates = np.array(
        [[-side / 2, -side / 2], [side / 2, -side / 2],
         [side / 2, side / 2], [-side / 2, side / 2]],
        dtype=float,
    )
    material = MindlinMaterial(young, poisson, side / 10.0)
    schemes: tuple[ShearScheme, ...] = ("full", "reduced", "mitc4")
    return CoreValidationReport(
        spectra={scheme: element_spectrum(coordinates, material, scheme) for scheme in schemes},
        constant_shear={
            scheme: constant_shear_patch_test(coordinates, shear_scheme=scheme)
            for scheme in schemes
        },
        thickness_scan=square_thickness_scan(
            side, young, poisson, q0, divisions=divisions, schemes=schemes
        ),
        distortion=square_distortion_comparison(
            side, MindlinMaterial(young, poisson, side / 100.0), q0,
            divisions=divisions, shear_scheme="mitc4"
        ),
    )

