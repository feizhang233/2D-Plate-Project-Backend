"""Headless PNG post-processing for plate analyses."""

from __future__ import annotations

import os
import tempfile
from threading import Lock
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .service import PlateAnalysis


AVAILABLE_PLOTS = ("deflection", "rotation", "moment", "stress_top")
_RENDER_LOCK = Lock()


def _plot_modules():
    cache = Path(tempfile.gettempdir()) / "mindlin-plate-matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.collections import PolyCollection
        from matplotlib.figure import Figure
        from matplotlib.tri import Triangulation
    except ImportError as exc:  # pragma: no cover - exercised only in minimal installs
        raise RuntimeError(
            "matplotlib is required for post-processing images"
        ) from exc
    return FigureCanvasAgg, PolyCollection, Figure, Triangulation


def _triangles(elements: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            triangle
            for element in elements
            for triangle in (element[[0, 1, 2]], element[[0, 2, 3]])
        ],
        dtype=np.int64,
    )


def _finish_figure(figure, axes, canvas, path: Path, title: str) -> None:
    axes.set_title(title)
    axes.set_xlabel("x")
    axes.set_ylabel("y")
    axes.set_aspect("equal", adjustable="box")
    figure.tight_layout()
    canvas.print_png(str(path))


def _nodal_plot(
    analysis: PlateAnalysis,
    values: np.ndarray,
    path: Path,
    title: str,
    colorbar_label: str,
    dpi: int,
) -> None:
    FigureCanvasAgg, PolyCollection, Figure, Triangulation = _plot_modules()
    figure = Figure(figsize=(7.2, 5.4), dpi=dpi)
    canvas = FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    mesh = analysis.mesh
    triangulation = Triangulation(
        mesh.nodes[:, 0], mesh.nodes[:, 1], _triangles(mesh.elements)
    )
    if len(mesh.elements) > 1 and float(np.ptp(values)) > np.finfo(float).eps:
        image = axes.tricontourf(triangulation, values, levels=20, cmap="viridis")
    else:
        image = axes.tripcolor(triangulation, values, shading="flat", cmap="viridis")
    axes.add_collection(
        PolyCollection(
            [mesh.nodes[element] for element in mesh.elements],
            facecolors="none",
            edgecolors="black",
            linewidths=0.35,
            alpha=0.45,
        )
    )
    colorbar = figure.colorbar(image, ax=axes, shrink=0.88)
    colorbar.set_label(colorbar_label)
    _finish_figure(figure, axes, canvas, path, title)


def _element_plot(
    analysis: PlateAnalysis,
    values: np.ndarray,
    path: Path,
    title: str,
    colorbar_label: str,
    dpi: int,
) -> None:
    FigureCanvasAgg, PolyCollection, Figure, _ = _plot_modules()
    figure = Figure(figsize=(7.2, 5.4), dpi=dpi)
    canvas = FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    polygons = [analysis.mesh.nodes[element] for element in analysis.mesh.elements]
    collection = PolyCollection(
        polygons,
        array=np.asarray(values),
        cmap="viridis",
        edgecolors="black",
        linewidths=0.45,
    )
    axes.add_collection(collection)
    axes.autoscale_view()
    colorbar = figure.colorbar(collection, ax=axes, shrink=0.88)
    colorbar.set_label(colorbar_label)
    _finish_figure(figure, axes, canvas, path, title)


def _equivalent_stress(stress: np.ndarray) -> float:
    sx, sy, txy = stress
    return float(np.sqrt(max(sx * sx - sx * sy + sy * sy + 3.0 * txy * txy, 0.0)))


def render_postprocess_images(
    analysis: PlateAnalysis,
    output_directory: str | Path,
    plots: Iterable[str] = AVAILABLE_PLOTS,
    dpi: int = 140,
) -> list[Path]:
    """Render the selected result fields and return their PNG paths."""

    requested = list(dict.fromkeys(str(plot).strip().lower() for plot in plots))
    unknown = sorted(set(requested) - set(AVAILABLE_PLOTS))
    if unknown:
        raise ValueError(f"unknown post-processing plots: {', '.join(unknown)}")
    if not 72 <= int(dpi) <= 300:
        raise ValueError("postprocess.dpi must be between 72 and 300")
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    nodal = analysis.solution.displacement.reshape(-1, 3)
    recovered = [item.response for item in analysis.recovered]
    renderers: dict[str, Callable[[Path], None]] = {
        "deflection": lambda path: _nodal_plot(
            analysis, nodal[:, 0], path, "Transverse deflection w", "w", int(dpi)
        ),
        "rotation": lambda path: _nodal_plot(
            analysis,
            np.linalg.norm(nodal[:, 1:3], axis=1),
            path,
            "Rotation magnitude",
            "sqrt(theta_x^2 + theta_y^2)",
            int(dpi),
        ),
        "moment": lambda path: _element_plot(
            analysis,
            np.asarray([np.linalg.norm(item.bending_moment) for item in recovered]),
            path,
            "Bending moment resultant magnitude",
            "||M||",
            int(dpi),
        ),
        "stress_top": lambda path: _element_plot(
            analysis,
            np.asarray([_equivalent_stress(item.stress_top) for item in recovered]),
            path,
            "Top-surface equivalent bending stress",
            "sigma_eq,top",
            int(dpi),
        ),
    }
    paths: list[Path] = []
    # Matplotlib keeps process-global font and style caches.  Serialize only
    # rendering; independent finite-element solves may still run concurrently.
    with _RENDER_LOCK:
        for name in requested:
            path = directory / f"{name}.png"
            renderers[name](path)
            paths.append(path)
    return paths
