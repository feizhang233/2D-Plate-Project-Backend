"""Minimal structured quadrilateral mesh utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class Mesh:
    nodes: FloatArray
    elements: IntArray

    def __post_init__(self) -> None:
        nodes = np.asarray(self.nodes, dtype=float)
        elements = np.asarray(self.elements, dtype=np.int64)
        if nodes.ndim != 2 or nodes.shape[1] != 2:
            raise ValueError("nodes must have shape (n,2)")
        if elements.ndim != 2 or elements.shape[1] != 4:
            raise ValueError("elements must have shape (m,4)")
        if elements.size and (elements.min() < 0 or elements.max() >= len(nodes)):
            raise ValueError("element connectivity contains an invalid node index")
        object.__setattr__(self, "nodes", nodes)
        object.__setattr__(self, "elements", elements)

    @property
    def ndof(self) -> int:
        return 3 * len(self.nodes)


@dataclass(frozen=True)
class BoundaryEdge:
    element_index: int
    local_edge: int
    nodes: tuple[int, int]


def exterior_edges(mesh: Mesh) -> list[BoundaryEdge]:
    """Return every topological exterior edge exactly once."""

    local_pairs = ((0, 1), (1, 2), (2, 3), (3, 0))
    occurrences: dict[tuple[int, int], list[BoundaryEdge]] = {}
    for element_index, connectivity in enumerate(mesh.elements):
        for local_edge, (start, end) in enumerate(local_pairs):
            node_pair = (int(connectivity[start]), int(connectivity[end]))
            key = tuple(sorted(node_pair))
            occurrences.setdefault(key, []).append(
                BoundaryEdge(element_index, local_edge, node_pair)
            )
    return [items[0] for items in occurrences.values() if len(items) == 1]


def select_exterior_edges(
    mesh: Mesh, predicate: Callable[[FloatArray], bool]
) -> list[BoundaryEdge]:
    """Select exterior edges by a predicate evaluated at each midpoint."""

    return [
        edge
        for edge in exterior_edges(mesh)
        if predicate(mesh.nodes[list(edge.nodes)].mean(axis=0))
    ]


def rectangular_mesh(
    length_x: float,
    length_y: float,
    elements_x: int,
    elements_y: int,
    transform: Callable[[float, float], tuple[float, float]] | None = None,
) -> Mesh:
    """Create a CCW Q4 mesh, optionally transforming every node."""

    if length_x <= 0.0 or length_y <= 0.0:
        raise ValueError("mesh lengths must be positive")
    if elements_x < 1 or elements_y < 1:
        raise ValueError("element counts must be positive")
    nodes: list[tuple[float, float]] = []
    for j in range(elements_y + 1):
        y = length_y * j / elements_y
        for i in range(elements_x + 1):
            x = length_x * i / elements_x
            nodes.append(transform(x, y) if transform is not None else (x, y))

    elements: list[tuple[int, int, int, int]] = []
    stride = elements_x + 1
    for j in range(elements_y):
        for i in range(elements_x):
            lower_left = j * stride + i
            lower_right = lower_left + 1
            upper_left = lower_left + stride
            upper_right = upper_left + 1
            elements.append((lower_left, lower_right, upper_right, upper_left))
    return Mesh(np.asarray(nodes, dtype=float), np.asarray(elements, dtype=np.int64))
