"""One checkpoint per project step; later examples run all earlier checkpoints."""

from __future__ import annotations

import numpy as np

from mindlin_plate import (
    MindlinMaterial,
    assemble_boundary_load,
    assemble_system,
    edge_angle,
    exterior_edges,
    global_dof_transform,
    hard_simply_supported_dofs,
    kinematic_matrices,
    local_edge_constraints,
    q4_consistent_load,
    q4_edge_consistent_load,
    q4_element_matrices,
    recover_element_centers,
    recover_element_gauss_points,
    rectangular_mesh,
    run_core_validation,
    select_exterior_edges,
    simply_supported_sinusoidal_solution,
    sinusoidal_load,
    split_nodal_dofs,
    solve_dirichlet,
    transform_system_to_local,
)


def _title(step: int, name: str) -> None:
    print(f"\n[Step {step}] {name}")


def _sample_material(thickness: float = 0.1) -> MindlinMaterial:
    return MindlinMaterial(young=210e9, poisson=0.3, thickness=thickness)


def checkpoint_1() -> None:
    _title(1, "连续体运动学、本构与解析基线")
    material = _sample_material()
    answer = simply_supported_sinusoidal_solution(1.0, material, 10e3)
    print(f"D = {material.bending_rigidity:.6e} N m")
    print(f"S = {material.shear_rigidity:.6e} N/m")
    print(f"Wb = {answer.bending_deflection:.6e} m")
    print(f"Ws = {answer.shear_deflection:.6e} m")
    print(f"W  = {answer.center_deflection:.6e} m")
    assert np.isclose(answer.center_deflection, 1.40984e-6, rtol=8e-6)


def _square_coordinates(length: float = 1.0) -> np.ndarray:
    half = length / 2.0
    return np.array(
        [[-half, -half], [half, -half], [half, half], [-half, half]], dtype=float
    )


def _nodal_field(coordinates: np.ndarray, field) -> np.ndarray:
    values = np.zeros(12, dtype=float)
    for i, (x, y) in enumerate(coordinates):
        values[3 * i : 3 * i + 3] = field(float(x), float(y))
    return values


def checkpoint_2() -> None:
    _title(2, "Q4 映射、B 矩阵与 patch 场")
    coords = _square_coordinates()
    dofs = _nodal_field(coords, lambda x, y: (x * y, y, x))
    max_shear = 0.0
    for xi, eta in [(-1 / np.sqrt(3), -1 / np.sqrt(3)), (0.0, 0.0), (0.4, -0.2)]:
        point, bending, shear = kinematic_matrices(coords, xi, eta)
        max_shear = max(max_shear, float(np.linalg.norm(shear @ dofs)))
        assert np.isclose(point.shape.sum(), 1.0)
        assert point.det_jacobian > 0.0
        assert np.allclose(bending @ dofs, [0.0, 0.0, 2.0], atol=1e-12)
    print(f"center det(J) = {kinematic_matrices(coords, 0, 0)[0].det_jacobian:.6f}")
    print(f"pure-twist max |gamma| = {max_shear:.3e}")


def checkpoint_3() -> None:
    _title(3, "单元刚度、一致载荷与刚体零模态")
    coords = _square_coordinates()
    material = _sample_material()
    matrices = q4_element_matrices(coords, material, "full")
    force = q4_consistent_load(coords, 10e3)
    eigenvalues = np.linalg.eigvalsh(matrices.total)
    tolerance = max(abs(eigenvalues[-1]) * 1e-10, 1e-8)
    zero_modes = int(np.count_nonzero(np.abs(eigenvalues) < tolerance))
    print(f"sum nodal transverse load = {force[0::3].sum():.6e} N")
    print(f"symmetry error = {np.linalg.norm(matrices.total - matrices.total.T):.3e}")
    print(f"numerical zero modes = {zero_modes}")
    assert np.isclose(force[0::3].sum(), 10e3)
    assert zero_modes == 3


def _center_node(mesh) -> int:
    center = 0.5 * (mesh.nodes.min(axis=0) + mesh.nodes.max(axis=0))
    return int(np.argmin(np.linalg.norm(mesh.nodes - center, axis=1)))


def _solve_square(thickness: float, scheme: str, divisions: int = 6):
    material = _sample_material(thickness)
    mesh = rectangular_mesh(1.0, 1.0, divisions, divisions)
    load = sinusoidal_load(10e3, 1.0, 1.0)
    stiffness, force = assemble_system(mesh, material, load, scheme)
    solution = solve_dirichlet(stiffness, force, hard_simply_supported_dofs(mesh))
    w = solution.displacement[3 * _center_node(mesh)]
    exact = simply_supported_sinusoidal_solution(1.0, material, 10e3).center_deflection
    return mesh, material, solution, w, exact


def checkpoint_4() -> None:
    _title(4, "全局组装、边界条件与线性求解")
    mesh, _, solution, w, exact = _solve_square(0.1, "full", divisions=4)
    residual = np.linalg.norm(
        solution.reactions[solution.free_dofs], ord=np.inf
    )
    print(f"mesh = {len(mesh.elements)} elements, {mesh.ndof} dofs")
    print(f"center w(full) = {w:.6e} m; analytical = {exact:.6e} m")
    print(f"free-dof residual = {residual:.3e} N")
    assert w > 0.0 and residual < 1e-4


def checkpoint_5() -> None:
    _title(5, "选择性减缩积分与剪切锁死扫描")
    ratios = [1e-1, 1e-2, 1e-3]
    last = None
    print("t/L       full/exact   reduced/exact")
    for ratio in ratios:
        _, _, _, w_full, exact = _solve_square(ratio, "full")
        _, _, _, w_reduced, _ = _solve_square(ratio, "reduced")
        last = (w_full / exact, w_reduced / exact)
        print(f"{ratio:1.0e}    {w_full/exact:10.4f}   {w_reduced/exact:12.4f}")
    assert last is not None and last[1] > 5.0 * last[0]


def checkpoint_6() -> None:
    _title(6, "MITC4 tying-point 剪切场")
    length = 1.0
    curvature = 2.0
    coords = _square_coordinates(length)
    material = _sample_material(0.01)
    dofs = _nodal_field(coords, lambda x, y: (0.5 * curvature * x * x, curvature * x, 0.0))
    energies: dict[str, float] = {}
    for scheme in ["full", "reduced", "mitc4"]:
        matrices = q4_element_matrices(coords, material, scheme)
        energies[scheme] = 0.5 * float(dofs @ matrices.shear @ dofs)
    _, _, _, w_mitc, exact = _solve_square(0.01, "mitc4")
    print("pure-bending shear energies:", ", ".join(f"{k}={v:.3e}" for k, v in energies.items()))
    print(f"t/L=1e-2 center w(MITC4)/analytical = {w_mitc/exact:.4f}")
    assert energies["full"] > 0.0
    assert abs(energies["mitc4"]) < energies["full"] * 1e-12


def checkpoint_7() -> None:
    _title(7, "结果恢复、表面应力与畸变网格")

    def skew(x: float, y: float) -> tuple[float, float]:
        # The boundary stays axis-aligned; only interior rows are shifted.
        envelope = 16.0 * x * (1.0 - x) * y * (1.0 - y)
        return x + 0.06 * envelope, y + 0.025 * envelope

    mesh = rectangular_mesh(1.0, 1.0, 4, 4, transform=skew)
    material = _sample_material(0.05)
    stiffness, force = assemble_system(
        mesh, material, sinusoidal_load(10e3, 1.0, 1.0), "mitc4"
    )
    solution = solve_dirichlet(stiffness, force, hard_simply_supported_dofs(mesh))
    recovered = recover_element_centers(mesh, solution.displacement, material, "mitc4")
    target = min(
        recovered,
        key=lambda item: np.linalg.norm(item.response.physical_point - np.array([0.5, 0.5])),
    ).response
    min_det = min(
        kinematic_matrices(mesh.nodes[e], 0.0, 0.0)[0].det_jacobian
        for e in mesh.elements
    )
    print(f"minimum center det(J) = {min_det:.6e}")
    print(f"near-center M = {np.array2string(target.bending_moment, precision=4)} N")
    print(f"near-center Q = {np.array2string(target.shear_force, precision=4)} N/m")
    print(f"top stress = {np.array2string(target.stress_top, precision=4)} Pa")
    assert min_det > 0.0
    assert np.allclose(target.stress_top, -target.stress_bottom)


def checkpoint_8() -> None:
    _title(8, "自然边界载荷、局部坐标约束与完整恢复")
    coords = _square_coordinates()
    edge_force = q4_edge_consistent_load(
        coords, edge=0, transverse_shear=10.0, moment=(2.0, -3.0)
    )
    totals = edge_force.reshape(4, 3).sum(axis=0)
    assert np.allclose(totals, [10.0, 2.0, -3.0])

    mesh = rectangular_mesh(1.0, 1.0, 2, 2)
    top_edges = select_exterior_edges(mesh, lambda midpoint: np.isclose(midpoint[1], 1.0))
    boundary_force = assemble_boundary_load(mesh, top_edges, transverse_shear=5.0)
    assert len(exterior_edges(mesh)) == 8
    assert len(top_edges) == 2
    assert np.isclose(boundary_force[0::3].sum(), 5.0)

    angle = edge_angle([0.0, 0.0], [1.0, 1.0])
    transform = global_dof_transform(4, {node: angle for node in range(4)})
    material = _sample_material()
    stiffness = q4_element_matrices(coords, material, "mitc4").total
    local_stiffness, _ = transform_system_to_local(stiffness, np.zeros(12), transform)
    local_dofs = np.linspace(-0.3, 0.4, 12)
    global_dofs = transform @ local_dofs
    assert np.isclose(
        global_dofs @ stiffness @ global_dofs,
        local_dofs @ local_stiffness @ local_dofs,
    )
    assert set(local_edge_constraints([0, 1], "symmetry")) == {2, 5}

    stiffness_global, force_global = assemble_system(
        mesh, material, sinusoidal_load(10e3, 1.0, 1.0), "mitc4"
    )
    solution = solve_dirichlet(
        stiffness_global, force_global, hard_simply_supported_dofs(mesh)
    )
    w, theta = split_nodal_dofs(solution.displacement)
    gauss_results = recover_element_gauss_points(
        mesh, solution.displacement, material, "mitc4"
    )
    print(f"single-edge resultant [V,mx,my] = {np.array2string(totals, precision=4)}")
    print(f"exterior/top edges = {len(exterior_edges(mesh))}/{len(top_edges)}")
    print(f"recovered nodal w/theta shapes = {w.shape}/{theta.shape}")
    print(f"Gauss-point result count = {len(gauss_results)}")
    assert len(gauss_results) == 4 * len(mesh.elements)


def checkpoint_9() -> None:
    _title(9, "完整验证闸门与 Project 收口")
    report = run_core_validation(divisions=6)
    print("scheme    zero modes   constant-shear patch")
    for scheme in ("full", "reduced", "mitc4"):
        print(
            f"{scheme:8s}  {report.spectra[scheme].zero_modes:5d}"
            f"        {'PASS' if report.constant_shear[scheme].passed else 'FAIL'}"
        )
    print("t/L       full/exact   reduced/exact   mitc4/exact")
    for row in report.thickness_scan:
        normalized = row.normalized_deflections
        print(
            f"{row.thickness_ratio:1.0e}    {normalized['full']:10.6f}"
            f"   {normalized['reduced']:12.6f}   {normalized['mitc4']:11.6f}"
        )
    distortion = report.distortion
    print(
        "distorted MITC4: normalized w="
        f"{distortion.distorted_normalized_deflection:.6f}, "
        f"change={100*distortion.relative_change:.3f}%, "
        f"min det(J)={distortion.minimum_gauss_det_jacobian:.6e}"
    )
    assert report.spectra["full"].has_only_rigid_zero_modes
    assert report.spectra["mitc4"].has_only_rigid_zero_modes
    # One-point shear integration intentionally remains an hourglass-prone baseline.
    assert report.spectra["reduced"].zero_modes > 3
    assert all(result.passed for result in report.constant_shear.values())
    last = report.thickness_scan[-1].normalized_deflections
    assert last["mitc4"] > 0.95 and last["full"] < 1e-4
    assert distortion.minimum_gauss_det_jacobian > 0.0
    assert distortion.relative_change < 0.05


CHECKPOINTS = [
    checkpoint_1,
    checkpoint_2,
    checkpoint_3,
    checkpoint_4,
    checkpoint_5,
    checkpoint_6,
    checkpoint_7,
    checkpoint_8,
    checkpoint_9,
]


def run_through(last_step: int) -> None:
    if not 1 <= last_step <= len(CHECKPOINTS):
        raise ValueError("last_step is out of range")
    print(f"累计示例：执行 Step 1-{last_step}")
    for checkpoint in CHECKPOINTS[:last_step]:
        checkpoint()
    print(f"\nStep 1-{last_step} 全部通过")
