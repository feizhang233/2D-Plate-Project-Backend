# Mindlin Plate Core & FastAPI

A small, verifiable finite-element library for quadrilateral plate analysis,
with a FastAPI service for JSON-based preprocessing, solving, result recovery,
and PNG post-processing.

## Interview study guide (中文)

如果目的是准备博士面试，而不是立即阅读实现细节，请从
[`START_HERE_博士面试.md`](START_HERE_博士面试.md) 开始。配套材料包括：

- [`docs/博士面试_数学核心速览.md`](docs/博士面试_数学核心速览.md)：用一条主线理解本项目的数学核心。
- [`docs/博士面试_高频问答.md`](docs/博士面试_高频问答.md)：用于口头复述和自测。
- [`docs/README.md`](docs/README.md)：项目文件地图和分层学习路线。
- [`output/pdf/README.md`](output/pdf/README.md)：两份完整算例 PDF 的区别和使用方式。

The numerical core supports two compatible 4-node, three-degree-of-freedom
plate formulations:

- **K / DKQ** — discrete Kirchhoff quadrilateral for thin plates; transverse
  shear deformation is not included.
- **M / Reissner-Mindlin Q4** — a shear-deformable plate model with `full`,
  `reduced`, and `mitc4` shear integration options.
- **Automatic selection** — chooses K or M from the plate planform and
  thickness, while still allowing an explicit method override.

## Features

- Structured rectangular meshes and custom Q4 meshes.
- Isotropic material definition with bending and shear constitutive matrices.
- Consistent uniform or sinusoidal transverse loads.
- Clamped, soft simply supported, hard simply supported, and custom
  displacement boundary conditions.
- Element-center recovery of curvature, shear strain, moments, shear forces,
  and top/bottom surface stresses.
- FastAPI endpoints, OpenAPI documentation, an editable JSON request template,
  and static PNG result artifacts.
- Formula-, element-, and system-level regression tests.

## Formulations and automatic selection

Each node has the following degrees of freedom:

$$
\mathbf{a}_i = [w_i, \theta_{xi}, \theta_{yi}]^T,
\qquad
\boldsymbol{\gamma} =
\begin{bmatrix}
w_{,x} - \theta_x \\
w_{,y} - \theta_y
\end{bmatrix}.
$$

The automatic selector uses the rotation-invariant minimum width of the plate
convex hull, $L_c$, rather than an element size. Its default rule is:

$$
\frac{t}{L_c} \leq \frac{1}{20} \Rightarrow \text{K (DKQ)},
\qquad
\frac{t}{L_c} > \frac{1}{20} \Rightarrow \text{M (MITC4)}.
$$

Set `plate_method` to `"K"` or `"M"` to override the selector, or adjust
`thinness_threshold` to change the cutoff.

## Installation

Python 3.10 or later is required.

```bash
git clone https://github.com/<your-account>/mindlin-plate-fastapi.git
cd mindlin-plate-fastapi

python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
```

Use `-e .` instead if you only need the library and API runtime dependencies.

## Python usage

```python
from mindlin_plate import MindlinMaterial, assemble_plate_system, rectangular_mesh

mesh = rectangular_mesh(length_x=2.0, length_y=1.0, elements_x=8, elements_y=4)
material = MindlinMaterial(young=210e9, poisson=0.3, thickness=0.02)
system = assemble_plate_system(mesh, material, load=10e3)

print(system.selection.method)           # "K"
print(system.selection.thickness_ratio)  # 0.02 / 1.0
stiffness, force = system.stiffness, system.force
```

`assemble_system` keeps the legacy M-method default. Use
`assemble_system(..., plate_method="auto")` to enable the same automatic
selection used by `assemble_plate_system`.

## FastAPI service

Start the local API server:

```bash
.venv/bin/mindlin-plate-api
```

The default server address is `http://127.0.0.1:8000`. Open the interactive
API documentation at `http://127.0.0.1:8000/docs`.

| Endpoint | Purpose |
| --- | --- |
| `GET /health` | Health check. |
| `GET /api/v1/templates/rectangular-plate` | Returns the editable request template. |
| `POST /api/v1/analyses` | Runs an analysis and creates PNG artifacts. |
| `GET /api/v1/analyses/{analysis_id}` | Returns an analysis created by the current server process. |

The request template is also available in the repository at
[`mindlin_plate/templates/rectangular_plate.json`](mindlin_plate/templates/rectangular_plate.json).

```bash
curl -s http://127.0.0.1:8000/api/v1/templates/rectangular-plate -o plate.json

curl -s -X POST http://127.0.0.1:8000/api/v1/analyses \
  -H 'Content-Type: application/json' \
  --data-binary @plate.json \
  -o result.json
```

`POST /api/v1/solve` is retained as an undocumented alias for
`POST /api/v1/analyses`.

### Request model

The JSON body contains the following top-level objects:

| Field | Description |
| --- | --- |
| `mesh` | A `rectangular` mesh (`length_x`, `length_y`, `elements_x`, `elements_y`) or a `custom` mesh (`nodes`, `elements`). |
| `material` | `young`, `poisson`, `thickness`, and optional `shear_correction`. |
| `load` | A `uniform` or `sinusoidal` transverse distributed load with `magnitude`. |
| `boundary` | `clamped`, `soft_simply_supported`, `hard_simply_supported`, or `custom` with `prescribed_dofs`. |
| `analysis` | `plate_method`, `shear_scheme`, and `thinness_threshold`. |
| `postprocess` | Requested plots and PNG `dpi`. |

The response contains the formulation decision, equilibrium summary, nodal
displacements and reactions, element-center recovered fields, and artifact
URLs. Generated images are served under `/artifacts/<analysis_id>/`.

### Post-processing images

Choose any subset of the following values in `postprocess.plots`:

- `deflection` — transverse deflection contour.
- `rotation` — rotation magnitude contour.
- `moment` — element-center bending resultant magnitude.
- `stress_top` — top-surface equivalent bending stress.

Artifacts are written to `outputs/plate-analyses` by default. Configure the
runtime with these environment variables when needed:

- `MINDLIN_PLATE_OUTPUT_DIR`
- `MINDLIN_PLATE_HOST`
- `MINDLIN_PLATE_PORT`

## Modelling notes

- Q4 nodes must be ordered counter-clockwise: lower-left, lower-right,
  upper-right, upper-left.
- The Jacobian determinant must be positive at all integration points.
- Bending is integrated with a $2 \times 2$ Gauss rule.
- `full` uses a $2 \times 2$ rule for the raw shear field and intentionally
  reproduces shear locking as a baseline.
- `reduced` uses one-point raw-shear integration and retains two extra
  zero-energy modes; use it only as a locking reference.
- `mitc4` interpolates covariant shear from four tying points and is the
  recommended M-method option.
- DKQ has zero transverse shear strain, shear force, and shear energy.
- The HTTP service uses a dense global stiffness matrix and limits each
  request to 2,500 degrees of freedom.
- Inputs must use one consistent system of units. The API `units` object is
  metadata only; it does not convert values.

## Verification and examples

Run the complete test suite:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

The `examples/` directory contains nine cumulative scripts. Each later step
reruns the checks from all preceding steps.

| Step | Added capability | Script |
| --- | --- | --- |
| 1 | Mindlin kinematics, material matrices, and sinusoidal reference solution | `examples/step_01_continuum.py` |
| 2 | Q4 shape functions, Jacobian, bending, and raw shear matrices | `examples/step_02_q4_kinematics.py` |
| 3 | Element stiffness, consistent loading, and energy | `examples/step_03_element.py` |
| 4 | Structured mesh, global assembly, constraints, and solve | `examples/step_04_global_solver.py` |
| 5 | Selective reduced integration and thickness scan | `examples/step_05_selective_integration.py` |
| 6 | MITC4 covariant shear and tying points | `examples/step_06_mitc4.py` |
| 7 | Moment, shear, stress recovery, and distorted-mesh checks | `examples/step_07_recovery_validation.py` |
| 8 | Natural edge loads, symmetry, inclined edges, and Gauss-point recovery | `examples/step_08_boundaries_postprocess.py` |
| 9 | Complete validation audit | `examples/step_09_complete_validation.py` |

## Project layout

```text
mindlin_plate/
  assembly.py       Global assembly and Dirichlet solver
  element.py        Q4, MITC4, and DKQ element operations
  theory.py         Planform characteristic length and K/M selection
  postprocess.py    Field recovery
  service.py        JSON-like input to core-analysis adapter
  plotting.py       Headless PNG rendering
  api.py            FastAPI application
  templates/        Request templates
examples/           Incremental validation examples
tests/              Regression tests
```
