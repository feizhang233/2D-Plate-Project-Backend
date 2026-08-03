# 薄板 K / 厚板 M 四边形板数学核心

本项目提供只依赖 NumPy、可逐步验证的四边形板有限元数学核心：

- **K 方法**：离散 Kirchhoff 四边形（DKQ），用于薄板，不计横向剪切变形；
- **M 方法**：Reissner-Mindlin Q4/MITC4，用于中厚板，计入横向剪切变形；
- **自动方法**：根据整个板平面形状的最小跨长与厚度之比选择 K 或 M。

两种方法共用自由度约定：

\[
\mathbf a_i=[w_i,\theta_{xi},\theta_{yi}]^T,
\qquad
\boldsymbol\gamma=
\begin{bmatrix}w_{,x}-\theta_x\\w_{,y}-\theta_y\end{bmatrix}.
\]

## 离散 K 方法数学核心

DKQ 不增加边中点自由度。对每条由节点 \(i,j\) 构成、长度为 \(L\) 的边，在边中点离散施加 Kirchhoff 约束。以 \(s,n\) 表示边的切向和法向：

\[
w_{,s}^{m}=\frac{3}{2L}(w_j-w_i)
-\frac14\left(w_{,s}^{i}+w_{,s}^{j}\right),
\qquad
w_{,n}^{m}=\frac12\left(w_{,n}^{i}+w_{,n}^{j}\right).
\]

四个角点和四个受约束边中点的斜率使用 Q8 Serendipity 场插值，得到

\[
\boldsymbol\kappa_K=
\begin{bmatrix}w_{,xx}&w_{,yy}&2w_{,xy}\end{bmatrix}^{T}
=\mathbf B_K\mathbf a_e,
\qquad
\mathbf K_K^e=\int_{A_e}\mathbf B_K^T\mathbf D_b\mathbf B_K\,dA.
\]

K 方法的横向剪切应变、剪力和剪切能严格为零；单元仍为 4 节点、每节点 3 自由度，可直接复用现有网格、装配、边界条件和后处理。

## K / M 自动选择

`plate_characteristic_length` 取板平面凸包的旋转不变最小宽度 \(L_c\)，因此板旋转、网格加密或长宽比变化不会误用单元尺寸作为判据。默认规则为：

\[
\frac{t}{L_c}\le\frac1{20}\Rightarrow K,
\qquad
\frac{t}{L_c}>\frac1{20}\Rightarrow M\;(\mathrm{MITC4}).
\]

阈值可通过 `thinness_threshold` 调整，也可用 `plate_method="K"` 或 `plate_method="M"` 强制指定：

```python
from mindlin_plate import (
    MindlinMaterial,
    assemble_plate_system,
    rectangular_mesh,
)

mesh = rectangular_mesh(2.0, 1.0, 8, 4)
material = MindlinMaterial(young=210e9, poisson=0.3, thickness=0.02)
system = assemble_plate_system(mesh, material, load=10e3)

print(system.selection.method)           # "K"
print(system.selection.thickness_ratio)  # 0.02 / 1.0
K, f = system.stiffness, system.force
```

低层 `assemble_system` 为兼容原有调用仍默认 M 方法；传入 `plate_method="auto"` 也可启用相同的自动判据。

## 9 个累计步骤

| Step | 新增数学能力 | 对应参考文档 | 累计示例 |
|---|---|---|---|
| 1 | Mindlin 运动学、材料矩阵、正弦载荷解析解 | 第 2-4、11 章 | `examples/step_01_continuum.py` |
| 2 | Q4 形函数、Jacobian、\(B_b\)、原始 \(B_s\) | 第 6-7 章 | `examples/step_02_q4_kinematics.py` |
| 3 | 弯曲/剪切刚度、一致载荷、单元能量 | 第 6.4-6.5、7.3 章 | `examples/step_03_element.py` |
| 4 | 结构化网格、全局组装、边界条件与求解 | 第 5、10 章 | `examples/step_04_global_solver.py` |
| 5 | 选择性减缩剪切积分与厚跨比扫描 | 第 8 章 | `examples/step_05_selective_integration.py` |
| 6 | MITC4 协变剪切、tying points | 第 9、13 章 | `examples/step_06_mitc4.py` |
| 7 | 弯矩、剪力、表面应力恢复与畸变网格检查 | 第 4.3、15 章 | `examples/step_07_recovery_validation.py` |
| 8 | 自然边界载荷、对称/倾斜边及 Gauss 点恢复 | 第 5、10 章 | `examples/step_08_boundaries_postprocess.py` |
| 9 | 全部验证闸门和 Project 完成审计 | 第 12-16 章 | `examples/step_09_complete_validation.py` |

每个示例都是累计示例。例如 Step 9 会依次运行 Step 1-9 的检查，而不是只运行最终验证片段。

## 运行

```bash
python3 examples/step_01_continuum.py
python3 examples/step_09_complete_validation.py
python3 -m unittest discover -s tests -v
```

如需使用指定 Python：

```bash
/path/to/python3 examples/step_09_complete_validation.py
```

## FastAPI 服务

HTTP 服务把网格前处理、K/M 数学核心、边界约束、线性求解、结果恢复和 PNG
后处理串成一次请求。安装并启动：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/mindlin-plate-api
```

启动后可访问：

- Swagger UI：`http://127.0.0.1:8000/docs`
- 健康检查：`GET /health`
- 前处理模板：`GET /api/v1/templates/rectangular-plate`
- 提交计算：`POST /api/v1/analyses`
- 查询本进程内的结果：`GET /api/v1/analyses/{analysis_id}`
- 后处理图像：响应中 `images[].url` 指向 PNG 静态文件

可直接提交随包提供的前处理模板
[`mindlin_plate/templates/rectangular_plate.json`](mindlin_plate/templates/rectangular_plate.json)：

```bash
curl -s http://127.0.0.1:8000/api/v1/templates/rectangular-plate \
  > plate.json
curl -s -X POST http://127.0.0.1:8000/api/v1/analyses \
  -H 'Content-Type: application/json' \
  --data-binary @plate.json \
  > result.json
```

前处理支持矩形结构网格和自定义 `nodes/elements` Q4 网格，均布或正弦分布面荷载，
固支、软简支、硬简支或自定义自由度约束。`plate_method` 可选 `auto/K/M`；
M 方法的 `shear_scheme` 可选 `full/reduced/mitc4`。单位字段仅作记录，所有输入仍须使用
同一套自洽单位。

响应包含计算摘要、节点位移/转角/反力、单元中心曲率/弯矩/剪力/上下表面应力，
并可生成以下图像：

- `deflection.png`：节点横向挠度云图；
- `rotation.png`：转角合量云图；
- `moment.png`：单元中心弯矩合量图；
- `stress_top.png`：上表面等效弯曲应力图。

默认图片目录为 `outputs/plate-analyses`。部署时可用
`MINDLIN_PLATE_OUTPUT_DIR`、`MINDLIN_PLATE_HOST` 和 `MINDLIN_PLATE_PORT`
修改输出路径、监听地址和端口。由于当前核心使用教学用途的稠密矩阵，服务层将单次任务限制为
2500 个自由度。

## 核心约定

- Q4 节点顺序：左下、右下、右上、左上；所有积分点都要求 `det(J) > 0`。
- 弯曲使用 \(2\times2\) Gauss 积分。
- `full`：原始剪切场使用 \(2\times2\) 积分，是会剪切锁死的基线。
- `reduced`：原始剪切场使用中心单点积分。
- `mitc4`：在四个 tying points 对协变剪切分量插值，再用 \(2\times2\) 积分。
- `K` / `dkq`：离散 Kirchhoff 薄板，仅积分弯曲刚度，使用 \(2\times2\) Gauss 积分。
- 单点 `reduced` 单元保留两个额外零能模态，仅作为锁死机理基线；正式方案使用通过零模态检查的 `mitc4`。
- 自然边界支持 \([\bar V,\bar m_x,\bar m_y]\)，倾斜边通过局部正交基底施加约束。
- 代码只构造数学方程；单位必须由调用者保持一致。

## 目录

- `mindlin_plate/kirchhoff.py`：DKQ 斜率插值、边约束和曲率矩阵。
- `mindlin_plate/theory.py`：板形状特征长度与 K/M 自动判据。
- `mindlin_plate/` 其余模块：材料、M 方法、装配、边界和后处理。
- `examples/`：9 个累计示例和共用检查器。
- `tests/`：公式级、单元级和全局级回归测试。

## 完成标准对应

- 完整积分 Q4 基准：`shear_scheme="full"`，可重现剪切锁死。
- 稳定厚薄板方案：`shear_scheme="mitc4"`，自由单元只有 3 个物理刚体零模态。
- Patch tests：刚体、纯扭曲、常剪切和参考纯弯模式均已自动化。
- 薄板极限：固定网格扫描 \(t/L=10^{-1},10^{-2},10^{-3},10^{-4}\)。
- 网格畸变：检查全部 Gauss 点 `det(J)>0`，并量化规则/畸变网格结果变化。
- 输出：节点撓度/转角、曲率、弯矩、剪力、上下表面弯曲应力和抛物线剪应力恢复。
- 边界：固支、软/硬简支、对称边、倾斜边局部自由度及分布边界剪力/弯矩。
