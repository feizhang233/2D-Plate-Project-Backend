# Mindlin Plate 教材题目与答案

本文档提取第 12 章中适合 Reissner-Mindlin 板有限元验证的题目，并将书中对应答案直接放在每道题下面。为便于复核，同时保留英文原文、书本页码和拆分 PDF 页码。

## 1. Example 12.1：四边固支方板中心集中力

### 1.1 位置

- 题目及 Figure 12-9：`第12章/12.4.pdf` 文件第 2 页，书本第 585 页。
- 书中解答及 Figure 12-10：`第12章/12.4.pdf` 文件第 3 页，书本第 586 页。

### 1.2 英文原文

> The problem of a square steel plate fixed along all four edges and subjected to a concentrated load at its center is shown in Figure 12-9. Determine the maximum vertical deflection of the plate.

### 1.3 中文题目及图中参数

一块正方形钢板沿四条边完全固支，并在板中心承受集中载荷。求板的最大竖向挠度。

Figure 12-9 给出的参数为：

| 项目 | 数值 |
|---|---:|
| 方板边长 $L$ | $20\ \mathrm{in.}$ |
| 板厚 $t$ | $0.1\ \mathrm{in.}$ |
| 中心集中力 $P$ | $100\ \mathrm{lb}$，向下 |
| 网格 | $2\times2$ |
| 边界 | 四边固支 |

书中解答采用：

$$
E=30\times10^6\ \mathrm{psi},\qquad \nu=0.3
$$

对应板弯曲刚度：

$$
D=\frac{Et^3}{12(1-\nu^2)}
=2.747\times10^3\ \mathrm{lb\,in.}
$$

### 1.4 书中答案

书中 $2\times2$ 网格计算得到中心最大竖向位移：

$$
w_{\max,h}=-0.07583\ \mathrm{in.}
$$

书中同时给出经典薄板解：

$$
w_{\max}
=0.0056\frac{PL^2}{D}
=-0.0815\ \mathrm{in.}
$$

书中说明：继续细化到 $4\times4$ 网格，计算结果应向经典解收敛。

### 1.5 后续验证设置

- 节点自由度顺序建议统一为 $[w,\theta_x,\theta_y]$。
- 四边固支应施加 $w=\theta_x=\theta_y=0$。
- $2\times2$ 网格必须在板中心设置节点，以便直接施加 $100\ \mathrm{lb}$ 集中力。
- 若以向上为正，中心载荷和挠度均应为负；若程序采用相反符号，只比较绝对值。
- 此题 $L/t=200$，属于薄板极限测试。MITC4 或选择性积分元素应随网格细化向 $-0.0815\ \mathrm{in.}$ 收敛。
- $2\times2$ 书中结果相对经典解的幅值误差约为 $6.96\%$。
- 集中载荷点附近存在应力奇异性，本题只适合严格比较中心挠度，不适合严格比较载荷点最大应力。

### 1.6 建议验收记录

| 网格 | 完整积分 Q4 | 选择性积分 Q4 | MITC4 | 参考值 |
|---|---:|---:|---:|---:|
| $2\times2$ | 待计算 | 待计算 | 待计算 | 书中有限元：$-0.07583\ \mathrm{in.}$ |
| $4\times4$ | 待计算 | 待计算 | 待计算 | 应向经典解收敛 |
| $8\times8$ | 待计算 | 待计算 | 待计算 | 经典解：$-0.0815\ \mathrm{in.}$ |

## 2. Problem 12.1：均布载荷下四边固支方板

### 2.1 位置

- 题目及 Figure P12-1：`第12章/Problems.pdf` 文件第 1 页，书本第 591 页。
- 书后答案：`第12章/Answers.pdf` 文件第 1 页，书本第 929 页。

### 2.2 英文原文

> A square steel plate (Figure P12-1) of dimensions 20 in. × 20 in. with thickness of 0.1 in. is clamped all around. The plate is subjected to a uniformly distributed loading of 1 lb/in². Using a 2 × 2 mesh and then a 4 × 4 mesh, determine the maximum deflection and maximum stress in the plate. Compare the finite element solution to the classical one in [1].

### 2.3 中文题目

一块 $20\ \mathrm{in.}\times20\ \mathrm{in.}$ 的正方形钢板，厚度为 $0.1\ \mathrm{in.}$，四周完全固支。板面承受 $1\ \mathrm{lb/in^2}$ 的均布载荷。分别使用 $2\times2$ 和 $4\times4$ 网格，求板的最大挠度和最大应力，并将有限元解与参考文献 [1] 的经典解比较。

结构化输入：

| 项目 | 数值 |
|---|---:|
| 方板边长 | $20\ \mathrm{in.}\times20\ \mathrm{in.}$ |
| 板厚 $t$ | $0.1\ \mathrm{in.}$ |
| 均布载荷 $q$ | $1\ \mathrm{lb/in^2}$ |
| 边界 | 四边固支 |
| 要求网格 | $2\times2$、$4\times4$ |
| 输出 | 最大挠度、最大应力 |

### 2.4 书后答案

答案页给出的是 $8\times8$ 网格结果：

$$
\delta_{\max}=0.0785\ \mathrm{in.}
$$

$$
\sigma_{vM}=7046\ \mathrm{psi}
$$

书中注明：上述结果与解析解相符。

### 2.5 重要说明

- 题目要求计算 $2\times2$ 和 $4\times4$ 网格，但书后只列出 $8\times8$ 网格最终答案。
- 题面只写“steel”，没有在该页重列 $E$ 和 $\nu$。在自动验收前必须确认教材软件所采用的钢材参数。
- 未确认材料参数前，建议先比较无因次挠度或网格收敛趋势，不要把 $0.0785\ \mathrm{in.}$ 设成严格单值断言。
- 均布载荷必须转换为一致节点载荷；不建议只把总载荷平均分给节点。
- 应力结果依赖 Gauss 点外推、节点平滑和上下表面选择，应记录具体恢复方法。

### 2.6 建议验收记录

| 网格 | $w_{\max}$ | 相对 $0.0785\ \mathrm{in.}$ 的误差 | $\sigma_{vM,\max}$ | 备注 |
|---|---:|---:|---:|---|
| $2\times2$ | 待计算 | 待计算 | 待计算 | 粗网格 |
| $4\times4$ | 待计算 | 待计算 | 待计算 | 题目要求 |
| $8\times8$ | 待计算 | 待计算 | 待计算 | 书后答案：$0.0785\ \mathrm{in.}$、$7046\ \mathrm{psi}$ |

## 3. Problem 12.8：均布压力下四周固支圆板

### 3.1 位置

- 题目及 Figure P12-8：`第12章/Problems.pdf` 文件第 3 页，书本第 593 页。
- 书后答案：`第12章/Answers.pdf` 文件第 1 页，书本第 929 页。

### 3.2 英文原文

> Determine the maximum deflection and maximum principal stress of the circular plate shown in Figure P12-8. The plate is subjected to a uniform pressure $p=50$ kPa and fixed along its outer edge. Let $E=200$ GPa, $\nu=0.3$, radius $r=500$ mm, and thickness $t=20$ mm.

### 3.3 中文题目

求 Figure P12-8 所示圆板的最大挠度和最大主应力。圆板外缘完全固支，板面承受 $50\ \mathrm{kPa}$ 的均布压力。

结构化输入：

| 项目 | 数值 |
|---|---:|
| 半径 $r$ | $500\ \mathrm{mm}$ |
| 厚度 $t$ | $20\ \mathrm{mm}$ |
| 弹性模量 $E$ | $200\ \mathrm{GPa}$ |
| 泊松比 $\nu$ | $0.3$ |
| 均布压力 $p$ | $50\ \mathrm{kPa}$ |
| 边界 | 圆周完全固支 |
| 输出 | 最大挠度、最大主应力 |

### 3.4 书后答案

$$
\delta_{\max}=0.3306\ \mathrm{mm}
$$

$$
\sigma_{\max}=22.73\ \mathrm{MPa}
$$

题目要求的是最大主应力，因此后续验证应将书后 $\sigma_{\max}$ 按最大主应力解释，而不是 von Mises 应力。

### 3.5 独立解析核对

对于线性 Kirchhoff 理论下的四周固支圆板：

$$
D=\frac{Et^3}{12(1-\nu^2)}
=146520.15\ \mathrm{N\,m}
$$

中心挠度为：

$$
w_0=\frac{pr^4}{64D}
=0.33325\ \mathrm{mm}
$$

该解析值与书后 $0.3306\ \mathrm{mm}$ 相差约 $0.80\%$，因此题目数据和答案具有良好自洽性。

### 3.6 后续验证设置

- 若用完整圆板建模，圆周上施加 $w=\theta_x=\theta_y=0$。
- 若采用四分之一对称模型，应正确施加两条对称边的法向转角约束，并避免把对称边误设成固支边。
- 曲边用直边 Q4 逼近时应做圆周离散收敛；几何误差可能先于单元误差占主导。
- 优先比较中心挠度，再比较最大主应力。
- 应记录最大主应力出现的位置、板的上/下表面和应力恢复方式。
- 对精细网格，建议中心挠度相对书后答案误差小于 $2\%$；应力误差目标可先设为 $5\%$，再随网格细化收紧。

### 3.7 建议验收记录

| 网格/圆周分段 | $w_{\max}$ | 挠度误差 | $\sigma_{1,\max}$ | 应力误差 | 备注 |
|---|---:|---:|---:|---:|---|
| 粗网格 | 待计算 | 待计算 | 待计算 | 待计算 | 检查几何误差 |
| 中等网格 | 待计算 | 待计算 | 待计算 | 待计算 | 检查收敛 |
| 精细网格 | 待计算 | 待计算 | 待计算 | 待计算 | 参考：$0.3306\ \mathrm{mm}$、$22.73\ \mathrm{MPa}$ |

## 4. 推荐执行顺序

1. 先做 Example 12.1，只验证中心挠度和薄板极限。
2. 再做 Problem 12.8，验证曲边网格、均布载荷、固支边界和应力恢复。
3. 最后做 Problem 12.1，确认材料参数后进行规则方板的挠度与应力回归。

三个题目的书中答案汇总：

| 题目 | 书中有限元/答案 | 解析或补充基准 |
|---|---|---|
| Example 12.1 | $2\times2$：$w_{\max}=-0.07583\ \mathrm{in.}$ | 经典解：$-0.0815\ \mathrm{in.}$ |
| Problem 12.1 | $8\times8$：$\delta_{\max}=0.0785\ \mathrm{in.}$，$\sigma_{vM}=7046\ \mathrm{psi}$ | 书中注明与解析解相符；材料参数需确认 |
| Problem 12.8 | $\delta_{\max}=0.3306\ \mathrm{mm}$，$\sigma_{1,\max}=22.73\ \mathrm{MPa}$ | 固支圆板解析挠度：$0.33325\ \mathrm{mm}$ |

