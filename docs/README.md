# 项目文件地图与学习路线

这个目录按“先理解，再验证，最后看实现”的顺序使用即可，不需要逐文件阅读。

## 第一层：面试复习

| 文件 | 用途 |
| --- | --- |
| [`../START_HERE_博士面试.md`](../START_HERE_博士面试.md) | 总入口；先建立一分钟整体认识 |
| [`博士面试_数学核心速览.md`](博士面试_数学核心速览.md) | 理解数学链、两种板理论和数值难点 |
| [`博士面试_高频问答.md`](博士面试_高频问答.md) | 练习口头表达并检查是否真正理解 |

## 第二层：看完整问题如何走通

| 位置 | 用途 |
| --- | --- |
| [`../examples/`](../examples/) | 从连续体理论到完整验证的 9 个递进示例 |
| [`../output/pdf/`](../output/pdf/) | 两个完整算例及详细推导 |
| [`../Mindlin_Plate_教材题目与答案.md`](../Mindlin_Plate_教材题目与答案.md) | 教材题目、答案、验证注意事项 |

最值得面试前运行的是：

```bash
.venv/bin/python examples/step_09_complete_validation.py
```

它会一次展示材料刚度、Q4 映射、单元零模态、整体求解、剪切锁死、MITC4、结果恢复和完整验证。

## 第三层：需要追问时再看源码

| 数学问题 | 主要文件 |
| --- | --- |
| K/M 理论选择 | `mindlin_plate/theory.py` |
| 材料、弯曲刚度、剪切刚度 | `mindlin_plate/material.py` |
| 中面运动学 | `mindlin_plate/kinematics.py` |
| Q4 形函数、Jacobian、B 矩阵、MITC4 | `mindlin_plate/q4.py` |
| DKQ 薄板运动学 | `mindlin_plate/kirchhoff.py` |
| 单元刚度、荷载、单元响应 | `mindlin_plate/element.py` |
| 整体组装、约束和线性求解 | `mindlin_plate/assembly.py` |
| 斜边局部约束 | `mindlin_plate/boundary.py` |
| 曲率、内力和应力恢复 | `mindlin_plate/postprocess.py` |
| 数值验证闸门 | `mindlin_plate/validation.py` |
| API 和绘图 | `mindlin_plate/service.py`、`api.py`、`plotting.py` |

## 复习检查点

完成复习后，应能独立回答：

- 为什么板问题可以从三维降到中面？
- 为什么薄板和中厚板需要不同处理？
- 为什么普通 Q4 Mindlin 单元在薄板下会剪切锁死？
- DKQ、全积分、减缩积分和 MITC4 的差别是什么？
- `Ke`、`fe`、`K`、`f`、`u` 分别从哪里来？
- 支座条件如何进入方程？反力如何计算？
- 位移结果如何变成弯矩、剪力和表面应力？
- 哪些验证能说明结果可信？哪些能力本项目没有？

如果这些问题都能用自己的话讲清楚，就已经具备面试所需的整体认知；不需要背诵每个矩阵元素。
