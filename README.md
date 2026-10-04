# magnetloss-calculation

一个用于**脉冲高场螺线管磁场与损耗计算**的纯 Python 开源程序包。

**作者：hery rainforeset**

本项目最初围绕 APEX Phase-I 脉冲高场磁体开展，但代码已经按物理层级拆分，可独立用于：

- 轴对称多匝螺线管磁场计算；
- 单根圆丝的有限磁扩散 / 非完全穿透计算；
- 1 丝、19 丝、304 丝二维显式有限元计算；
- 二级子缆损耗 surrogate；
- 整磁体逐匝损耗映射；
- 单脉冲绝热温升反馈计算。

整个仓库均为 **Python**，不依赖 MATLAB、COMSOL、GetDP 或其他商业有限元软件。

---

## 1. 当前参考磁体

默认参考模型为 APEX Phase-I 脉冲螺线管：

| 参数 | 数值 |
|---|---:|
| Double Pancake 数量 | 9 |
| Single Pancake 数量 | 18 |
| 每个 SP 径向匝数 | 9 |
| 总匝数 | 162 |
| 线圈内半径 | 125 mm |
| 线圈外半径 | 215.7 mm |
| 有效轴向长度 | 332.7 mm |
| 电缆包络尺寸 | 8.3 mm × 16.9 mm |
| 峰值电流 | 47.29 kA |
| 电流波形 | 3 ms 上升 / 1 ms 平顶 / 3 ms 下降 |
| 单丝直径 | 0.2 mm |
| 每匝总丝数 | 2432 |
| 当前有效 RRR 模型 | 3000 |

当前多级导体拓扑为：

```text
0.2 mm 单丝
    -> 19 丝一级子缆
    -> 16 根一级子缆 = 304 丝二级子缆
    -> 8 根二级子缆 = 2432 丝三级总缆
```

默认假设所有细丝之间**电绝缘**，因此当前版本不包含丝间耦合电流损耗。

---

## 2. 计算模型层级

### 2.1 有限矩形截面 Biot–Savart 磁场模型

每一匝宏观导体都不是简化成零截面圆环，而是视为一个具有有限径向和轴向尺寸的矩形电流密度区域。

程序使用 Gauss–Legendre 积分，将一个矩形截面匝离散成大量微小圆环，并对圆环磁场进行数值积分。

这样可以避免在导体自身位置使用零截面圆环模型时出现的自场奇异性。

主要代码：

```text
src/apexloss/field.py
```

当前参考结果：

```text
Ipk = 47.29 kA

中心场 B0                 ≈ 20.19523 T
最大匝中心磁场 |B|max     ≈ 20.63009 T
```

---

### 2.2 单丝有限磁扩散 / 非完全穿透模型

`SingleStrandDiffusion` 用于计算圆形金属细丝在横向脉冲磁场中的有限磁扩散过程。

该模型求解横向 `m = 1` 模式的时间域磁扩散问题，能够描述：

- 磁场向细丝内部扩散；
- 非完全穿透效应；
- 脉冲结束后的磁扩散尾项；
- 动态电阻率对涡流损耗的影响。

主要代码：

```text
src/apexloss/strand_diffusion.py
```

它可以与完全穿透近似直接比较：

```math
Q_{\rm CP}
=
\int
\frac{\pi a^4}{4\rho(B,T)}
\left(\frac{dB}{dt}\right)^2 dt
```

当前冻结的单丝回归算例：

```text
单丝直径                 = 0.2 mm
峰值磁场                 = 20.619192 T
波形                     = 3/1/3 ms
径向网格                 = 160
时间步长                 = 1 us

Q(0–7 ms)                = 0.5810908 J/m
脉冲后扩散尾项           = 0.00723958 J/m
有限扩散总损耗           = 0.5883304 J/m
完全穿透参考值           = 0.6083368 J/m
```

原始独立 Python 脚本保留在：

```text
legacy/APEX_strand_dynamic_rho_diffusion.py
```

用于结果追溯与回归验证。

---

### 2.3 1 / 19 / 304 丝二维显式 FEM

局部多丝模型采用二维标量磁矢势 `A_z` 形式，并使用 P1 三角形线性有限元。

当前内置三种几何：

- 单根 0.2 mm 圆丝；
- 19 丝 `1+6+12` 一级子缆；
- 304 丝二级子缆。

主要代码：

```text
src/apexloss/multifilament_fem.py
```

对于 field-only 工况，每根丝都独立施加零净输运电流约束，从而只保留单丝内部涡流。

当前模型的关键假设是：

```math
Q_{\rm coupling} \approx 0
```

因此总损耗近似写为：

```math
Q_{\rm total}
\approx
Q_{\rm transport}
+
Q_{\rm intra-filament}
```

当前二维 FEM 使用嵌入式结构三角网格，并非严格贴体网格，因此圆丝边界仍存在一定几何离散误差。

---

### 2.4 304 丝二级子缆 surrogate

为了避免对整个磁体中数十万根细丝进行直接有限元建模，程序将 304 丝二级子缆的完整 FEM 结果参数化为局部损耗 surrogate。

当前数据库按以下变量建立：

```text
Bpk
T
theta_B
```

其中方向相关性表示为：

```math
Q(\theta)
=
Q_{xx}\cos^2\theta
+
Q_{yy}\sin^2\theta
+
2Q_{xy}\sin\theta\cos\theta
```

低场区域不直接对 `Q(B)` 线性插值，而是对 `Q/B²` 插值，从而保持：

```math
Q \propto B^2
```

的低场物理标度。

主要代码：

```text
src/apexloss/surrogate.py
```

---

## 3. 整磁体损耗映射

整磁体模型采用：

```text
162 匝
× 每匝 2×2 Gauss 点
= 648 个宏观局部磁场点
```

整个计算链为：

```text
162 匝有限截面 Biot–Savart
        ↓
每个 Gauss 点得到 Br、Bz
        ↓
304 丝二级子缆 surrogate
        ↓
× 8 根二级子缆
        ↓
逐匝积分
        ↓
整磁体损耗
```

对于当前几何，2×2、3×3 和 4×4 截面 Gauss 采样在总损耗上已经基本收敛，因此默认采用 2×2。

主要代码：

```text
src/apexloss/whole_magnet.py
```

---

## 4. 热反馈模型

当前提供单脉冲绝热热反馈模型。

在脉冲过程中：

```math
T \uparrow
\Rightarrow
\rho(T,B) \uparrow
```

这会产生两个相反效应：

```text
温度升高 -> 电阻率升高 -> 涡流降低
温度升高 -> 电阻率升高 -> 输运焦耳热增加
```

因此实际单脉冲损耗不能简单使用全程 4.2 K 的等温损耗。

当前热模型不包含：

- 冷却剂；
- 匝间传热；
- 绝缘层热阻；
- 结构材料热容；
- 重复脉冲热恢复。

---

## 5. 当前参考结果

以下数值用于当前版本的回归验证，并不代表已经完成最终工程误差评估。

| 物理量 | 当前参考值 |
|---|---:|
| 中心场 | 20.19523 T |
| 最大匝中心磁场 | 20.63009 T |
| 固定 4.2 K 整磁体本征细丝涡流损耗 | 约 73.4 kJ/pulse |
| 旧 648 点单丝 spatial-RVE 热反馈涡流损耗 | 18.8466 kJ/pulse |
| 304 丝 reduced-order 热反馈涡流损耗 | 约 19.9–20.3 kJ/pulse |
| 当前输运焦耳热 | 约 2.3 kJ/pulse |
| 脉冲结束平均温度 | 约 35.8 K |
| 脉冲结束局部最高温度 | 约 47 K |

不同局部模型之间的差异被保留为**模型形式不确定度**，而不是通过人为调参消除。

---

## 6. 安装

建议 Python 3.10 或更高版本。

```bash
git clone https://github.com/zxcvbnm22003Hyl/magnetloss-calculation.git
cd magnetloss-calculation

python -m pip install -e .
```

如果需要运行测试：

```bash
python -m pip install -e .[dev]
pytest -q
```

---

## 7. 命令行使用

### 7.1 计算 162 匝磁场

```bash
apexloss field --out turn_fields.csv
```

---

### 7.2 单丝有限磁扩散

使用当前材料表：

```bash
apexloss strand \
  --Bpk 20.619192 \
  --temperature 4.2 \
  --nr 160 \
  --dt-us 1
```

使用历史回归算例中的 4.2 K Kohler / RRR 模型：

```bash
apexloss strand \
  --Bpk 20.619192 \
  --nr 160 \
  --dt-us 1 \
  --legacy-kohler
```

---

### 7.3 304 丝二维 FEM

```bash
apexloss secondary-fem \
  --Bpk 20.619192 \
  --temperature 4.2 \
  --h-mm 0.03 \
  --dt-ms 0.25
```

显式 304 丝 FEM 的计算量显著高于 surrogate 查询。

---

### 7.4 整磁体固定温度损耗

```bash
apexloss whole-magnet \
  --mode fixed \
  --temperature 4.2 \
  --target-order 2 \
  --outdir output_fixed
```

---

### 7.5 整磁体绝热热反馈

```bash
apexloss whole-magnet \
  --mode thermal \
  --temperature 4.2 \
  --target-order 2 \
  --dt-us 20 \
  --outdir output_thermal
```

---

## 8. Python API 示例

### 8.1 磁场

```python
from apexloss import MagnetGeometry, FiniteTurnBiotSavart

field = FiniteTurnBiotSavart(
    MagnetGeometry(),
    current_A=47.29e3,
    source_order=20,
)

print(field.center_field_T())

turns = field.turn_center_fields()
print(turns.head())
```

---

### 8.2 单丝有限穿透

```python
from apexloss import SingleStrandDiffusion
from apexloss.materials import LegacyKohlerRRR

model = SingleStrandDiffusion(
    diameter_m=0.2e-3,
    radial_cells=160,
)

result = model.simulate_trapezoid(
    B_peak_T=20.619192,
    dt_s=1e-6,
    tail_s=23e-3,
    resistivity=LegacyKohlerRRR(),
)

print(result.Q_total_J_per_m)
```

---

### 8.3 304 丝显式 FEM

```python
from apexloss.multifilament_fem import MultifilamentFEM

fem = MultifilamentFEM.secondary304(
    mesh_h_m=0.03e-3
)

result = fem.run_trapezoid(
    B_peak_T=20.619192,
    temperature_K=4.2,
    dt_s=0.25e-3,
)

print(result.energy_J_per_m)
```

---

## 9. 仓库结构

```text
src/apexloss/
    field.py
        有限矩形截面 Biot–Savart

    strand_diffusion.py
        单丝有限磁扩散 / 非完全穿透

    multifilament_fem.py
        1 / 19 / 304 丝二维显式 FEM

    surrogate.py
        304 丝损耗 surrogate

    whole_magnet.py
        162 匝整磁体损耗映射

    materials.py
        rho(T,B)、热容、焓

    geometry.py
        参考磁体几何

    waveforms.py
        脉冲波形

    data/
        材料与 surrogate 数据

examples/
    可直接运行的示例

legacy/
    历史 Python 模型，用于结果追溯

reference_results/
    当前冻结的数值回归结果

tests/
    单元测试与回归测试

configs/
    默认参考参数

docs/
    数值方法、模型假设和验证说明
```

---

## 10. 数值验证

可以运行：

```bash
python scripts/validate_reference.py
```

用于验证当前代码是否能够复现冻结基准。

仓库同时配置了 GitHub Actions，用于自动运行测试。

---

## 11. 当前模型的主要假设与限制

当前版本需要特别注意以下边界：

- 默认细丝之间完全电绝缘；
- 不包含丝间耦合电流损耗；
- 不包含二级/三级子缆之间的有限接触电阻；
- v0.1.0 暂不包含真实纽绞几何；
- 暂不包含 Helicoidal Transformation；
- 不包含结构金属涡流；
- 不包含接头和引线损耗；
- 热模型仅针对单脉冲绝热过程；
- 当前高纯 Al `rho(T,B)` 仍属于设计级材料模型；
- 最终工程计算应使用实际拉丝、绝缘和成缆后的实测材料参数；
- 304 丝 FEM 使用非贴体嵌入式网格；
- reduced-order 热 surrogate 不能完全替代具有内部磁扩散状态记忆的直接瞬态 FEM。

进一步说明见：

- [数值方法](docs/numerical_methods.md)
- [模型假设](docs/model_assumptions.md)
- [验证基准](docs/validation.md)
- [中文简要说明](docs/README_zh-CN.md)

---

## 12. 后续计划

计划继续加入：

- 更高精度的贴体网格 FEM；
- 有限丝间接触电阻；
- inter-strand coupling loss；
- 多级纽绞结构；
- Helicoidal Transformation；
- 绞距参数扫描；
- measured `rho(B,T)` 数据接口；
- 重复脉冲热恢复与冷却模型；
- 更一般的脉冲波形；
- 更通用的磁体几何输入。

---

## 13. 引用

如果该程序用于论文、报告或其他学术工作，请引用本 GitHub 仓库对应版本；后续如有相关方法论文，可同时引用论文。

---

## 14. License

BSD-3-Clause。

Copyright © 2026 **hery rainforeset**.
