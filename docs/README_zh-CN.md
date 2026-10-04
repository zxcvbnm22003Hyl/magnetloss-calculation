# APEX-PulseLoss 中文说明

这是一个纯 Python 的脉冲高场磁体损耗计算仓库，主要包含：

- 162 匝有限矩形截面 Biot-Savart 磁场计算；
- 0.2 mm 圆丝单丝有限磁扩散/非完全穿透数值模型；
- 1 丝、19 丝和 304 丝二维显式 FEM；
- 304 丝二级子缆损耗 surrogate；
- 162 匝、每匝 2x2 Gauss 点的整磁体损耗映射；
- 单脉冲绝热温升反馈。

默认假设芯丝之间电绝缘，因此不包含丝间耦合电流损耗。v0.1.0 也暂不包含真实纽绞或 Helicoidal Transformation。

原始单丝有限穿透脚本和旧的 648 点 spatial-RVE Python 程序保留在 `legacy/`，用于结果追溯。
