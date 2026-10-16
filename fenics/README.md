# FEniCS 三维连续体示例

**本版状态：未执行三维求解，不能使用 1D 测试结果替代本目录验证。**
这些脚本使用 legacy FEniCS/dolfin/mshr/ufl，不适用于直接替换成 FEniCSx。
安装说明参考 https://fenicsproject.org/download/archive/ 。

本版修改：

- `hollow_box_beam_static.py` 使用单一中心空腔构造闭口箱形截面（壁厚 50.8 mm）；原来的两个侧空腔形成开口截面，且一处坐标反向。
- 重力体力为 `rho*g` [N/m³]，在体积分中不乘总体积。
- 顶面局部面力方向定义为负 z，与重力同向。这是本版明确采用的载荷工况。
- 瞬态脚本先输出静态初始形状 `t=0`，然后输出每一步实际时间；历史位移由释放后的初始加速度作 Taylor 初始化。
- 三维瞬态采用原有隐式位移差分式，可能有数值耗散；不是 1D 求解器的 Newmark 平均加速度实现。

运行时应先检查网格截面、总重力与外载积分、支反力平衡，并做网格/时间步收敛。
`PROJECT_SUMMARY.md` 保留原背景说明，其中对原实现的描述不覆盖上述修改。

本目录的脚本与仓库主体是**两套独立实现**，不共享代码：

| | 实现 | 依赖 |
|---|---|---|
| 仓库主体 `src/fem_beam.py` | 一维梁理论 | NumPy、Matplotlib |
| 本目录 | 三维连续体有限元 | FEniCS / dolfin、mshr、ufl |

因此单独放置，避免与 `examples/` 下基于 `fem_beam` 的脚本混淆。
运行本目录脚本需另行安装 FEniCS。

| 文件 | 内容 |
|---|---|
| `static_cantilever_deflection.py` | 三维悬臂梁在分布面力下的静态挠度 |
| `hollow_box_beam_static.py` | 空心箱梁在自重与面载下的静态变形 |
| `transient_beam_vibration.py` | 静态预载后的悬臂梁瞬态振动 |
| `PROJECT_SUMMARY.md` | 原项目说明 |
