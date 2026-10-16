# 本地验证记录：Beam Vibration Lab 0.2.0

验证时间 UTC：2026-10-09T11:53:12.112389+00:00

环境：Python 3.12.14；Windows-11-10.0.26200-SP0。

本记录描述实际本地运行。GitHub Actions 已配置，尚未在远程执行。

## 执行结果

- 48 项测试通过，0 失败，0 错误；包括 9 项原有测试。
- 7 个算例脚本通过（6 个原算例，1 个扩展算例）。
- 静力、模态、谐响应、瞬态共 4 个 JSON 配置分析通过。
- 在仓库目录之外调用已安装的 beam-fem 命令通过，版本 0.2.0。
- src、examples、scripts、fenics 中的 Python 源码语法编译通过。
- Legacy FEniCS 三维求解未执行；语法检查不能替代后端运行验证。

## 新增实验结果

均布载荷 -20 N/m；钢梁 L=1 m，b=50 mm，h=10 mm，30 单元。

| 支承 | 最大挠度绝对值 / mm | 一阶频率 / Hz |
|---|---:|---:|
| 悬臂 | 3.000000 | 8.153807 |
| 简支 | 0.312500 | 22.888085 |
| 两端固支 | 0.062500 | 51.884738 |

释放 -1 N 静态预载、积分 0.5 s、2000 步：无阻尼机械能最大相对漂移 2.526e-10。
阻尼使一、二阶各具有 2% 阻尼比；最终机械能比 E(t)/E(0)=0.354212。
Rayleigh 系数 alpha_m=1.767273769 1/s，alpha_k=0.000107441180 s。

谐响应使用 1 N 端部谐波力、800 个频率点，检查零频静力极限、复数平衡残差与单模态共振解析响应。

## 可核查证据

- [运行清单](validation_evidence/manifest.json)
- [JUnit 测试结果](validation_evidence/pytest.xml)
- [逐步骤报告](validation_evidence/report.md)
- [扩展实验原始数值](../results/extended_summary.json)
- [谐响应频率、幅值和相位 CSV](../results/harmonic_response.csv)

## 复现环境

主要科学计算依赖复用了先前已验证环境中的相同版本，并复制到本项目独立虚拟环境；构建工具单独安装。
本地已验证 editable 安装。标准 wheel/source-distribution 构建与 wheel 安装检查随交付包另附 build_check.json。
精确依赖版本见 requirements-tested.txt；跨平台 CI 使用 pyproject.toml 中的兼容版本范围。

## 限制

采用小变形梁理论与密集矩阵；结果不代表三维实体孔洞应力、剪切效应或非线性设计校核。
原材料算例存在较大挠度，保留用于参数趋势展示。极细网格的频率误差会受舍入与条件数影响。
原仓库未指定许可证，本次不替原作者授权，也尚未创建或推送用户的远程仓库。

## 安装包复核

源码发行包与 wheel 构建通过；wheel 包含 NOTICE 来源声明。将 wheel 安装到另一虚拟环境后，在仓库外确认模块从 site-packages 载入；独立复制测试再运行，48 项通过。安装后 beam-fem 静力配置运行通过。

[构建与安装检查](validation_evidence/build_check.json) · [wheel 安装后测试 XML](validation_evidence/wheel-pytest.xml)
