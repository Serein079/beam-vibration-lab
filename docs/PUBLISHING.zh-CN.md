# 放到自己的 GitHub

## 先明确来源和许可

基线是 yanyana117/finite-element-vibration-simulation，commit `8eb430f`。
该基线没有 LICENSE，本版没有替原作者授权。请保留 NOTICE、来源链接和原作者署名。
如果这就是你自己的原仓库，可由你决定许可证；如果原仓库属于其他作者，独立公开再分发前应获得相应许可。
GitHub Fork 可以保留平台上的来源关系，但不能用它代替为整个项目指定新许可证所需的权利。
参考：[GitHub 官方许可证说明](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository)。

当前交付是本地改造版与发布包，没有替你创建远程仓库或推送。

## 独立仓库发布步骤（许可确认后）

建议名称 `beam-vibration-lab`，简介：
`Euler–Bernoulli beam FEM with multiple supports, harmonic response, energy validation and reproducible CLI workflows.`
建议 topics：`finite-element-method`、`beam`、`structural-dynamics`、`python`、`numerical-methods`。

先在你的 GitHub 创建空仓库，避免生成冲突的 README 和 LICENSE。
以下命令在项目目录执行，`YOUR_USERNAME` 替换为实际账号：

```powershell
cd D:\CAE_Projects\beam-vibration-lab
git status
git add .
git commit -m "Extend beam solver with supports, harmonic response and validation"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/beam-vibration-lab.git
git push -u origin main
```

本地已初始化时无需再次 `git init`。提交前配置你自己的 Git 用户名与邮箱。
如 origin 已存在，先执行 `git remote -v` 查看，再根据目标修改；不要直接推送到原作者仓库。
`git status --ignored` 可检查 `.venv`、缓存、构建产物和 runs 是否被忽略。

若选择在 GitHub 上 Fork 原仓库后继续开发，先克隆自己的 Fork，在新分支中应用本版改动，保留原 Git 历史与本版 NOTICE。
不要用发布包中的全新仓库元数据覆盖 Fork 历史。

## 验证与展示

```powershell
.\.venv\Scripts\python.exe -m pip install -e '.[dev]'
.\.venv\Scripts\python.exe scripts/run_all.py --output runs/prepublish
.\.venv\Scripts\python.exe -m build --sdist --wheel
```

推送后查看 Actions：配置会在 Windows/Linux、Python 3.11/3.12/3.13 上安装、运行测试和算例、构建包，并上传运行记录。
本地验证通过不等于所有远程矩阵已经通过。成功后可以建立 `v0.2.0` 发布，并附源码包、构建包和验证记录；许可证状态未解决前不要宣称项目采用某种开源许可证。

作品集介绍应聚焦你实际新增的内容：支承条件推广、支反力与力矩、频域分析、阻尼拟合、能量验证、参数校验、CLI 与工程化。
可用三张扩展图解释“边界条件如何影响刚度和频率”“为何共振幅值取决于阻尼”“自由振动能量如何变化”。
