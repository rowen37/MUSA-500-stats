# HW3 分工方案 —— 逻辑回归：酒驾车祸的预测因素

> 课程：CPLN 671 / MUSA 501　|　作业：Homework 3
> 数据：`assignment3/Logistic Regression Data.csv`（43,364 条事故，酒驾 2,485 条，占 5.73%）
> 说明：`assignment3/HW 3 - Instructions.pdf`
> 工具：R（需要 `gmodels`、`pROC` 包）
> 仓库：https://github.com/rowen37/MUSA-500-stats

## 总体流程

```
任务一（探索性分析）─表格/数据─▶ 任务二（建模与评价）─全部结果─▶ 任务三（论文式报告）
```

三人**串行接力**：后一人的输入是前一人的交付物。每完成一步就 commit + push。

---

## 任务一：数据准备 + 探索性分析

**负责人：＿＿＿＿**

对应说明第 1、2 步。脚本：`code/01_exploratory.R`

1. `setwd`；`install.packages` / `library` 加载 `gmodels` 等；`read.csv` 导入数据
2. `table` + `prop.table` 统计 `DRINKING_D` 的频数与占比
3. 对 7 个二元自变量（`FATAL_OR_M`、`OVERTURNED`、`CELL_PHONE`、`SPEEDING`、`AGGRESSIVE`、`DRIVER1617`、`DRIVER65PLUS`）：
   - `CrossTable` 交叉表（记录每个变量在 DRINKING_D=0/1 两类中取 1 的 N 和 %，以及总计）
   - 卡方检验：`CrossTable(x, mydata$DRINKING_D, prop.r=FALSE, prop.t=FALSE, prop.chisq=FALSE, chisq=TRUE)`，报告**不含 Yates 校正**的结果
4. 连续变量 `PCTBACHMOR`、`MEDHHINC`：
   - `tapply(x, mydata$DRINKING_D, mean)` 和 `tapply(x, mydata$DRINKING_D, sd)` 分组均值/标准差
   - 独立样本 t 检验：`t.test(x ~ mydata$DRINKING_D)`（记录 t 值、df、p 值）
5. 全部 9 个自变量的 Pearson 相关系数矩阵，判断多重共线性（|r| > 0.8 视为严重）

**交付物（放 `output/`）：**
- `table_01_dv_distribution.csv` —— 因变量频数/占比
- `table_02_crosstab_chisq.csv` —— 交叉表汇总（N、%、χ²、df、p 值），格式对照说明 p.4
- `table_03_continuous_ttest.csv` —— 均值、SD、t 值、df、p 值，对照说明 p.5
- `table_04_correlation_matrix.csv` —— 相关系数矩阵 + 多重共线性结论

**验收：** 每个二元变量都能报出"非酒驾中取1的N/%、酒驾中取1的N/%、合计"，卡方与 t 检验结论（p<0.05 与否）齐全。

---

## 任务二：逻辑回归建模 + 模型评价

**负责人：＿＿＿＿**　（依赖任务一的交付物）

对应说明第 3 步。脚本：`code/02_logistic_models.R`

1. **完整模型**（9 个自变量）：`glm(DRINKING_D ~ FATAL_OR_M + OVERTURNED + CELL_PHONE + SPEEDING + AGGRESSIVE + DRIVER1617 + DRIVER65PLUS + PCTBACHMOR + MEDHHINC, family = binomial, data = mydata)`
   - `summary` 输出
   - 每个变量算**优势比 `exp(β)` 和 95% CI `exp(confint())`**，按课件"Merging Odds Ratios"语法合并成"系数 + p 值 + OR + CI"矩阵
   - 阈值 0.02 / 0.03 / 0.05 / 0.07 / 0.08 / 0.09 / 0.10 / 0.15 / 0.20 / 0.50 下的**灵敏度、特异度、误分类率**，标出误分类率最低的阈值
   - **ROC 曲线**（`pROC`，导出 PNG），用"距左上角最小距离"法确定最优阈值，并算 **AUC**
2. **简化模型**（去掉 `PCTBACHMOR`、`MEDHHINC`，只用 7 个二元变量），同样输出系数/OR/CI 矩阵
3. 模型比较：`AIC(mylogit1, mylogit2)`，AIC 低者为优

**交付物（放 `output/`）：**
- `table_05_full_model_or.csv` —— 完整模型合并矩阵
- `table_06_cutoff_metrics.csv` —— 10 档阈值的灵敏度/特异度/误分类率（高亮最低误分类）
- `fig_01_roc_curve.png` —— ROC 曲线（图上标注最优阈值与 AUC）
- `table_07_reduced_model_or.csv` —— 简化模型合并矩阵
- `table_08_aic_comparison.csv` —— 两模型 AIC 对比

**验收：** 每个变量有 β、p、OR、CI；两个模型的显著性差异说清楚；最优阈值有两种口径（最低误分类 vs ROC 左上角距离）并对比。

---

## 任务三：撰写论文式报告

**负责人：＿＿＿＿**　（依赖任务一、二全部交付物）

对照说明最后的 REPORT OUTLINE，期刊格式：

1. **Introduction**（约 2 段）：酒驾车祸问题背景与费城背景；为什么这些变量可能与酒驾相关；说明用 R 做逻辑回归
2. **Methods**（3-4 页，只讲方法不放结果）：
   - OLS 不适用于二元因变量的原因；odds / odds ratio 概念与解读
   - 用**实际变量名**写出 logit 形式 `ln(p/(1-p)) = ...` 和 logistic 形式 `p = 1/(1+exp(-(...)))` 的回归方程，解释每一项、logit 与 logistic 函数性质
   - 每个变量的假设检验（H0/Ha）、Wald 统计量及其分布；为什么看 OR 而不是 β
   - 模型评价：伪 R² 的局限、AIC 定义（越低越好）、灵敏度/特异度/误分类率怎么算、ŷ 如何得到、为什么试不同阈值、ROC 曲线与最优阈值方法、AUC 分级标准（excellent/good/fair/poor/failing）
   - 逻辑回归假设（哪些 OLS 假设保留、哪些不适用）
   - 探索性分析（卡方检验、t 检验）的方法说明与假设
3. **Results**（3-5 页）：放入任务一、二的全部表格与图，逐一解读——因变量分布、交叉表+卡方结论、t 检验结论、相关矩阵与多重共线性、假设满足情况；完整模型结果（逐变量显著性 + OR 解读）；阈值表与最优阈值对比；ROC 与 AUC 解读；简化模型对比与 AIC 结论
4. **Discussion**（1-2 段）：哪些是强/弱预测因子、方向是否符合预期、**是否应改用 Paul Allison 的稀有事件方法**（关键线索：因变量中 1 仅占 5.73%）、局限性

**交付物：** `output/HW3_report.pdf`（及 LaTeX/Word 源文件），格式参照 assignment1/2 的提交包结构，最后更新 `README.txt`。

---

## 协作约定

- **文件结构**（与 HW1/HW2 一致）：`code/` 放脚本，`output/` 放表格和图，数据不动
- **命名**：脚本按 `01_`、`02_` 编号；表格 `table_XX_说明.csv`；图 `fig_XX_说明.png`
- **读数据注意**：文件名含空格，R 里写 `read.csv("Logistic Regression Data.csv")`
- **进度**：每完成一个任务立即 `git add` + `commit` + `push`，commit message 注明任务编号（如 `HW3 task1: exploratory analysis tables`）
- 截止日期：＿＿＿＿
