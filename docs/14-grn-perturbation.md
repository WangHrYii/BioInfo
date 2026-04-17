# 第14章 基因调控网络与虚拟基因敲除

> 不做实验，用算法预测"敲掉一个基因会发生什么"

## 14.1 基因调控网络（GRN）基础

- 什么是 GRN？为什么要推断调控关系？
- 转录因子 → 靶基因的调控逻辑

## 14.2 GRN 推断方法

- SCENIC / pySCENIC：单细胞水平 GRN
- GENIE3 / GRNBoost2：基于随机森林的推断
- dorothea：转录因子-靶基因先验知识库
- EGRN（增强子驱动的 GRN）

## 14.3 蛋白活性推断

- VIPER 算法原理：从基因表达推断蛋白活性
- decoupleR：统一的活性推断框架
- 实战：推断转录因子活性、通路活性

## 14.4 虚拟基因敲除（in silico perturbation）

- 什么是虚拟基因敲除？与真实实验的关系
- CellOracle：基于 GRN 的虚拟 KO（预测细胞命运改变）
  - 构建 GRN → 模拟扰动 → 预测细胞状态变化
  - 可视化：扰动向量场
- scGen：基于变分自编码器（VAE）的扰动预测
- GEARS：基于图神经网络的组合扰动预测
- CPA（Compositional Perturbation Autoencoder）

## 14.5 Perturb-seq / CROP-seq 数据分析

- 实验原理：单细胞水平的基因扰动
- 数据分析流程（Pertpy）
- 计算预测 vs 实验验证对比

## 14.6 实战：用 CellOracle 进行虚拟转录因子敲除

- 数据准备
- GRN 构建
- 模拟 KO
- 结果解读与可视化
