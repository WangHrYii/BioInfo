# 第12章 单细胞转录组分析（scRNA-seq）

## 12.1 单细胞测序原理（10X Genomics 为例）

（此章节内容待完善）

## 12.2 数据预处理

- Cell Ranger 流程
- STARsolo（开源替代）
- 质控：过滤低质量细胞与双细胞（doublet，Scrublet/DoubletFinder）

## 12.3 标准分析流程（Seurat / Scanpy）

- 归一化（LogNormalize, SCTransform, scran）
- 特征选择（高变基因 HVGs）
- 降维：PCA → UMAP / t-SNE
- 聚类（Leiden / Louvain）
- 细胞类型注释（手动 marker vs 自动注释 SingleR / CellTypist）
- 差异表达基因（Marker genes）

## 12.4 多样本整合

- 批次效应问题
- Harmony
- scVI / scANVI
- Scanorama
- BBKNN
- 整合方法评估与选择

## 12.5 轨迹分析与拟时序

- Monocle3
- RNA Velocity（scVelo / velocyto）
- PAGA（图抽象分析）
- Palantir

## 12.6 细胞通讯分析

- CellChat
- CellPhoneDB
- NicheNet
- 配体-受体互作可视化

## 12.7 基因调控网络推断

- SCENIC / pySCENIC
- 转录因子活性评估
- 调控模块（regulon）鉴定

## 12.8 单细胞 CNV 推断

- inferCNV
- CopyKAT
- 肿瘤细胞 vs 正常细胞区分

## 12.9 单细胞多组学

- scATAC-seq（ArchR / Signac）
- CITE-seq（蛋白 + 转录组）
- Multiome（ATAC + RNA 联合分析）
- 单细胞多组学整合（WNN, MOFA+）
