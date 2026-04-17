# 第11章 表观基因组分析

## 11.1 ChIP-seq 分析流程

- Peak calling（MACS2 / MACS3）
- 差异结合分析（DiffBind）
- Motif 分析（HOMER, MEME）
- 可视化（deepTools：heatmap, profile plot）

## 11.2 ATAC-seq 分析流程

- 数据预处理特殊注意点（Tn5 偏好、线粒体去除）
- 开放染色质区域鉴定
- 核小体定位分析
- 转录因子足迹分析（Footprinting）
- 差异可及性分析

## 11.3 CUT&Tag / CUT&Run

- 与 ChIP-seq 的区别
- 数据分析流程（低背景噪音处理）

## 11.4 DNA 甲基化分析

- Bismark 比对与甲基化提取
- 差异甲基化位点（DMS）与区域（DMR）鉴定
- 甲基化与基因表达的关联分析

## 11.5 Hi-C / 3D 基因组分析

- TAD（拓扑关联域）鉴定
- A/B compartment 分析
- 染色质环（Loop）检测
- 工具：HiC-Pro, Juicer, cooler
