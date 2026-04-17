# 生物信息学学习大纲

---

## 第一部分 生物学基础

### 第 1 章 从 DNA 到蛋白质

- [ ] 1.1 细胞是什么 — 生命的基本单位
- [ ] 1.2 DNA：生命的源代码 — 双螺旋、碱基配对、基因组/染色体/基因的关系
- [ ] 1.3 中心法则：DNA → RNA → 蛋白质 — 转录、翻译
- [ ] 1.4 基因的开关：基因表达调控 — 启动子、增强子、转录因子、表观遗传、非编码 RNA
- [ ] 1.5 突变与变异 — SNP、Indel、CNV、结构变异、突变后果
- [ ] 1.6 免疫系统基础 — 先天/适应性免疫、T/B/NK 细胞、肿瘤微环境
- [ ] 1.7 信号通路基础 — PI3K/AKT、MAPK、Wnt、Notch、JAK/STAT

### 第 2 章 测序技术

- [ ] 2.1 什么是测序
- [ ] 2.2 一代测序（Sanger）
- [ ] 2.3 二代测序（Illumina）— 边合成边测序、单端/双端、读长/深度/覆盖度
- [ ] 2.4 三代测序 — PacBio SMRT / Oxford Nanopore、长读长优势
- [ ] 2.5 常见测序类型一览 — WGS/WES、RNA-seq/scRNA-seq/空间、ChIP-seq/ATAC-seq/CUT&Tag、蛋白质组/代谢组
- [ ] 2.6 实验设计要点 — 生物学重复 vs 技术重复、测序深度选择、样本量

---

## 第二部分 编程基础

### 第 3 章 Linux / 命令行

- [X] ~~已掌握~~（日常开发环境，跳过）

### 第 4 章 Python 编程

- [X] ~~已掌握~~（主力语言，跳过）
- [ ] 4.7 Biopython 入门 — 序列操作、翻译、读取 GenBank/FASTA、调用 NCBI（**仅此小节需要学**）

### 第 5 章 R 语言

- [ ] 5.1 R 与 RStudio 安装
- [ ] 5.2 R 基础语法 — 向量、矩阵、data.frame、factor、list
- [ ] 5.3 数据操作 — dplyr (filter/select/mutate/group_by)、tidyr (pivot)
- [ ] 5.4 数据可视化 — ggplot2 图层语法、散点图/箱线图/热图/火山图、配色
- [ ] 5.5 Bioconductor 入门 — DESeq2、edgeR、clusterProfiler、GenomicRanges
- [ ] 5.6 R Markdown / Quarto — 可重复性分析报告

---

## 第三部分 核心生信分析

### 第 6 章 文件格式大全

- [ ] 6.1 序列文件：FASTA、FASTQ
- [ ] 6.2 比对文件：SAM、BAM、CRAM
- [ ] 6.3 变异文件：VCF、BED
- [ ] 6.4 注释文件：GFF、GTF
- [ ] 6.5 表达矩阵：counts matrix、TPM、FPKM、CPM
- [ ] 6.6 其他格式：BigWig、HDF5 (h5ad)、MEX (10X)、loom
- [ ] 6.7 格式转换实战

### 第 7 章 数据获取与质控

- [ ] 7.1 公共数据库 — GEO/SRA、ENCODE、TCGA/GTEx/ICGC、sra-tools/wget/aspera、GDC Portal
- [ ] 7.2 原始数据质控 — FastQC、MultiQC、Trimmomatic / fastp

### 第 8 章 序列比对

- [ ] 8.1 为什么要比对
- [ ] 8.2 参考基因组获取与索引构建
- [ ] 8.3 DNA 比对：BWA / BWA-MEM2、Bowtie2
- [ ] 8.4 RNA 比对：HISAT2、STAR（处理可变剪接）
- [ ] 8.5 比对后处理 — SAMtools 排序/索引/统计、Picard 去重、IGV 可视化

### 第 9 章 转录组分析（RNA-seq）

- [ ] 9.1 RNA-seq 分析完整流程概览
- [ ] 9.2 基因表达定量 — featureCounts/HTSeq、Salmon/Kallisto（无比对定量）、TPM vs FPKM vs Counts
- [ ] 9.3 差异表达分析 — DESeq2 实战、edgeR、limma-voom、火山图/MA 图、多组比较与批次效应校正
- [ ] 9.4 功能富集分析 — GO (BP/CC/MF)、KEGG、GSEA、GSVA、clusterProfiler、MSigDB
- [ ] 9.5 共表达网络 — WGCNA 原理与实战、模块-表型关联、Hub 基因
- [ ] 9.6 可变剪接分析 — rMATS、SUPPA2、剪接事件可视化
- [ ] 9.7 融合基因检测 — STAR-Fusion、Arriba
- [ ] 9.8 非编码 RNA 分析 — miRNA、lncRNA、circRNA (CIRI2/CIRCexplorer)、ceRNA 网络

### 第 10 章 基因组变异分析

- [ ] 10.1 Germline 变异检测 — GATK Best Practices、HaplotypeCaller、DeepVariant
- [ ] 10.2 Somatic 变异检测 — Mutect2、VarScan2、Strelka2
- [ ] 10.3 结构变异与 CNV — CNVkit、Manta/DELLY、GISTIC2
- [ ] 10.4 变异注释与过滤 — ANNOVAR、VEP、ClinVar/dbSNP/gnomAD、致病性预测 (SIFT/PolyPhen-2/CADD)
- [ ] 10.5 变异可视化 — maftools、突变特征 (SigProfiler)、TMB/MSI

### 第 11 章 表观基因组分析

- [ ] 11.1 ChIP-seq — Peak calling (MACS2/3)、差异结合 (DiffBind)、Motif (HOMER/MEME)、deepTools 可视化
- [ ] 11.2 ATAC-seq — 预处理（Tn5 偏好/线粒体去除）、开放染色质鉴定、核小体定位、TF 足迹、差异可及性
- [ ] 11.3 CUT&Tag / CUT&Run — 与 ChIP-seq 的区别、低背景噪音处理
- [ ] 11.4 DNA 甲基化 — Bismark 比对、DMS/DMR 鉴定、甲基化与表达关联
- [ ] 11.5 Hi-C / 3D 基因组 — TAD、A/B compartment、染色质环、HiC-Pro/Juicer/cooler

### 第 12 章 单细胞转录组（scRNA-seq）

- [ ] 12.1 单细胞测序原理（10X Genomics）
- [ ] 12.2 数据预处理 — Cell Ranger、STARsolo、质控/过滤双细胞 (Scrublet/DoubletFinder)
- [ ] 12.3 标准分析流程 — 归一化 (SCTransform/scran)、HVGs、PCA→UMAP/t-SNE、聚类 (Leiden)、细胞注释 (SingleR/CellTypist)、Marker genes
- [ ] 12.4 多样本整合 — Harmony、scVI/scANVI、Scanorama、BBKNN、整合方法评估
- [ ] 12.5 轨迹分析与拟时序 — Monocle3、RNA Velocity (scVelo)、PAGA、Palantir
- [ ] 12.6 细胞通讯 — CellChat、CellPhoneDB、NicheNet、配体-受体可视化
- [ ] 12.7 基因调控网络 — SCENIC/pySCENIC、TF 活性评估、regulon 鉴定
- [ ] 12.8 单细胞 CNV 推断 — inferCNV、CopyKAT、肿瘤 vs 正常区分
- [ ] 12.9 单细胞多组学 — scATAC-seq (ArchR/Signac)、CITE-seq、Multiome、WNN/MOFA+

### 第 13 章 空间转录组

- [ ] 13.1 技术原理 — Visium/Slide-seq/Stereo-seq (测序) vs MERFISH/seqFISH/CosMx/Xenium (成像)
- [ ] 13.2 数据分析流程 — Space Ranger、Scanpy+Squidpy、Seurat 空间模块、SPATA2
- [ ] 13.3 空间可变基因 — SpatialDE、SPARK、Moran's I
- [ ] 13.4 空间域识别 — BayesSpace、SpaGCN、STAGATE、stLearn
- [ ] 13.5 与 scRNA-seq 整合（反卷积）— Cell2location、RCTD、SPOTlight、Tangram
- [ ] 13.6 空间细胞通讯 — COMMOT、SpaTalk
- [ ] 13.7 空间轨迹与动态分析

---

## 第四部分 进阶分析技能

### 第 14 章 基因调控网络与虚拟基因敲除

- [ ] 14.1 GRN 基础 — 什么是基因调控网络、TF→靶基因调控逻辑
- [ ] 14.2 GRN 推断方法 — SCENIC/pySCENIC、GENIE3/GRNBoost2、dorothea、EGRN
- [ ] 14.3 蛋白活性推断 — VIPER 算法、decoupleR 框架、TF/通路活性推断实战
- [ ] 14.4 虚拟基因敲除 — CellOracle (GRN→扰动→细胞命运)、scGen (VAE)、GEARS (GNN)、CPA
- [ ] 14.5 Perturb-seq / CROP-seq 数据分析 — 实验原理、Pertpy 分析、计算预测 vs 实验验证
- [ ] 14.6 实战：CellOracle 虚拟 TF 敲除 — 数据准备→GRN 构建→模拟 KO→结果解读

### 第 15 章 肿瘤免疫与临床生信

- [ ] 15.1 肿瘤免疫微环境 — ESTIMATE、CIBERSORT/CIBERSORTx、xCell/MCP-counter/EPIC、TIMER2.0、免疫检查点基因
- [ ] 15.2 肿瘤新抗原预测 — HLA 分型 (OptiType)、NetMHCpan/pVACseq、免疫原性评估
- [ ] 15.3 生存分析 — Kaplan-Meier、Log-rank、Cox 回归、最佳截断值
- [ ] 15.4 预后模型构建 — LASSO Cox、风险评分、列线图 (Nomogram)、校准曲线、DCA、训练/验证集
- [ ] 15.5 分子分型 — ConsensusClusterPlus、NMF 分型、亚型间差异分析
- [ ] 15.6 药物敏感性 — GDSC/CCLE、pRRophetic/oncoPredict、药物-靶点网络

### 第 16 章 遗传学与群体遗传

- [ ] 16.1 GWAS — 原理与统计模型、PLINK、Manhattan/QQ 图、位点注释
- [ ] 16.2 孟德尔随机化（MR）— 工具变量选择、TwoSampleMR/MR-PRESSO、敏感性分析
- [ ] 16.3 多基因风险评分（PRS）— PRSice-2、LDpred2、临床应用
- [ ] 16.4 群体遗传学基础 — PCA 群体结构、连锁不平衡、选择压力

### 第 17 章 机器学习与深度学习

- [X] ~~17.1 监督学习基础~~ — 已掌握（RF/SVM/XGBoost/模型评估）
- [X] ~~17.2 无监督学习~~ — 已掌握（聚类/降维）
- [ ] 17.3 特征选择与生物标志物筛选 — LASSO/Elastic Net、Boruta、RFE（**生信特有应用**）
- [X] ~~17.4 深度学习入门~~ — 已掌握（CNN/RNN/Transformer/VAE/GNN）
- [ ] 17.5 生物大语言模型 — ESM/ProtTrans、DNABERT/Enformer、scGPT/Geneformer/scFoundation（**重点：可发挥 LLM 背景**）

### 第 18 章 多组学整合分析

- [ ] 18.1 为什么要多组学整合 — 单组学局限、早期/中期/晚期整合策略
- [ ] 18.2 转录组 + 表观基因组整合 — RNA-seq+ATAC-seq 关联、表达与甲基化关联、增强子-基因推断
- [ ] 18.3 多组学因子分析 — MOFA/MOFA+、iCluster/SNF、DIABLO (mixOmics)
- [ ] 18.4 单细胞多组学整合 — WNN、MOFA+、scAI/GLUE
- [ ] 18.5 蛋白质组与代谢组 — 差异蛋白/代谢物、MetaboAnalyst、与转录组整合

---

## 第五部分 工程化与可视化

### 第 19 章 工作流管理

- [ ] 19.1 为什么需要工作流管理（可重复性）
- [ ] 19.2 Snakemake 入门与实战（**Python 语法，可快速上手**）
- [ ] 19.3 Nextflow 入门与 nf-core 社区流程 — nf-core/rnaseq、nf-core/sarek
- [ ] 19.4 Docker / Singularity 环境打包
- [ ] 19.5 HPC 与云计算 — SLURM、AWS/GCloud、Terra/Galaxy

### 第 20 章 生信图表百科

> 40+ 种图表，可作工具书按需查阅，不必逐个学

- [ ] 20.1 出版级图表规范 — 分辨率/格式/字体/配色原则/多面板拼接
- [ ] 20.2 基础统计图 — 柱状图、箱线图、散点图、密度图、韦恩图
- [ ] 20.3 转录组图 — 火山图、MA 图、热图、PCA 图、富集气泡图/柱状图、GSEA 山峰图、WGCNA 模块图
- [ ] 20.4 基因组与变异图 — Manhattan 图、QQ 图、Oncoplot、突变特征谱、CNV 图、Lollipop 图
- [ ] 20.5 表观基因组图 — Peak 信号热图/Profile、基因组浏览器轨道、甲基化分布、Motif Logo
- [ ] 20.6 单细胞图 — UMAP/t-SNE、Dot Plot、Feature Plot、轨迹图、RNA Velocity、细胞通讯图、inferCNV 热图、SCENIC regulon 热图
- [ ] 20.7 空间转录组图 — 空间表达图、空间聚类图、反卷积饼图
- [ ] 20.8 临床与生存图 — KM 曲线、森林图、列线图、校准曲线、ROC、DCA
- [ ] 20.9 网络与关系图 — PPI 网络、GRN 图、ceRNA 网络
- [ ] 20.10 特殊图表 — Circos、UpSet、桑基图、弦图、瀑布图、河流图/鱼图、蜂巢图
- [ ] 20.11 虚拟扰动图 — 扰动向量场、蛋白活性热图

---

## 第六部分 论文复现实战

### 实战一 经典 RNA-seq 差异分析复现

- [ ] 原始数据下载 → 质控 → 比对 → 定量 → 差异分析 → 富集分析 → 图表复现
- [ ] 重点：读懂论文 Methods、找到对应数据
- [ ] 前置：第 7-9 章

### 实战二 TCGA 数据挖掘

- [ ] TCGA 数据获取 → 突变全景 → 生存分析 → 免疫浸润 → 预后模型 → 药物敏感性
- [ ] 重点：利用公共数据库做一篇生信文章
- [ ] 前置：第 9-10 章、第 15 章

### 实战三 单细胞转录组复现

- [ ] 数据获取 → 质控 → 聚类 → 注释 → 轨迹 → 细胞通讯 → GRN → 图表复现
- [ ] 重点：大规模单细胞数据处理、细胞类型注释
- [ ] 前置：第 12 章

### 实战四 空间转录组复现

- [ ] 数据获取 → 空间质控 → 空间聚类 → 空间可变基因 → scRNA 整合反卷积 → 空间通讯
- [ ] 重点：空间数据特殊处理、空间模式展示
- [ ] 前置：第 13 章

### 实战五 多组学整合复现

- [ ] 多组学数据获取 → 各组学独立分析 → 整合分析 → 调控关系推断
- [ ] 重点：不同组学数据关联讲完整生物学故事
- [ ] 前置：第 11 章、第 18 章

### 实战六 虚拟基因敲除与调控网络

- [ ] 单细胞数据 → GRN 构建 → 虚拟 KO → 细胞命运预测 → 实验验证对比
- [ ] 重点：用计算方法提出可验证的生物学假说
- [ ] 前置：第 12 章、第 14 章
