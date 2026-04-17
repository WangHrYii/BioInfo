# 技能树导图 —— 一张图读懂生物信息学全貌

> 在开始学习之前，先用几张技能树看清楚：生信到底在干什么？每个技术解决什么问题？它们之间有什么关系？
> 不用背下来，随时回来查阅就好。

---

## 技能树总览

生物信息学 = **生物知识** + **编程能力** + **数据**。三者缺一不可。

```mermaid
graph TD
    ROOT["生物信息学 Bioinformatics"] --> BIO["生物知识"]
    ROOT --> CODE["编程能力"]
    ROOT --> DATA["数据来源"]

    BIO --> BIO1["DNA / RNA / 蛋白质"]
    BIO --> BIO2["中心法则"]
    BIO --> BIO3["基因调控"]
    BIO --> BIO4["突变与变异"]
    BIO --> BIO5["测序技术原理"]

    CODE --> CODE1["Linux 命令行"]
    CODE --> CODE2["Python"]
    CODE --> CODE3["R 语言"]
    CODE --> CODE4["Shell 脚本"]
    CODE --> CODE5["Conda 环境管理"]

    DATA --> DATA1["公共数据库"]
    DATA --> DATA2["自己测序的数据"]
    DATA1 --> DATA3["GEO / SRA"]
    DATA1 --> DATA4["TCGA / GTEx"]
    DATA1 --> DATA5["ENCODE"]
```

**怎么理解这张图？**

- **左边**是你需要的生物学知识——不需要背教科书，但需要理解"数据代表什么意思"
- **中间**是你需要的编程工具——生信分析全靠命令行和代码
- **右边**是数据来源——你可以用公共数据库里的数据做分析，很多论文就是靠挖掘公共数据发表的

---

## 核心技术分支：我测了数据，能分析什么？

这是最重要的一张图。根据你手上的数据类型，对应不同的分析方向。

### DNA 数据（基因组）

```mermaid
graph LR
    DNA["DNA 数据"] --> WGS["WGS / WES"]
    DNA --> BS["BS-seq / WGBS"]

    WGS --> VAR["变异检测"]
    VAR --> VAR1["遗传病致病变异筛选"]
    VAR --> VAR2["肿瘤体细胞突变分析"]
    VAR --> VAR3["突变特征 Mutational Signatures"]
    VAR --> VAR4["TMB / MSI"]

    BS --> METH["DNA 甲基化分析"]
    METH --> METH1["甲基化水平定量"]
    METH --> METH2["差异甲基化区域 DMR"]
```

- **WGS**（全基因组测序）/ **WES**（全外显子测序）：找 DNA 上的变异（SNP、Indel、CNV、结构变异）
- **BS-seq / WGBS**：看 DNA 上哪些位置被甲基化了（表观遗传调控）

### RNA 数据（转录组）

```mermaid
graph LR
    RNA["RNA 数据"] --> BULK["Bulk RNA-seq"]
    RNA --> SC["scRNA-seq"]
    RNA --> SP["空间转录组"]

    BULK --> B1["差异表达 DEGs"]
    BULK --> B2["功能富集 GO/KEGG/GSEA"]
    BULK --> B3["共表达网络 WGCNA"]
    BULK --> B4["可变剪接 / 融合基因"]
    BULK --> B5["免疫浸润分析"]

    SC --> S1["聚类与细胞注释"]
    SC --> S2["轨迹分析 / RNA Velocity"]
    SC --> S3["细胞通讯 CellChat"]
    SC --> S4["基因调控网络 SCENIC"]
    SC --> S5["CNV 推断 inferCNV"]

    SP --> SP1["空间可变基因"]
    SP --> SP2["空间域识别"]
    SP --> SP3["空间细胞通讯"]
    SP --> SP4["与 scRNA-seq 整合"]
```

转录组是生信分析中**论文最多**的方向：
- **Bulk RNA-seq**：最成熟，一团细胞的平均表达
- **scRNA-seq**：近年最火，单个细胞分辨率
- **空间转录组**：最前沿，保留组织空间信息

### 表观基因组数据

```mermaid
graph LR
    EPI["表观基因组数据"] --> CHIP["ChIP-seq"]
    EPI --> ATAC["ATAC-seq"]
    EPI --> CUT["CUT&Tag / CUT&Run"]

    CHIP --> C1["Peak calling"]
    CHIP --> C2["差异结合分析"]
    CHIP --> C3["Motif 分析"]

    ATAC --> A1["开放染色质鉴定"]
    ATAC --> A2["TF 足迹分析"]
    ATAC --> A3["scATAC-seq"]

    CUT --> CUT1["低细胞量表观分析"]
```

### 蛋白质组 / 代谢组 / 多组学

```mermaid
graph LR
    OTHER["其他组学"] --> PROT["蛋白质组"]
    OTHER --> META["代谢组"]
    OTHER --> MULTI["多组学整合"]

    PROT --> P1["差异蛋白分析"]
    PROT --> P2["通路富集"]

    META --> M1["差异代谢物"]
    META --> M2["代谢通路分析"]

    MULTI --> MU1["转录组 + 表观基因组"]
    MULTI --> MU2["多组学因子分析 MOFA"]
    MULTI --> MU3["单细胞多组学"]
```

---

## 临床/转化分析技能树

如果你的研究方向偏临床，这些是你要重点关注的分析：

### 肿瘤生信

```mermaid
graph TD
    TUMOR["肿瘤生信"] --> T1["TCGA / ICGC 数据挖掘"]
    TUMOR --> T2["生存分析"]
    TUMOR --> T3["肿瘤免疫微环境"]
    TUMOR --> T4["免疫浸润定量"]
    TUMOR --> T5["新抗原预测"]
    TUMOR --> T6["药物敏感性预测"]
    TUMOR --> T7["分子分型"]
    TUMOR --> T8["预后模型"]

    T2 --> T2A["Kaplan-Meier"]
    T2 --> T2B["Cox 回归"]
    T4 --> T4A["CIBERSORT"]
    T4 --> T4B["TIMER / xCell"]
    T8 --> T8A["LASSO"]
    T8 --> T8B["列线图 Nomogram"]
```

> 肿瘤生信是目前**发文章数量最多**的方向，很多"生信文章"本质上就是 TCGA 数据挖掘 + 免疫浸润 + 生存分析 + 预后模型。

### 遗传学与药物基因组学

```mermaid
graph TD
    GENE["遗传学分析"] --> G1["GWAS"]
    GENE --> G2["孟德尔随机化 MR"]
    GENE --> G3["PRS 多基因风险评分"]

    DRUG["药物基因组学"] --> D1["药物靶点预测"]
    DRUG --> D2["药物重定位"]
    DRUG --> D3["虚拟筛选"]
```

---

## 基因调控网络与虚拟扰动技能树

这是一个比较前沿的方向——**不做实验，用算法预测"敲掉一个基因会发生什么"**。

```mermaid
graph TD
    GRN["基因调控网络 GRN"] --> GRN1["SCENIC / pySCENIC"]
    GRN --> GRN2["GENIE3 / GRNBoost2"]
    GRN --> GRN3["dorothea"]
    GRN --> GRN4["CellOracle"]

    PERT["虚拟基因敲除"] --> P1["VIPER 蛋白活性推断"]
    PERT --> P2["CellOracle 虚拟 KO"]
    PERT --> P3["scGen 深度学习扰动"]
    PERT --> P4["GEARS 图神经网络"]
    PERT --> P5["Perturb-seq 数据验证"]

    CAUSAL["因果推断"] --> CA1["因果网络推断"]
    CAUSAL --> CA2["孟德尔随机化"]
    CAUSAL --> CA3["贝叶斯网络"]
```

核心思路：
1. 先从数据中**推断出基因之间的调控关系**（GRN）
2. 然后在计算机里**模拟敲除某个基因**，看细胞状态会怎么变
3. 最后与真实实验（Perturb-seq / CRISPR）对比验证

---

## 进阶与工程化技能树

当你掌握了基础分析后，这些进阶能力会让你更强：

```mermaid
graph TD
    ADV["进阶能力"] --> ML["机器学习 / 深度学习"]
    ADV --> VIZ["高级可视化"]
    ADV --> ENG["工程化"]

    ML --> ML1["随机森林 / SVM / XGBoost"]
    ML --> ML2["CNN / RNN / Transformer"]
    ML --> ML3["图神经网络 GNN"]
    ML --> ML4["生物大语言模型"]
    ML4 --> ML4A["蛋白质: ESM / ProtTrans"]
    ML4 --> ML4B["DNA: DNABERT / Enformer"]
    ML4 --> ML4C["单细胞: scGPT / Geneformer"]

    VIZ --> VIZ1["ComplexHeatmap"]
    VIZ --> VIZ2["Circos 圈图"]
    VIZ --> VIZ3["交互式: Plotly / Shiny"]
    VIZ --> VIZ4["基因组浏览器 IGV"]

    ENG --> ENG1["工作流: Snakemake / Nextflow"]
    ENG --> ENG2["容器: Docker / Singularity"]
    ENG --> ENG3["HPC / SLURM"]
    ENG --> ENG4["云计算 AWS / GCloud"]
```

---

## 技能树阅读指南：推荐学习路线

不知道从哪开始？按这个顺序走：

```mermaid
graph TD
    S1["第1站: 生物学基础"] --> S2["第2站: 编程基础"]
    S2 --> S3["第3站: 通用技能"]
    S3 --> S4["第4站: RNA-seq 分析"]
    S4 --> S5["第5站: 选一个方向深入"]
    S5 --> S6["第6站: 高级分析"]
    S6 --> S7["第7站: 论文复现实战"]

    S1 -.- S1D["第1-2章 不懂生物就是在盲分析"]
    S2 -.- S2D["第3-5章 磨好工具再上路"]
    S3 -.- S3D["第6-8章 文件格式、质控、比对"]
    S4 -.- S4D["第9章 最常用，论文数量最多"]
    S6 -.- S6D["第14-17章 调控网络、临床分析、ML"]
    S7 -.- S7D["第21-26章 学以致用"]
```

**第5站：选一个方向深入**

| 方向 | 对应章节 | 热度 |
|------|---------|------|
| 肿瘤方向 | 第10章（变异）+ 第15章（肿瘤免疫） | ★★★★★ 发文最多 |
| 单细胞方向 | 第12章（scRNA-seq） | ★★★★★ 最火技术 |
| 空间组学 | 第13章（空间转录组） | ★★★★ 最前沿 |
| 表观方向 | 第11章（ChIP/ATAC-seq） | ★★★★ |
| 遗传学方向 | 第16章（GWAS / MR / PRS） | ★★★★ |
| 虚拟扰动方向 | 第14章（GRN / 虚拟 KO） | ★★★ 新兴方向 |

### 几点建议

1. **不要试图一口气学完所有分支**——先走通一条主线，再横向拓展
2. **RNA-seq 是必学的**——它是最基础也是应用最广的分析，很多其他分析（免疫浸润、WGCNA、GSEA）都建立在 RNA-seq 数据之上
3. **边学边做**——每学完一章就找一个公共数据集练手，光看不练是学不会的
4. **善用公共数据**——GEO、TCGA 里有海量免费数据，完全可以用来练习和发论文
5. **这份技能树会随着你的学习不断更新**——当你学完某个分支后再回来看，会有完全不同的理解

---

> **准备好了吗？让我们从第一站开始！** → [第1章 生命的语言：从 DNA 到蛋白质](01-dna-to-protein.md)
