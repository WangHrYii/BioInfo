# 教程配套示例

所有命令从项目根目录执行。默认将生成数据放入 `demo-data/`，分析输出放入 `results/`。不要把生成的示例数据与真实研究数据混放，也不要把人工数值当作实验结论。

| 示例 | 对应章节 | 依赖与数据 | 运行方式 | 默认输出 |
| --- | --- | --- | --- | --- |
| [生成小数据](make_demo_data.py ':ignore') | 6–8、10–11、19 | Python 3 标准库；人工数据，固定种子 42 | `python examples/make_demo_data.py --out demo-data` | 小参考、双端 FASTQ、独立演示 GTF/BED/VCF/计数表 |
| [双端统计](fastq_summary.py ':ignore') | 19 | Python 3 标准库；已生成 FASTQ | 见第 19 章工作流，或使用下方命令 | 每样本 TSV |
| [RNA-seq](rnaseq_airway.R ':ignore') | 9 | R 4.4 / Bioconductor 3.20 教学基线；airway、DESeq2、clusterProfiler、org.Hs.eg.db | `Rscript examples/rnaseq_airway.R` | `results/rnaseq/` |
| [单细胞](scrna_pbmc.py ':ignore') | 12 | Python 3.11 教学基线，包版本见本章；首次联网下载公开 PBMC3k | `python examples/scrna_pbmc.py` | `results/scrna/` |
| [空间统计](spatial_demo.py ':ignore') | 13 | Python 3、NumPy、Matplotlib；人工网格 | `python examples/spatial_demo.py` | `results/spatial/` |
| [图表图例](plot_gallery.py ':ignore') | 20 | Python 3、NumPy、Matplotlib；人工数据和演示统计 | `python examples/plot_gallery.py` | `results/figures/`，三套 PNG/PDF |
| [Snakemake](workflow/Snakefile ':ignore') | 19 | Linux/WSL、Python 3、Snakemake 8.25.5 基线 | `snakemake --snakefile examples/workflow/Snakefile --cores 2` | `results/workflow/` |
| [Nextflow](workflow/main.nf ':ignore') | 19 | Linux/WSL、兼容 Java、Nextflow 24.10 基线、Python 3 | `nextflow run examples/workflow/main.nf` | `results/nextflow/` |

不依赖工作流管理器，也可以直接检查一个样本：

```bash
python examples/fastq_summary.py --sample control --r1 demo-data/control_R1.fastq.gz --r2 demo-data/control_R2.fastq.gz --out results/workflow/control.tsv
```

每个样本应为 20 对、40 条 reads、3000 个碱基。control 的 Q30 比例是 1，treated 是 0.99。FASTQ 是相同人工参考抽出的片段；VCF 与计数表是独立的格式示例，未由这些片段推导。两个 FASTQ 样本不用于差异表达检验。

真实 airway 和 PBMC3k 脚本使用公开数据入口，仍受包版本、网络与平台影响；各章节明确了分析假设、输出和教学范围。空间统计与图表脚本完全使用人工数据，可离线生成图。网站中的 `assets/tutorial/` 图例由相应脚本生成，不能作为任何生物学发现的证据。

## 本次验证记录（2026-09-10）

- 在 Windows / Python 3.8.5 下运行数据生成、双端统计、空间统计和绘图脚本；检查正常结果与双端记录数/名称不匹配时报错。章节内 7 段 Python 练习均执行通过。
- 在独立 Windows / Python 3.10 环境中，以第 12 章列出的 Scanpy/AnnData/NumPy/SciPy/igraph/leidenalg 版本运行 PBMC3k 脚本，保留 2638 个细胞、13656 个基因，生成质控、聚类、marker 和 h5ad 输出。教学安装建议 Python 3.11，本次实际验证为 3.10。
- 使用浏览器检查总览、第 2 章、第 6–13 章、第 19–20 章；18 个 Mermaid 图可解析，配套图片正常加载；检查 117 处本地链接目标存在。
- 15 段 Bash 命令块通过语法解析。当前机器未配置 R 和完整 Linux 生信环境，未实际执行 airway R 脚本、FastQC/fastp、比对/变异命令、Snakemake/Nextflow、容器及集群流程；语法检查不等于这些流程运行成功。
