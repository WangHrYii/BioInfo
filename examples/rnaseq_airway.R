# Run from the repository root: Rscript examples/rnaseq_airway.R
# Required Bioconductor packages: airway, DESeq2, clusterProfiler, org.Hs.eg.db.
suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
out <- file.path("results", "rnaseq")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
data("airway", package = "airway")
se <- airway
se$dex <- relevel(factor(se$dex), ref = "untrt")
se$cell <- factor(se$cell)
stopifnot(identical(colnames(assay(se)), rownames(as.data.frame(colData(se)))))
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
tab <- as.data.frame(res)
tab$gene_id <- rownames(tab)
tab <- tab[order(tab$padj, na.last = TRUE), ]
write.csv(tab, file.path(out, "differential_expression.csv"), row.names = FALSE)
write.csv(counts(dds, normalized = TRUE), file.path(out, "normalized_counts.csv"))
write.csv(as.data.frame(colData(dds)), file.path(out, "sample_metadata.csv"))
vsd <- vst(dds, blind = FALSE)
pdf(file.path(out, "pca.pdf"), width = 7, height = 5)
print(plotPCA(vsd, intgroup = c("dex", "cell")))
dev.off()
pdf(file.path(out, "ma.pdf"), width = 7, height = 5)
plotMA(res, ylim = c(-5, 5), main = "Treated versus untreated")
dev.off()
valid <- is.finite(tab$log2FoldChange) & !is.na(tab$padj)
plot_tab <- tab[valid, ]
selected <- plot_tab$padj < 0.05 & abs(plot_tab$log2FoldChange) >= 1
pdf(file.path(out, "volcano.pdf"), width = 7, height = 5)
plot(plot_tab$log2FoldChange, -log10(pmax(plot_tab$padj, 1e-300)),
     pch = 16, cex = 0.35, col = ifelse(selected, "#D55E00", "#999999"),
     xlab = "log2 fold change (treated / untreated)", ylab = "-log10 adjusted p-value")
abline(v = c(-1, 1), h = -log10(0.05), lty = 2, col = "#666666")
dev.off()

# Over-representation analysis: tested, mappable genes form the background.
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(dds),
                                columns = "ENTREZID", keytype = "ENSEMBL")
mapping <- unique(mapping[!is.na(mapping$ENTREZID), c("ENSEMBL", "ENTREZID")])
tested <- rownames(res)[!is.na(res$padj)]
up <- rownames(res)[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >= 1]
background <- unique(mapping$ENTREZID[mapping$ENSEMBL %in% tested])
up_ids <- unique(mapping$ENTREZID[mapping$ENSEMBL %in% up])
if (length(up_ids) > 0 && length(background) > length(up_ids)) {
  go <- enrichGO(gene = up_ids, universe = background, OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  write.csv(as.data.frame(go), file.path(out, "go_bp_up.csv"), row.names = FALSE)
} else {
  message("No eligible upregulated genes for GO analysis; no GO result written.")
}
saveRDS(dds, file.path(out, "dds.rds"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("Results saved to ", normalizePath(out))
