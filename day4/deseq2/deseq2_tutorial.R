# Install packages (first time only)
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "airway", "apeglm", "clusterProfiler",
                       "org.Hs.eg.db", "ggplot2", "pheatmap"))

library(DESeq2)
library(airway)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(pheatmap)
# ── 1. Load the airway dataset ─────────────────────────────────────
data("airway")
se <- airway                     # SummarizedExperiment object

# What does the data look like?
colData(se)                      # sample metadata
dim(assay(se))                   # genes × samples
summary(se)

# ── 2. Create DESeq2 object ───────────────────────────────────────
dds <- DESeqDataSet(se, design = ~ cell + dex)
# design: ~ cell + dex means:
# account for cell line differences, test effect of dexamethasone (dex)

# Set the reference level (control = "untrt")
dds$dex <- relevel(dds$dex, ref = "untrt")

# ── 3. Pre-filtering ──────────────────────────────────────────────
# Remove genes with very few counts (speeds up analysis, improves FDR)
keep <- rowSums(counts(dds) >= 10) >= 3   # keep if ≥10 counts in ≥3 samples
dds <- dds[keep, ]
nrow(dds)                        # how many genes remain?

# ── 4. Run DESeq2 ─────────────────────────────────────────────────
# This one function does: normalisation + dispersion estimation + statistical testing
dds <- DESeq(dds)

# What happened? Inspect the steps:
sizeFactors(dds)                 # library size factors per sample
head(dispersions(dds))           # gene-level dispersion estimates

# ── 5. Extract Results ────────────────────────────────────────────
res <- results(dds,
               contrast = c("dex", "trt", "untrt"),  # treated vs untreated
               alpha = 0.05)                           # FDR threshold

summary(res)
# Shows: how many up, down, filtered, outliers

head(res[order(res$padj), ])     # top DE genes by adjusted p-value

# Key columns in results:
# baseMean      = average normalised count across all samples
# log2FoldChange = log2(treated / untreated) — positive = upregulated in treated
# lfcSE         = standard error of log2FC
# stat          = Wald test statistic
# pvalue        = raw p-value
# padj          = Benjamini-Hochberg adjusted p-value (USE THIS)

# ── 6. LFC Shrinkage ─────────────────────────────────────────────
library(apeglm)
res_shrunk <- lfcShrink(dds,
                        coef = "dex_trt_vs_untrt",  # the contrast
                        type = "apeglm")             # recommended shrinkage method
# Use res_shrunk for volcano plots and GSEA ranking

# ── 7. Volcano Plot ───────────────────────────────────────────────
res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)
res_df$significant <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_colour_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano plot: dexamethasone treatment",
       x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_minimal()

# ── 8. MA Plot ────────────────────────────────────────────────────
# MA plot: log2FC (y) vs mean expression (x)
# Good genes to scrutinise: high expression + big fold change
plotMA(res_shrunk, ylim = c(-4, 4))

# ── 9. PCA Plot ───────────────────────────────────────────────────
vsd <- vst(dds, blind = FALSE)   # variance-stabilising transformation
plotPCA(vsd, intgroup = c("dex", "cell"))
# You should see treated and untreated samples separate clearly

# ── 10. Heatmap of top DE genes ───────────────────────────────────
# Get top 30 DE genes
top30 <- head(order(res$padj, na.last = NA), 30)
mat <- assay(vsd)[top30, ]
mat <- mat - rowMeans(mat)       # centre by gene mean

annotation_col <- as.data.frame(colData(vsd)[, c("cell", "dex")])
pheatmap(mat,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 8,
         main = "Top 30 DE genes — dex vs untrt")

# ── 11. GO Analysis with clusterProfiler ──────────────────────────
library(clusterProfiler)
library(org.Hs.eg.db)

# Get significant genes — upregulated in treated
sig_up <- rownames(res)[which(!is.na(res$padj) &
                                res$padj < 0.05 &
                                res$log2FoldChange > 1)]

# Convert Ensembl IDs to Entrez IDs (required for clusterProfiler)
sig_entrez <- bitr(sig_up,
                   fromType = "ENSEMBL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

# Run GO over-representation analysis
go_bp <- enrichGO(gene         = sig_entrez$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",           # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)            # show gene symbols

head(go_bp)
dotplot(go_bp, showCategory = 15, title = "GO Biological Process — upregulated in dex")
barplot(go_bp, showCategory = 15)

# ── 12. GSEA ─────────────────────────────────────────────────────
# Build ranked gene list: all genes, ranked by shrunk log2FC
all_genes <- res_shrunk$log2FoldChange
names(all_genes) <- rownames(res_shrunk)
all_genes <- all_genes[!is.na(all_genes)]
all_genes_sorted <- sort(all_genes, decreasing = TRUE)


# Merge the gene map with your fold changes
all_genes_df <- data.frame(
  ENSEMBL = names(all_genes_sorted),
  LFC     = all_genes_sorted
)

gene_map <- bitr(names(all_genes_sorted),
                 fromType = "ENSEMBL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

merged <- merge(gene_map, all_genes_df, by = "ENSEMBL")

# Remove duplicates — keep the entry with highest absolute LFC per Entrez ID
merged <- merged[order(abs(merged$LFC), decreasing = TRUE), ]
merged <- merged[!duplicated(merged$ENTREZID), ]

# Build the ranked vector — no duplicates now
ranked <- merged$LFC
names(ranked) <- merged$ENTREZID
ranked <- sort(ranked, decreasing = TRUE)

# Verify no duplicates remain
sum(duplicated(names(ranked)))   # must be 0

# Now re-run gseGO
gsea_bp <- gseGO(geneList      = ranked,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 minGSSize     = 15,
                 maxGSSize     = 500,
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 verbose       = FALSE)

gseaplot2(gsea_bp, geneSetID = 1, title = gsea_bp$Description[1])

# ── 13. Save results ─────────────────────────────────────────────
write.csv(as.data.frame(res_shrunk),
          "DESeq2_dex_vs_untrt_results.csv",
          row.names = TRUE)
write.csv(as.data.frame(go_bp),
          "GO_BP_upregulated.csv",
          row.names = FALSE)


