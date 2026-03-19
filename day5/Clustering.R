
# 1. Load libraries
library(DESeq2)
library(airway)

# 2. Load data + create DESeq object
data("airway")
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, ref = "untrt")

# 3. Filter + run DESeq2
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

# 4. VST — needed before any PCA, heatmap, or clustering
vsd <- vst(dds, blind = FALSE)

# 5. Now this works
mat <- assay(vsd)

##Ch-6 Clustering

# On VST-transformed DESeq2 data
mat <- assay(vsd)                        # genes × samples matrix
mat_scaled <- t(scale(t(mat)))           # scale genes to mean 0, sd 1

# Compute distance between samples
dist_samples <- dist(t(mat_scaled), method = "euclidean")
hclust_samples <- hclust(dist_samples, method = "ward.D2")
plot(hclust_samples, main = "Sample dendrogram")

# Create the annotation from the vsd object itself
annotation_col <- as.data.frame(colData(vsd)[, c("cell", "dex")])

# Heatmap with both gene and sample clustering
pheatmap(mat_scaled,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Clustered heatmap — VST data")

###k-means clustering

set.seed(42)                             # for reproducibility
km <- kmeans(mat_scaled, centers = 4, nstart = 25)
# nstart = 25: try 25 random starts, keep the best

km$cluster                               # which cluster each gene belongs to
km$centers                               # centroid coordinates
km$withinss                              # within-cluster sum of squares per cluster
km$tot.withinss                          # total — minimise this

# Add cluster labels to your gene table
gene_clusters <- data.frame(
  gene = rownames(mat_scaled),
  cluster = km$cluster
)
table(gene_clusters$cluster)             # how many genes per cluster?





