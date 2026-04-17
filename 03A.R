# 03A_qc_and_deseq2_pipeline.R
# QC + DESeq2 helpers: tabla counts, VST, PCA, biplot, correlación (1 - r)

suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

project_root <- "C:/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"
setwd(project_root)

dir.create("deseq2_out", showWarnings = FALSE, recursive = TRUE)

dds_path <- "deseq2_out/dds.rds"
stopifnot(file.exists(dds_path))
dds <- readRDS(dds_path)

# --- mapping pretty_id ---
mapping_path <- "deseq2_out/gene_id_mapping.tsv"
if (!file.exists(mapping_path)) {
  gene_ids <- rownames(dds)
  mapping <- tibble(
    gene_id = gene_ids,
    pretty_id = paste0("gene-IR01_", sprintf("%05d", seq_along(gene_ids)))
  )
  write_tsv(mapping, mapping_path)
}
map <- read_tsv(mapping_path, show_col_types = FALSE)

# --- counts table (raw counts) ---
cts <- counts(dds, normalized = FALSE)
cts_df <- as.data.frame(cts) %>%
  rownames_to_column("gene_id") %>%
  left_join(map, by = "gene_id") %>%
  mutate(pretty_id = coalesce(pretty_id, gene_id)) %>%
  select(pretty_id, all_of(colnames(cts)))

write_tsv(cts_df, "deseq2_out/sample_gene_counts.tsv")

# --- VST (blind) ---
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)
saveRDS(vst_mat, "deseq2_out/vst_assay.rds")

# --- PCA using top 500 variable genes ---
rv <- matrixStats::rowVars(vst_mat)
top500 <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
mat500 <- t(vst_mat[top500, , drop = FALSE])

pca500 <- prcomp(mat500, scale. = FALSE, center = TRUE)
coldata <- as.data.frame(colData(dds)) %>% rownames_to_column("sample_id")

pca_samples <- as.data.frame(pca500$x) %>%
  rownames_to_column("sample") %>%
  left_join(coldata, by = "sample")

expl_var <- (pca500$sdev^2) / sum(pca500$sdev^2)
attr(pca_samples, "explained_variance") <- expl_var
saveRDS(pca_samples, "deseq2_out/pca_samples_500var.rds")

# --- Biplot genes: top 50 variable genes among top500 ---
top50 <- top500[seq_len(min(50, length(top500)))]
genes50 <- rownames(vst_mat)[top50]

loadings <- as.data.frame(pca500$rotation[, 1:2, drop = FALSE]) %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id %in% genes50) %>%
  left_join(map, by = "gene_id") %>%
  mutate(pretty_id = coalesce(pretty_id, gene_id)) %>%
  select(gene_id, pretty_id, PC1, PC2)

scores <- pca_samples %>% select(sample, PC1, PC2)
saveRDS(list(scores = scores, loadings = loadings, expl_var = expl_var),
        "deseq2_out/pca_biplot_50genes.rds")

# --- Pearson correlation between samples using VST ---
cor_mat <- cor(vst_mat, method = "pearson")
cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("sample1") %>%
  pivot_longer(-sample1, names_to = "sample2", values_to = "pearson_r")

write_tsv(cor_df, "deseq2_out/sample_cor_pearson.tsv")

# Requested transform: 1 - r (standard and defined for [-1,1])
cor_df2 <- cor_df %>%
  mutate(value = 1 - pearson_r) %>%
  select(sample1, sample2, value)

write_tsv(cor_df2, "deseq2_out/sample_cor_1minusR.tsv")

message("OK: QC objects saved in deseq2_out/")


###########################################################################

