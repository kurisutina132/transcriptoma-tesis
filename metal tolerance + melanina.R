# 03A_qc_and_deseq2_pipeline.R
# Genera artefactos: VST, PCA, biplot, correlación, y tablas DESeq2 por contraste.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

project_root <- "C:/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"
setwd(project_root)

# ----------------------------
# Helpers (reusables)
# ----------------------------
clean_sample_label <- function(x) {
  x <- as.character(x)
  sub("_.*$", "", x)
}

preferred_order <- c("C1_me","C2_me","C3_me","T1_as","T2_as","T3_as","T4_cr","T5_cr","T6_cr")

condition_from_label <- function(lbl) {
  s2 <- tolower(lbl)
  dplyr::case_when(
    str_detect(s2, "^c") ~ "Control",
    str_detect(s2, "as") ~ "As",
    str_detect(s2, "cr") ~ "Cr",
    TRUE ~ "Unknown"
  )
}

# ----------------------------
# Inputs base
# ----------------------------
dds <- readRDS("deseq2_out/dds.rds")

# Map opcional (si no existe, igual corre DESeq2)
map_path <- "deseq2_out/gene_id_mapping.tsv"
map <- if (file.exists(map_path)) read_tsv(map_path, show_col_types = FALSE) else tibble()

# coldata canon
coldata <- as.data.frame(colData(dds)) %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_label = clean_sample_label(sample_id),
         condition_plot = condition_from_label(sample_label))

# ----------------------------
# 1) Guardar VST assay si no existe (recomendado para heatmaps)
# ----------------------------
if (!file.exists("deseq2_out/vst_assay.rds")) {
  vsd <- vst(dds, blind = TRUE)
  saveRDS(assay(vsd), "deseq2_out/vst_assay.rds")
}

# ----------------------------
# 2) Contraste T vs C: crear DEG_T_vs_C.tsv si no existe
# ----------------------------
if (!file.exists("deseq2_out/DEG_T_vs_C.tsv")) {
  sample_ids <- colnames(dds)
  sample_label <- clean_sample_label(sample_ids)
  
  TC <- ifelse(grepl("^C", sample_label, ignore.case = TRUE), "C", "T")
  TC <- factor(TC, levels = c("C", "T"))
  colData(dds)$TC <- TC
  
  dds_tc <- dds
  design(dds_tc) <- ~ TC
  dds_tc <- DESeq(dds_tc)
  
  res_tc <- results(dds_tc, contrast = c("TC", "T", "C")) |>
    as.data.frame() |>
    rownames_to_column("gene_id")
  
  write_tsv(res_tc, "deseq2_out/DEG_T_vs_C.tsv")
  saveRDS(dds_tc, "deseq2_out/dds_tc.rds")
}

# ----------------------------
# 3) Verificación rápida de archivos DE (no crea, solo reporta)
# ----------------------------
deg_files <- list.files("deseq2_out", pattern = "^DEG_.*\\.tsv$", full.names = FALSE)
message("DEG files encontrados:\n - ", paste(deg_files, collapse = "\n - "))

# Nota: As_vs_Control, Cr_vs_Control, As_vs_Cr deben existir ya desde tu pipeline original.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(DT)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

project_root <- "C:/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"
setwd(project_root)

dds <- readRDS("deseq2_out/dds.rds")
map <- read_tsv("deseq2_out/gene_id_mapping.tsv", show_col_types = FALSE)

clean_sample_label <- function(x) sub("_.*$", "", as.character(x))
preferred_order <- c("C1_me","C2_me","C3_me","T1_as","T2_as","T3_as","T4_cr","T5_cr","T6_cr")

condition_from_label <- function(lbl) {
  s2 <- tolower(lbl)
  dplyr::case_when(
    str_detect(s2, "^c") ~ "Control",
    str_detect(s2, "as") ~ "As",
    str_detect(s2, "cr") ~ "Cr",
    TRUE ~ "Unknown"
  )
}

coldata <- as.data.frame(colData(dds)) %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(
    sample_label = clean_sample_label(sample_id),
    condition_plot = condition_from_label(sample_label)
  )

# ---------- Tablas ----------
read_deseq2_tbl <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    left_join(map, by = "gene_id") %>%
    mutate(pretty_id = coalesce(pretty_id, gene_id)) %>%
    transmute(
      `Identificación genética` = pretty_id,
      baseMean = baseMean,
      log2FoldChange = log2FoldChange,
      lfcSE = lfcSE,
      stat = stat,
      pvalue = pvalue,
      padj = padj
    )
}

deseq2_dt <- function(path) {
  DT::datatable(
    read_deseq2_tbl(path),
    rownames = FALSE,
    extensions = c("Buttons", "Scroller"),
    options = list(
      dom="Bfrtip",
      buttons=c("csv","excel"),
      pageLength=10,
      scrollX=TRUE,
      scrollY=350,
      scroller=TRUE
    )
  )
}

summarize_deg <- function(path, padj_cut = 0.05, lfc_cut = 0.3) {
  x <- read_tsv(path, show_col_types = FALSE)
  sig <- x %>% filter(!is.na(padj), padj < padj_cut)
  up <- sig %>% filter(log2FoldChange >= lfc_cut)
  down <- sig %>% filter(log2FoldChange <= -lfc_cut)
  
  tibble(
    genes_en_tabla = nrow(x),
    DEGs_padj_0.05 = nrow(sig),
    sobreexpresados = nrow(up),
    subexpresados = nrow(down)
  )
}

# ---------- MA ----------
ma_plot <- function(path, title) {
  df <- read_tsv(path, show_col_types = FALSE)
  ggplot(df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(alpha = 0.4, size = 1) +
    scale_x_log10() +
    geom_hline(yintercept = 0, color = "grey30") +
    labs(title = title, x = "baseMean (log10)", y = "log2 fold change") +
    theme_bw()
}

# ---------- Volcano 4 clases ----------
volcano_plot <- function(path, title, p_cut = 0.05, lfc_cut = 0.3, use_padj = TRUE) {
  df <- read_tsv(path, show_col_types = FALSE) %>%
    left_join(map, by = "gene_id") %>%
    mutate(pretty_id = coalesce(pretty_id, gene_id))
  
  pcol <- if (use_padj && "padj" %in% names(df)) "padj" else "pvalue"
  
  df <- df %>%
    mutate(
      p_use = .data[[pcol]],
      neglog10p = -log10(pmax(p_use, 1e-300)),
      sig_p = !is.na(p_use) & p_use < p_cut,
      sig_lfc = abs(log2FoldChange) >= lfc_cut,
      cat = case_when(
        sig_p & sig_lfc ~ "p-value and log2FC",
        sig_lfc & !sig_p ~ "Log2 FC",
        sig_p & !sig_lfc ~ "P-value",
        TRUE ~ "NS"
      )
    ) %>%
    filter(is.finite(log2FoldChange), is.finite(neglog10p))
  
  cols <- c(
    "NS" = "grey70",
    "Log2 FC" = "#2ca02c",
    "P-value" = "#1f77b4",
    "p-value and log2FC" = "#d62728"
  )
  
  ggplot(df, aes(log2FoldChange, neglog10p, color = cat)) +
    geom_point(alpha = 0.75, size = 1.3) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey30") +
    scale_color_manual(values = cols, breaks = names(cols)) +
    labs(title = title, x = "Log2 fold change", y = expression(-log[10](p)), color = "") +
    theme_bw()
}

# ---------- Heatmap ComplexHeatmap (top500 + top20, clusters separados) ----------
heatmap_degs_top500_complex_blocks <- function(
    path,
    title,
    n_show = 500,
    n_mark_total = 20,
    n_clusters = 2,
    padj_cut = 0.05,
    anno_fontsize = 8,
    heatmap_height_cm = 22,
    gap_mm = 16
) {
  vst <- readRDS("deseq2_out/vst_assay.rds")
  
  res <- read_tsv(path, show_col_types = FALSE) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    left_join(map, by = "gene_id") %>%
    mutate(pretty_id = coalesce(pretty_id, gene_id))
  
  show_tbl <- res %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = n_show)
  
  mat <- vst[show_tbl$gene_id, , drop = FALSE]
  rownames(mat) <- show_tbl$pretty_id
  
  colnames(mat) <- clean_sample_label(colnames(mat))
  common <- intersect(preferred_order, colnames(mat))
  if (length(common) >= 2) mat <- mat[, common, drop = FALSE]
  
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  
  # column annotation C/T
  col_group <- ifelse(grepl("^C", colnames(mat_z), ignore.case = TRUE), "C", "T")
  col_group <- factor(col_group, levels = c("C","T"))
  ha <- HeatmapAnnotation(
    Condition = col_group,
    col = list(Condition = c(C = "#7FC97F", T = "#BEAED4")),
    annotation_name_side = "right"
  )
  
  col_fun <- colorRamp2(c(-4, 0, 4), c("#3B4CC0", "white", "#B40426"))
  
  ht <- Heatmap(
    mat_z,
    name = "Z-score",
    col = col_fun,
    top_annotation = ha,
    cluster_rows = TRUE,
    row_split = n_clusters,             # genera slices 1..k
    cluster_columns = TRUE,
    show_row_names = FALSE,             # etiquetas van por anno_mark
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 11),
    column_names_rot = 90,
    gap = unit(gap_mm, "mm"),
    height = unit(heatmap_height_cm, "cm"),
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 13, fontface = "bold")
  )
  
  ht_drawn <- draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(4, 20, 4, 6), "mm")
  )
  
  # cluster membership from slice orders
  ord_list <- row_order(ht_drawn)
  pretty_ids <- rownames(mat_z)
  cl_by_gene <- rep(NA_integer_, length(pretty_ids))
  names(cl_by_gene) <- pretty_ids
  for (i in seq_along(ord_list)) cl_by_gene[ pretty_ids[ ord_list[[i]] ] ] <- i
  
  # pick top genes per cluster to reduce overlap
  pool <- res %>%
    filter(padj < padj_cut) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    distinct(pretty_id, .keep_all = TRUE) %>%
    filter(pretty_id %in% pretty_ids) %>%
    mutate(cluster = cl_by_gene[pretty_id])
  
  per_cl <- max(1, floor(n_mark_total / n_clusters))
  marked <- pool %>% group_by(cluster) %>% slice_head(n = per_cl) %>% ungroup()
  if (nrow(marked) < n_mark_total) {
    extra <- pool %>% anti_join(marked, by = "pretty_id") %>%
      slice_head(n = n_mark_total - nrow(marked))
    marked <- bind_rows(marked, extra)
  }
  
  mark_genes <- marked$pretty_id
  
  ra <- rowAnnotation(
    `Top 20` = anno_mark(
      at = mark_genes,
      labels = mark_genes,
      labels_gp = gpar(fontsize = anno_fontsize),
      link_width = unit(14, "mm"),
      link_gp = gpar(col = "grey35", lwd = 0.7),
      padding = unit(1, "mm")
    ),
    width = unit(10, "cm")
  )
  
  draw(
    ht + ra,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(4, 40, 4, 6), "mm")
  )
  
  grid.text(
    title,
    x = unit(0.5, "npc"),
    y = unit(0.98, "npc"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
}


#T vs C (T/C)

# 4.X Tratamiento vs Control (T vs C)

summarize_deg("deseq2_out/DEG_T_vs_C.tsv")
deseq2_dt("deseq2_out/DEG_T_vs_C.tsv")

ma_plot("deseq2_out/DEG_T_vs_C.tsv", "MA plot – T vs C")

volcano_plot(
  "deseq2_out/DEG_T_vs_C.tsv",
  "Volcano – T vs C (padj<0.05; |log2FC|>0.3)"
)

# Recomiendo fig.width grande en el chunk del Rmd
heatmap_degs_top500_complex_blocks(
  "deseq2_out/DEG_T_vs_C.tsv",
  "Heatmap (VST z-score) – Top 500 + Top 20 marcados – T vs C",
  n_show = 500, n_mark_total = 20, n_clusters = 2,
  anno_fontsize = 8, heatmap_height_cm = 22, gap_mm = 16
)


#As vs Cr (As/Cr)

# 4.X As vs Cr

summarize_deg("deseq2_out/DEG_As_vs_Cr.tsv")
deseq2_dt("deseq2_out/DEG_As_vs_Cr.tsv")

ma_plot("deseq2_out/DEG_As_vs_Cr.tsv", "MA plot – As vs Cr")

volcano_plot(
  "deseq2_out/DEG_As_vs_Cr.tsv",
  "Volcano – As vs Cr (padj<0.05; |log2FC|>0.3)"
)

heatmap_degs_top500_complex_blocks(
  "deseq2_out/DEG_As_vs_Cr.tsv",
  "Heatmap (VST z-score) – Top 500 + Top 20 marcados – As vs Cr",
  n_show = 500, n_mark_total = 20, n_clusters = 2,
  anno_fontsize = 8, heatmap_height_cm = 22, gap_mm = 16
)

#Cr vs Control (Cr/control)

# 4.X Cr vs Control

summarize_deg("deseq2_out/DEG_Cr_vs_Control.tsv")
deseq2_dt("deseq2_out/DEG_Cr_vs_Control.tsv")

ma_plot("deseq2_out/DEG_Cr_vs_Control.tsv", "MA plot – Cr vs Control")

volcano_plot(
  "deseq2_out/DEG_Cr_vs_Control.tsv",
  "Volcano – Cr vs Control (padj<0.05; |log2FC|>0.3)"
)

heatmap_degs_top500_complex_blocks(
  "deseq2_out/DEG_Cr_vs_Control.tsv",
  "Heatmap (VST z-score) – Top 500 + Top 20 marcados – Cr vs Control",
  n_show = 500, n_mark_total = 20, n_clusters = 2,
  anno_fontsize = 8, heatmap_height_cm = 22, gap_mm = 16
)


