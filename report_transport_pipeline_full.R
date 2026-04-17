suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

proj_dir <- "/mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"

# full DESeq2 (background)
full_deseq_files <- list(
  As_vs_Control = file.path(proj_dir, "deseq2_out", "DEG_As_vs_Control.tsv"),
  Cr_vs_Control = file.path(proj_dir, "deseq2_out", "DEG_Cr_vs_Control.tsv"),
  As_vs_Cr      = file.path(proj_dir, "deseq2_out", "DEG_As_vs_Cr.tsv")
)

transport_subcats_dir <- file.path(proj_dir, "gene_lists", "transport_subcats")
out_dir <- file.path(proj_dir, "reports", "transport_full")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- c("As_vs_Control", "Cr_vs_Control", "As_vs_Cr")
subcats <- c("Metal", "MDR", "Nutrient", "Other")

# fixed palette (as requested)
palette_fixed <- c(
  Background = "grey80",
  Metal      = "#d73027",  # red
  MDR        = "#4575b4",  # blue
  Nutrient   = "#1a9850",  # green
  Other      = "#762a83"   # purple/violet
)

# volcano thresholds
thr_padj <- 0.05
thr_lfc  <- 1

safe_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) return(NULL)
  hit
}

read_subcat <- function(comp, subcat) {
  f <- file.path(transport_subcats_dir, sprintf("%s_Transport_%s.tsv", comp, subcat))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)

  idcol <- safe_col(dt, c("pretty_id", "gene_id"))
  if (is.null(idcol)) stop(paste("No pretty_id/gene_id in", f))
  setnames(dt, idcol, "id")

  dt[, subcat := subcat]
  dt[, comparison := comp]
  dt[, log2FoldChange := as.numeric(log2FoldChange)]
  dt[, padj := as.numeric(padj)]
  dt
}

make_overlay <- function(comp) {
  dt <- rbindlist(lapply(subcats, function(s) read_subcat(comp, s)), fill=TRUE)
  if (nrow(dt) == 0) return(NULL)
  unique(dt[, .(id, subcat)])
}

save_plot_dual <- function(filename_stem, plot, w, h) {
  ggsave(file.path(out_dir, paste0(filename_stem, ".png")), plot, width=w, height=h, dpi=300)
  ggsave(file.path(out_dir, paste0(filename_stem, ".pdf")), plot, width=w, height=h, device = cairo_pdf)
}

# ---------- 1) Conteos por subcat ----------
counts <- rbindlist(lapply(comparisons, function(comp) {
  rbindlist(lapply(subcats, function(s) {
    dt <- read_subcat(comp, s)
    if (is.null(dt)) return(NULL)
    data.table(comparison=comp, subcat=s, n=nrow(dt))
  }))
}), fill=TRUE)

fwrite(counts, file.path(out_dir, "transport_counts_by_subcat.tsv"), sep="\t", quote=FALSE)

p_counts <- ggplot(counts, aes(x=subcat, y=n, fill=subcat)) +
  geom_col(color="grey20") +
  facet_wrap(~comparison, scales="free_y") +
  theme_bw(base_size=12) +
  theme(legend.position="none") +
  scale_fill_manual(values=palette_fixed[names(palette_fixed) != "Background"]) +
  labs(title="Transport DEGs by subcategory", x=NULL, y="# genes")

save_plot_dual("transport_counts_by_subcat", p_counts, w=10, h=4)

# ---------- 2) Top up/down (por subcat) ----------
write_top_tables <- function(dt, comp, subcat, top_n=10) {
  if (is.null(dt) || nrow(dt)==0) return(invisible(NULL))
  up <- dt[order(-log2FoldChange)][1:min(top_n, .N)]
  dn <- dt[order(log2FoldChange)][1:min(top_n, .N)]
  fwrite(up, file.path(out_dir, sprintf("%s_%s_top_up.tsv", comp, subcat)), sep="\t", quote=FALSE)
  fwrite(dn, file.path(out_dir, sprintf("%s_%s_top_down.tsv", comp, subcat)), sep="\t", quote=FALSE)
}
for (comp in comparisons) for (s in subcats) write_top_tables(read_subcat(comp, s), comp, s, 10)

# ---------- 3) Distribuciones log2FC ----------
all_transport <- rbindlist(lapply(comparisons, function(comp) {
  rbindlist(lapply(subcats, function(s) read_subcat(comp, s)), fill=TRUE)
}), fill=TRUE)

for (comp in comparisons) {
  dtc <- all_transport[comparison==comp]
  if (nrow(dtc)==0) next

  p_hist <- ggplot(dtc, aes(x=log2FoldChange, fill=subcat)) +
    geom_histogram(bins=40, alpha=0.65, position="identity") +
    theme_bw(base_size=12) +
    scale_fill_manual(values=palette_fixed[names(palette_fixed) != "Background"]) +
    labs(title=paste0(comp, " — log2FC distribution (Transport)"),
         x="log2FoldChange", y="Count")
  save_plot_dual(sprintf("%s_log2fc_hist", comp), p_hist, w=8, h=5)

  p_box <- ggplot(dtc, aes(x=subcat, y=log2FoldChange, fill=subcat)) +
    geom_boxplot(outlier.size=0.6) +
    theme_bw(base_size=12) +
    theme(legend.position="none") +
    scale_fill_manual(values=palette_fixed[names(palette_fixed) != "Background"]) +
    labs(title=paste0(comp, " — log2FC by subcategory (Transport)"), x=NULL, y="log2FoldChange")
  save_plot_dual(sprintf("%s_log2fc_box", comp), p_box, w=7, h=5)
}

# ---------- 4) MA + Volcano con background full ----------
for (comp in comparisons) {
  full_f <- full_deseq_files[[comp]]
  if (!file.exists(full_f)) stop(paste("Missing full DESeq2 file:", full_f))
  full <- fread(full_f)

  idcol   <- safe_col(full, c("pretty_id","gene_id","gene","id"))
  lfc_col  <- safe_col(full, c("log2FoldChange","log2FC","logFC"))
  padj_col <- safe_col(full, c("padj","FDR","adj.P.Val"))
  base_col <- safe_col(full, c("baseMean","base_mean","mean"))

  if (is.null(idcol) || is.null(lfc_col) || is.null(padj_col) || is.null(base_col)) {
    stop(paste0("In ", full_f, " missing one of: id/log2FoldChange/padj/baseMean. Columns are: ",
                paste(names(full), collapse=", ")))
  }

  setnames(full, c(idcol, lfc_col, padj_col, base_col), c("id","log2FoldChange","padj","baseMean"))
  full[, log2FoldChange := as.numeric(log2FoldChange)]
  full[, padj := as.numeric(padj)]
  full[, baseMean := as.numeric(baseMean)]

  overlay <- make_overlay(comp)
  if (!is.null(overlay)) full <- merge(full, overlay, by="id", all.x=TRUE)
  full[, highlight := ifelse(is.na(subcat), "Background", as.character(subcat))]
  full[, highlight := factor(highlight, levels=c("Background", subcats))]

  full[, neglog10padj := -log10(pmax(padj, 1e-300))]
  full[, log10baseMean := log10(pmax(baseMean, 1e-6))]

  # Volcano thresholds
  y_thr <- -log10(thr_padj)

  p_volc <- ggplot(full, aes(x=log2FoldChange, y=neglog10padj)) +
    geom_point(data=full[highlight=="Background"], color=palette_fixed[["Background"]], size=0.55, alpha=0.45) +
    geom_point(data=full[highlight!="Background"], aes(color=highlight), size=1.1, alpha=0.9) +
    geom_hline(yintercept = y_thr, linetype="dashed", color="grey30", linewidth=0.4) +
    geom_vline(xintercept = c(-thr_lfc, thr_lfc), linetype="dashed", color="grey30", linewidth=0.4) +
    theme_bw(base_size=12) +
    scale_color_manual(values=palette_fixed) +
    labs(title=paste0(comp, " — Volcano (Transport highlighted)"),
         x="log2FoldChange", y="-log10(padj)", color="Category")
  save_plot_dual(sprintf("%s_volcano_transport", comp), p_volc, w=7.5, h=6)

  # MA plot (with same coloring)
  p_ma <- ggplot(full, aes(x=log10baseMean, y=log2FoldChange)) +
    geom_point(data=full[highlight=="Background"], color=palette_fixed[["Background"]], size=0.55, alpha=0.45) +
    geom_point(data=full[highlight!="Background"], aes(color=highlight), size=1.1, alpha=0.9) +
    geom_hline(yintercept = 0, color="grey30", linewidth=0.4) +
    theme_bw(base_size=12) +
    scale_color_manual(values=palette_fixed) +
    labs(title=paste0(comp, " — MA plot (Transport highlighted)"),
         x="log10(baseMean)", y="log2FoldChange", color="Category")
  save_plot_dual(sprintf("%s_MA_transport", comp), p_ma, w=7.5, h=6)
}

cat("Done. Wrote tables + plots to:", out_dir, "\n")
