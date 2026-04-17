#!/usr/bin/env bash
set -euo pipefail

out="reports/transport_full/transport_report.Rmd"
mkdir -p "$(dirname "$out")"
: > "$out"

w() { printf "%s\n" "$1" >> "$out"; }

# YAML
w "---"
w 'title: "Transport DEGs report (Metal / MDR / Nutrient / Other)"'
w 'date: "`r format(Sys.Date(), '\''%Y-%m-%d'\'')`"'
w "output:"
w "  html_document:"
w "    theme: readable"
w "    toc: true"
w "    toc_float: true"
w "    df_print: paged"
w "---"
w ""

# setup
w '```{r setup, include=FALSE}'
w 'knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)'
w 'userlib <- path.expand("~/R/library")'
w 'if (dir.exists(userlib)) .libPaths(c(userlib, .libPaths()))'
w 'library(data.table)'
w 'library(knitr)'
w 'library(DT)'
w 'proj_dir <- "/mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"'
w 'out_dir  <- file.path(proj_dir, "reports", "transport_full")'
w 'comparisons <- c("As_vs_Control", "Cr_vs_Control", "As_vs_Cr")'
w 'subcats <- c("Metal", "MDR", "Nutrient", "Other")'
w 'thr_padj <- 0.05'
w 'thr_lfc  <- 1'
w 'img <- function(name) file.path(out_dir, name)'
w 'tsv <- function(name) file.path(out_dir, name)'
w 'counts_file <- tsv("transport_counts_by_subcat.tsv")'
w 'stopifnot(file.exists(counts_file))'
w 'counts <- fread(counts_file)'
w 'counts_wide <- dcast(counts, comparison ~ subcat, value.var="n", fill=0)'
w 'fmt_p <- function(x) { if (is.na(x)) return(NA_character_); if (x < 0.001) return(format(x, scientific=TRUE, digits=2)); format(round(x,4), nsmall=4) }'
w 'pick_col <- function(dt, candidates) { hit <- candidates[candidates %in% names(dt)][1]; if (is.na(hit)) return(NULL); hit }'
w 'make_gene_label <- function(dt) {'
w '  idcol <- pick_col(dt, c("id","pretty_id","gene_id"))'
w '  desccol <- pick_col(dt, c("description","product","annotation","gene_name","Gene","name"))'
w '  if (is.null(idcol)) idcol <- names(dt)[1]'
w '  if (!is.null(desccol)) paste0(dt[[idcol]], " — ", dt[[desccol]]) else as.character(dt[[idcol]])'
w '}'
w '```'
w ""

# Executive summary
w "## Executive summary"
w ""
w '```{r}'
w 'exec <- copy(counts)'
w 'exec[, n := as.integer(n)]'
w 'exec_tot <- exec[, .(total_transport = sum(n, na.rm=TRUE)), by=comparison]'
w 'exec2 <- merge(exec, exec_tot, by="comparison", all.x=TRUE)'
w 'exec2[, pct := ifelse(total_transport > 0, round(100*n/total_transport, 1), 0)]'
w 'setorder(exec2, comparison, -n)'
w 'topcat <- exec2[, .SD[1], by=comparison]'
w 'summary_tbl <- merge(exec_tot, topcat[, .(comparison, top_subcat=subcat, top_n=n, top_pct=pct)], by="comparison")'
w 'DT::datatable(summary_tbl, rownames=FALSE, options=list(pageLength=10))'
w '```'
w ""

# Counts
w "## 1) Counts by subcategory"
w ""
w '```{r}'
w 'DT::datatable(counts_wide, rownames=FALSE, options=list(pageLength=10))'
w '```'
w ""
w '```{r}'
w 'if (file.exists(img("transport_counts_by_subcat.png"))) knitr::include_graphics("transport_counts_by_subcat.png")'
w '```'
w ""

# Plots (use markdown image embedding)
w "## 2) Plots per comparison"
w ""
w '```{r, results="asis"}'
w 'for (comp in comparisons) {'
w '  cat("\n\n### ", comp, "\n\n", sep="")'
w '  sub <- counts[comparison == comp]'
w '  if (nrow(sub) > 0) {'
w '    tot <- sum(sub$n, na.rm=TRUE)'
w '    cat("**Summary:** ", tot, " Transport DEGs. Volcano thresholds: `padj < ", thr_padj, "` and `|log2FC| ≥ ", thr_lfc, "`.\n\n", sep="")'
w '  }'
w '  plots <- c(sprintf("%s_volcano_transport.png", comp), sprintf("%s_MA_transport.png", comp), sprintf("%s_log2fc_hist.png", comp), sprintf("%s_log2fc_box.png", comp))'
w '  for (p in plots) if (file.exists(img(p))) cat("![](", p, ")\n\n", sep="")'
w '}'
w '```'
w ""

# Auto-generated results text
w "## 3) Auto-generated Results text (descriptive)"
w ""
w '```{r, results="asis"}'
w 'topN <- 4'
w 'for (comp in comparisons) {'
w '  cat("\n\n### ", comp, " — Results text\n\n", sep="")'
w '  for (s in subcats) {'
w '    up_f <- sprintf("%s_%s_top_up.tsv", comp, s)'
w '    dn_f <- sprintf("%s_%s_top_down.tsv", comp, s)'
w '    if (!file.exists(tsv(up_f)) && !file.exists(tsv(dn_f))) next'
w '    cat("#### ", s, "\n\n", sep="")'
w '    if (file.exists(tsv(up_f))) cat("- Top up TSV: [", up_f, "](", up_f, ")\n", sep="")'
w '    if (file.exists(tsv(dn_f))) cat("- Top down TSV: [", dn_f, "](", dn_f, ")\n\n", sep="")'
w '    if (file.exists(tsv(up_f))) {'
w '      up <- fread(tsv(up_f)); up <- up[order(-as.numeric(log2FoldChange))]; up <- head(up, topN)'
w '      labs <- make_gene_label(up); lfc <- as.numeric(up$log2FoldChange); pad <- as.numeric(up$padj)'
w '      items <- paste0(labs, " (log2FC=", round(lfc,2), "; padj=", vapply(pad, fmt_p, character(1)), ")")'
w '      cat("**Most up-regulated (top ", topN, "):** ", paste(items, collapse="; "), ".\n\n", sep="")'
w '    }'
w '    if (file.exists(tsv(dn_f))) {'
w '      dn <- fread(tsv(dn_f)); dn <- dn[order(as.numeric(log2FoldChange))]; dn <- head(dn, topN)'
w '      labs <- make_gene_label(dn); lfc <- as.numeric(dn$log2FoldChange); pad <- as.numeric(dn$padj)'
w '      items <- paste0(labs, " (log2FC=", round(lfc,2), "; padj=", vapply(pad, fmt_p, character(1)), ")")'
w '      cat("**Most down-regulated (top ", topN, "):** ", paste(items, collapse="; "), ".\n\n", sep="")'
w '    }'
w '  }'
w '}'
w '```'
w ""

# Interactive tables
w "## 4) Top tables (interactive)"
w "Use the table search box (e.g., TRI12, ABCG, P-type, ZIP, Ctr) to filter."
w ""
w '```{r, results="asis"}'
w 'for (comp in comparisons) {'
w '  cat("\n\n### ", comp, "\n\n", sep="")'
w '  for (s in subcats) {'
w '    up_f <- sprintf("%s_%s_top_up.tsv", comp, s)'
w '    dn_f <- sprintf("%s_%s_top_down.tsv", comp, s)'
w '    if (!file.exists(tsv(up_f)) && !file.exists(tsv(dn_f))) next'
w '    cat("\n\n#### ", s, "\n\n", sep="")'
w '    if (file.exists(tsv(up_f))) cat("- Top up TSV: [", up_f, "](", up_f, ")\n", sep="")'
w '    if (file.exists(tsv(dn_f))) cat("- Top down TSV: [", dn_f, "](", dn_f, ")\n\n", sep="")'
w '    if (file.exists(tsv(up_f))) { up <- fread(tsv(up_f)); print(DT::datatable(up, rownames=FALSE, options=list(pageLength=10, scrollX=TRUE))); cat("\n\n") }'
w '    if (file.exists(tsv(dn_f))) { dn <- fread(tsv(dn_f)); print(DT::datatable(dn, rownames=FALSE, options=list(pageLength=10, scrollX=TRUE))); cat("\n\n") }'
w '  }'
w '}'
w '```'

echo "Wrote $out"
