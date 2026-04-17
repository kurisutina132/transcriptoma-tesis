suppressPackageStartupMessages({
  library(data.table)
})

proj_dir <- "/mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"

deg_dir  <- file.path(proj_dir, "deseq2_out")
map_file <- file.path(proj_dir, "deseq2_out", "gene_id_mapping.tsv")
ann_file <- file.path(proj_dir, "integrated_annotations", "master_annotations_pretty_id.tsv")

out_dir <- file.path(deg_dir, "annotated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# load shared tables once
map <- fread(map_file)
ann <- fread(ann_file)

map[, protein_id := sub("^lcl\\|", "", gene_id)]
map <- unique(map[, .(protein_id, pretty_id)])

# files to annotate: significant only
deg_files <- list.files(deg_dir, pattern = "_sig_padj0\\.05\\.tsv$", full.names = TRUE)

if (length(deg_files) == 0) {
  stop("No files matching *_sig_padj0.05.tsv found in deseq2_out/")
}

for (f in deg_files) {
  message("Annotating: ", basename(f))
  deg <- fread(f)

  if (!("gene_id" %in% names(deg))) {
    warning("Skipping (no gene_id column): ", basename(f))
    next
  }

  deg[, protein_id := sub("^lcl\\|", "", gene_id)]

  # add pretty_id
  deg2 <- merge(deg, map, by = "protein_id", all.x = TRUE)

  # join annotations by pretty_id
  out <- merge(deg2, ann, by = "pretty_id", all.x = TRUE)

  # reorder columns (best-effort)
  front <- c("pretty_id","gene_id","protein_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  front <- intersect(front, names(out))
  setcolorder(out, c(front, setdiff(names(out), front)))

  out_file <- file.path(out_dir, basename(f))
  fwrite(out, out_file, sep = "\t", quote = FALSE, na = "")

  message("  wrote: ", out_file,
          " | rows=", nrow(out),
          " | with_pretty_id=", sum(!is.na(out$pretty_id)),
          " | with_ips_go=", if ("ips_go" %in% names(out)) sum(!is.na(out$ips_go)) else NA)
}

message("Done. Output dir: ", out_dir)
