suppressPackageStartupMessages({
  library(data.table)
})

proj_dir <- "/mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"

deg_file <- file.path(proj_dir, "deseq2_out", "DEG_As_vs_Control_sig_padj0.05.tsv")
mapping  <- file.path(proj_dir, "deseq2_out", "gene_id_mapping.tsv")
master   <- file.path(proj_dir, "integrated_annotations", "master_annotations_pretty_id.tsv")

out_dir <- file.path(proj_dir, "deseq2_out", "annotated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# read
deg <- fread(deg_file)
map <- fread(mapping)
ann <- fread(master)

# normalize IDs
deg[, protein_id := sub("^lcl\\|", "", gene_id)]
map[, protein_id := sub("^lcl\\|", "", gene_id)]

# add pretty_id to DEG
deg2 <- merge(deg, map[, .(protein_id, pretty_id)], by="protein_id", all.x=TRUE)

# join annotations by pretty_id
deg_annot <- merge(deg2, ann, by="pretty_id", all.x=TRUE)

# move columns to front
front <- c("pretty_id","gene_id","protein_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
front <- intersect(front, names(deg_annot))
setcolorder(deg_annot, c(front, setdiff(names(deg_annot), front)))

out_file <- file.path(out_dir, basename(deg_file))
fwrite(deg_annot, out_file, sep="\t", quote=FALSE, na="")

cat("Wrote:", out_file, "\n")
cat("Rows:", nrow(deg_annot), "\n")
cat("Rows with pretty_id:", sum(!is.na(deg_annot$pretty_id)), "\n")
cat("Rows with IPS GO:", sum(!is.na(deg_annot$ips_go)), "\n")
