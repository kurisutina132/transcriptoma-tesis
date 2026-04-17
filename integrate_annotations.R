# integrate_annotations.R
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

proj_dir <- "/mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"

interpro_tsv <- file.path(proj_dir, "interproscan_out", "interproscan.tsv")
mapping_tsv  <- file.path(proj_dir, "deseq2_out", "gene_id_mapping.tsv")
emapper_ann  <- "/home/kurisutina/exophiala/trancriptoma/emapper_out/exophiala_emapper_20260309_1917.emapper.annotations"

out_dir <- file.path(proj_dir, "integrated_annotations")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

collapse_unique <- function(x) {
  x <- x[!is.na(x) & x != "-" & x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(unique(x), collapse = "|")
}
split_pipe <- function(x) {
  x <- x[!is.na(x) & x != "-" & x != ""]
  if (length(x) == 0) return(character(0))
  unlist(strsplit(x, "\\|", fixed = FALSE), use.names = FALSE)
}

# 1) mapping gene_id -> pretty_id
map <- fread(mapping_tsv)
stopifnot(all(c("gene_id","pretty_id") %in% names(map)))
map[, protein_id := sub("^lcl\\|", "", gene_id)]
map <- unique(map[, .(protein_id, pretty_id)])

# 2) InterProScan aggregate
ips <- fread(interpro_tsv, sep = "\t", header = FALSE, quote = "", fill = TRUE)
setnames(ips, c("protein_id","md5","len","analysis","sig_acc","sig_desc",
                "start","end","score","status","date","ipr","ipr_desc","go","pathway"))

ips_agg <- ips[, .(
  ips_ipr      = collapse_unique(ipr),
  ips_ipr_desc = collapse_unique(ipr_desc),
  ips_go       = collapse_unique(unique(unlist(lapply(go, split_pipe)))),
  ips_pathway  = collapse_unique(unique(unlist(lapply(pathway, split_pipe))))
), by = protein_id]

ips_pretty <- merge(map, ips_agg, by = "protein_id", all.x = TRUE)

ips_pretty_agg <- ips_pretty[, .(
  ips_ipr      = collapse_unique(unique(unlist(lapply(ips_ipr, split_pipe)))),
  ips_ipr_desc = collapse_unique(unique(unlist(lapply(ips_ipr_desc, split_pipe)))),
  ips_go       = collapse_unique(unique(unlist(lapply(ips_go, split_pipe)))),
  ips_pathway  = collapse_unique(unique(unlist(lapply(ips_pathway, split_pipe))))
), by = pretty_id]

# 3) eggNOG-mapper annotations
# IMPORTANT: emapper files have metadata lines starting with '##' and a header line starting with '#query'.
# We must drop only '##' lines, not the '#query' header.
em_lines <- readLines(emapper_ann, warn = FALSE)
em_lines2 <- em_lines[!startsWith(em_lines, "##")]

em <- fread(text = paste(em_lines2, collapse = "\n"), sep = "\t", header = TRUE, quote = "", fill = TRUE)

# detect query column
qcol <- NULL
for (cand in c("#query", "query", "Query")) {
  if (cand %in% names(em)) { qcol <- cand; break }
}
if (is.null(qcol)) stop(paste("No query column found. Columns are:", paste(names(em), collapse=", ")))

setnames(em, qcol, "query_id")
em[, protein_id := sub("^lcl\\|", "", query_id)]

keep_cols <- intersect(c(
  "protein_id",
  "seed_ortholog","evalue","score",
  "eggNOG_OGs","max_annot_lvl","COG_category","Description",
  "Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC",
  "CAZy","BiGG_Reaction","PFAMs"
), names(em))

em2 <- unique(em[, ..keep_cols])

em_pretty <- merge(map, em2, by = "protein_id", all.x = TRUE)

cols_to_collapse <- setdiff(names(em_pretty), c("pretty_id"))
em_pretty_agg <- em_pretty[, as.list(lapply(.SD, function(v) collapse_unique(unique(v)))),
                           by = pretty_id, .SDcols = cols_to_collapse]

# 4) merge
master <- merge(ips_pretty_agg, em_pretty_agg, by = "pretty_id", all = TRUE)

out_file <- file.path(out_dir, "master_annotations_pretty_id.tsv")
fwrite(master, out_file, sep = "\t", quote = FALSE, na = "")

cat("Wrote:", out_file, "\n")
cat("Rows:", nrow(master), "\n")
cat("Has IPS GO:", sum(!is.na(master$ips_go)), "\n")
if ("GOs" %in% names(master)) cat("Has eggNOG GO:", sum(!is.na(master$GOs)), "\n")
