suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

proj_dir <- "/mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"
in_dir   <- file.path(proj_dir, "gene_lists")
out_dir  <- file.path(proj_dir, "gene_lists", "transport_subcats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- c("As_vs_Control", "Cr_vs_Control", "As_vs_Cr")

read_transport <- function(comp) {
  f <- file.path(in_dir, paste0(comp, "_Transport.tsv"))
  if (!file.exists(f)) stop(paste("Missing:", f))
  fread(f)
}

write_tbl <- function(dt, comp, name) {
  f <- file.path(out_dir, paste0(comp, "_Transport_", name, ".tsv"))
  fwrite(dt, f, sep="\t", quote=FALSE, na="")
  f
}

p_transport_word <- regex("\\b(transporter|transport|permease|symporter|antiporter|channel|pump|efflux|import|export)\\b", ignore_case=TRUE)

# metal words (no 2-letter symbols)
p_metal_word <- regex(paste0("\\b(", paste(c(
  "copper","cu\\(2\\+\\)","cupro",
  "zinc","zn\\(2\\+\\)",
  "iron","fe\\(2\\+\\)","fe\\(3\\+\\)","ferric","ferrous","siderophore",
  "manganese","manganous",
  "arsenic","arsenite","arsenate","\\bars\\b",
  "chromate","chromium",
  "nickel","cobalt","cadmium"
), collapse="|"), ")\\b"), ignore_case=TRUE)

# MDR keywords: broadened
p_mdr_word <- regex(paste0("\\b(", paste(c(
  "drug","multidrug","mdr","efflux","xenobiotic","toxin",
  "pleiotropic drug resistance","pdr",
  "d(ha)?1","d(ha)?2",          # DHA1/DHA2 families often MDR (annotations vary)
  "trichothecene","tri12",
  "qac","emr"
), collapse="|"), ")\\b"), ignore_case=TRUE)

# nutrient-ish words (to avoid labeling nutrient MFS as MDR)
p_nutrient_word <- regex(paste0("\\b(", paste(c(
  "amino acid","sulfate","phosphate","sugar","hexose","glucose",
  "nitrate","ammonium","peptide","oligopeptide","urea","purine","pyrimidine","nucleobase"
), collapse="|"), ")\\b"), ignore_case=TRUE)

# PFAM signals
pf_metal <- regex(paste(c(
  "HMA","E1-E2_ATPase","Cation_ATPase","Cation_efflux","CDF","ZIP","Nramp",
  "Ftr1","Fet3","Ion_trans","Ctr","Transp_cyt_pur"
), collapse="|"), ignore_case=TRUE)

pf_abc <- regex(paste(c("ABC_", "ABC_tran", "ABC2_", "PDR_CDR"), collapse="|"), ignore_case=TRUE)
pf_mfs <- regex("MFS_", ignore_case=TRUE)

# explicit TF exclusion (any bZIP)
pf_tfs_exclude <- regex("\\bbZIP(_\\d+)?\\b", ignore_case=TRUE)

for (comp in comparisons) {
  dt <- read_transport(comp)
  dt[, text := paste(Description, PFAMs, ips_ipr_desc, ips_pathway, sep=" | ")]

  dt[, pfam_is_metal := (!is.na(PFAMs) & str_detect(PFAMs, pf_metal))]
  dt[, pfam_is_transporterish := (!is.na(PFAMs) & str_detect(PFAMs, regex("MFS_|ABC_|AA_permease|Sulfate_transp|Ion_trans|Ctr|Zip|Nramp|Cation_efflux|Transp_", ignore_case=TRUE)))]

  # Metal
  dt[, is_metal := pfam_is_metal |
                  (str_detect(text, p_metal_word) & str_detect(text, p_transport_word) & pfam_is_transporterish)]
  dt[, is_metal := is_metal & !( !is.na(PFAMs) & str_detect(PFAMs, pf_tfs_exclude) )]

  # MDR:
  # - ABC/PDR with MDR keywords OR
  # - MFS with MDR keywords (and NOT obviously nutrient unless MDR keywords present)
  dt[, is_abc_like := (!is.na(PFAMs) & str_detect(PFAMs, pf_abc))]
  dt[, is_mfs_like := (!is.na(PFAMs) & str_detect(PFAMs, pf_mfs))]

  dt[, mdr_kw_hit := str_detect(text, p_mdr_word)]

  dt[, is_mdr := (is_abc_like & mdr_kw_hit) |
                 (is_mfs_like & mdr_kw_hit) |
                 # catch cases where description explicitly says MDR but PFAM might be missing
                 (mdr_kw_hit & str_detect(text, p_transport_word) & !str_detect(text, p_nutrient_word))]

  # Nutrient: nutrient word + transport word, excluding metal/mdr
  dt[, is_nutr := (str_detect(text, p_nutrient_word) & str_detect(text, p_transport_word)) & !is_mdr & !is_metal]

  metal <- dt[is_metal == TRUE]
  mdr   <- dt[is_mdr == TRUE & is_metal == FALSE]
  nutr  <- dt[is_nutr == TRUE & is_metal == FALSE & is_mdr == FALSE]
  other <- dt[!(pretty_id %chin% c(metal$pretty_id, mdr$pretty_id, nutr$pretty_id))]

  write_tbl(metal, comp, "Metal")
  write_tbl(mdr, comp, "MDR")
  write_tbl(nutr, comp, "Nutrient")
  write_tbl(other, comp, "Other")
}

cat("Done. Wrote transport subcategories to:", out_dir, "\n")
