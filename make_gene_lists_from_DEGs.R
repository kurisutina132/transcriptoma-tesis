suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

proj_dir <- "/mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"
in_dir   <- file.path(proj_dir, "deseq2_out", "annotated")
out_dir  <- file.path(proj_dir, "gene_lists")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- c(
  As_vs_Control = file.path(in_dir, "DEG_As_vs_Control_sig_padj0.05.tsv"),
  Cr_vs_Control = file.path(in_dir, "DEG_Cr_vs_Control_sig_padj0.05.tsv"),
  As_vs_Cr      = file.path(in_dir, "DEG_As_vs_Cr_sig_padj0.05.tsv")
)

# ---------- helpers ----------
split_multi <- function(x) {
  if (is.na(x) || x == "" || x == "-") return(character(0))
  x <- gsub("\\s+", "", x)
  unlist(strsplit(x, "[\\|,]"), use.names = FALSE)
}
has_any_go <- function(go_str, targets) {
  gos <- split_multi(go_str)
  any(gos %in% targets)
}
get_go_unified <- function(dt) {
  ips <- if ("ips_go" %in% names(dt)) dt$ips_go else NA_character_
  eg  <- if ("GOs" %in% names(dt)) dt$GOs else NA_character_

  out <- character(nrow(dt))
  for (i in seq_len(nrow(dt))) {
    gos <- c(split_multi(ips[i]), split_multi(eg[i]))
    gos <- unique(gos[nzchar(gos)])
    out[i] <- if (length(gos)) paste(gos, collapse=",") else NA_character_
  }
  out
}
txt_all <- function(dt) paste(dt$Description, dt$PFAMs, dt$ips_ipr_desc, dt$ips_pathway, sep=" | ")

write_list <- function(dt, out_file) {
  out <- dt[, .(
    pretty_id, gene_id, log2FoldChange, padj,
    Description, PFAMs, ips_ipr, ips_ipr_desc, ips_pathway, go_all
  )]
  fwrite(out, out_file, sep="\t", quote=FALSE, na="")
}

# ---------- DDR ----------
GO_DDR <- c("GO:0006281","GO:0006974","GO:0006302","GO:0006310","GO:0000723","GO:0000724")
KW_DDR <- regex(paste0("\\b(", paste(c(
  "dna repair","rad51","rad52","rad54","mre11","rad50","nbs1","xrs2",
  "ku70","ku80","ligase iv","checkpoint","atr","atm",
  "mismatch repair","mlh","msh","nucleotide excision","base excision",
  "photolyase","recq","ap endonuclease","exonuclease","homologous recombination"
), collapse="|"), ")\\b"), ignore_case = TRUE)

# ---------- Transport (strict) ----------
KW_TRANSPORT_STRONG <- regex(paste0("\\b(", paste(c(
  "transporter","permease","symporter","antiporter","efflux","import","export","uptake",
  "abc transporter","mfs","major facilitator","solute carrier","ion channel"
), collapse="|"), ")\\b"), ignore_case=TRUE)

PFAM_TRANSPORT <- regex(paste(c(
  "MFS_", "Aa_trans", "Sugar_tr", "Na_H_Exchanger",
  "ZIP", "Nramp", "Ftr1", "Fet3", "CDF", "Cation_efflux",
  "AcrB", "RND", "Ammonium_transp"
), collapse="|"), ignore_case=TRUE)

GO_TRANSPORT_STRICT <- c("GO:0005215","GO:0022857","GO:0015075")

# ---------- ETC ----------
GO_ETC <- c("GO:0022900","GO:0006119","GO:0006120","GO:0006121","GO:0006122")
KW_ETC <- regex(paste0("\\b(", paste(c(
  "electron transport","oxidative phosphorylation","respiratory chain",
  "cytochrome c oxidase","nadh dehydrogenase","succinate dehydrogenase",
  "ubiquinol","bc1","complex i","complex ii","complex iii","complex iv"
), collapse="|"), ")\\b"), ignore_case=TRUE)

# ---------- ATPases broad (reference bucket) ----------
GO_ATPASE <- c("GO:0016887","GO:0042626")
KW_ATPASE_BROAD <- regex(paste0("\\b(", paste(c(
  "atpase","atp synthase","atp-binding cassette"
), collapse="|"), ")\\b"), ignore_case=TRUE)

# ---------- ABC transporters ----------
PFAM_ABC <- regex(paste(c(
  "ABC_", "ABC_tran", "ABC2_", "ABC_membrane", "PDR_CDR", "Pdr5"
), collapse="|"), ignore_case=TRUE)
KW_ABC <- regex("\\babc transporter\\b|\\batp-binding cassette\\b", ignore_case=TRUE)

# ---------- PVFF ATPases (STRICT) ----------
# Keywords: specific to PVFF pumps and ATP synthase/V-ATPase
KW_PVFF <- regex(paste0("\\b(", paste(c(
  # P-type
  "p-type atpase","p type atpase","e1-e2 atpase","e1e2 atpase",
  "cation-transporting atpase","plasma membrane h\\+\\-atpase",
  "sarcoplasmic reticulum ca\\(2\\+\\) atpase","serca",

  # V-type
  "v-type atpase","v-?atpase","vacuolar atpase","vacuolar h\\+\\-atpase",
  "vacuolar atp synthase","vacuolar atpase subunit",

  # F-type / ATP synthase
  "f-type atpase","f\\-type atp synthase","atp synthase subunit",
  "fof1","f1fo","f0f1"
), collapse="|"), ")\\b"), ignore_case=TRUE)

# PFAM: allow partial matches capturing ATP-synt_A/B/C, etc. + known P-type/V-type labels you observed
PFAM_PVFF <- regex(paste(c(
  # P-type
  "P_type_ATPase","P-type_ATPase","E1E2_ATPase","E1-E2_ATPase","Cation_ATPase","HMA",

  # V-ATPase / assembly factors
  "V_ATPase","V-ATPase","Vha","Vma","Vma12",

  # ATP synthase / FoF1
  "ATP-synt","ATP_synt","ATP_synth","ATPsyn","AtpF","AtpI","OSCP","ATPase_I"
), collapse="|"), ignore_case=TRUE)

summary <- list()

for (comp in names(inputs)) {
  f <- inputs[[comp]]
  if (!file.exists(f)) stop(paste("Missing:", f))

  dt <- fread(f)
  dt[, go_all := get_go_unified(.SD)]
  dt[, text_all := txt_all(.SD)]

  # DDR / Transport / ETC
  dt[, is_DDR := mapply(has_any_go, go_all, MoreArgs=list(targets=GO_DDR)) | str_detect(text_all, KW_DDR)]

  dt[, is_Transport := mapply(has_any_go, go_all, MoreArgs=list(targets=GO_TRANSPORT_STRICT)) |
                      str_detect(text_all, KW_TRANSPORT_STRONG) |
                      (!is.na(PFAMs) & str_detect(PFAMs, PFAM_TRANSPORT)) ]

  dt[, is_ETC := mapply(has_any_go, go_all, MoreArgs=list(targets=GO_ETC)) | str_detect(text_all, KW_ETC)]

  # Broad ATPases (reference)
  dt[, is_ATPases_broad := mapply(has_any_go, go_all, MoreArgs=list(targets=GO_ATPASE)) | str_detect(text_all, KW_ATPASE_BROAD)]

  # ABC
  dt[, is_ABC := str_detect(text_all, KW_ABC) | (!is.na(PFAMs) & str_detect(PFAMs, PFAM_ABC))]

  # PVFF strict: must match PVFF keyword OR PVFF PFAM; exclude ABC
  dt[, pvff_kw_hit := str_detect(text_all, KW_PVFF)]
  dt[, pvff_pfam_hit := (!is.na(PFAMs) & str_detect(PFAMs, PFAM_PVFF))]
  dt[, is_ATPases_PVFF := (pvff_kw_hit | pvff_pfam_hit) & !is_ABC]

  # PVFF dropped: broad ATPases but not PVFF strict (and not ABC)
  dt[, is_PVFF_dropped := is_ATPases_broad & !is_ATPases_PVFF & !is_ABC]

  # output files
  f_ddr    <- file.path(out_dir, paste0(comp, "_DDR.tsv"))
  f_tr     <- file.path(out_dir, paste0(comp, "_Transport.tsv"))
  f_etc    <- file.path(out_dir, paste0(comp, "_ETC.tsv"))
  f_atp    <- file.path(out_dir, paste0(comp, "_ATPases.tsv"))
  f_pvff   <- file.path(out_dir, paste0(comp, "_ATPases_PVFF.tsv"))
  f_abc    <- file.path(out_dir, paste0(comp, "_ABC_transporters.tsv"))
  f_drop   <- file.path(out_dir, paste0(comp, "_PVFF_dropped.tsv"))

  write_list(dt[is_DDR == TRUE], f_ddr)
  write_list(dt[is_Transport == TRUE], f_tr)
  write_list(dt[is_ETC == TRUE], f_etc)
  write_list(dt[is_ATPases_broad == TRUE], f_atp)
  write_list(dt[is_ATPases_PVFF == TRUE], f_pvff)
  write_list(dt[is_ABC == TRUE], f_abc)
  write_list(dt[is_PVFF_dropped == TRUE], f_drop)

  summary[[comp]] <- data.table(
    comparison = comp,
    n_total = nrow(dt),
    n_DDR = sum(dt$is_DDR, na.rm=TRUE),
    n_Transport = sum(dt$is_Transport, na.rm=TRUE),
    n_ETC = sum(dt$is_ETC, na.rm=TRUE),
    n_ATPases_broad = sum(dt$is_ATPases_broad, na.rm=TRUE),
    n_ATPases_PVFF = sum(dt$is_ATPases_PVFF, na.rm=TRUE),
    n_ABC = sum(dt$is_ABC, na.rm=TRUE),
    n_PVFF_dropped = sum(dt$is_PVFF_dropped, na.rm=TRUE)
  )
}

summary_dt <- rbindlist(summary, use.names=TRUE, fill=TRUE)
fwrite(summary_dt, file.path(out_dir, "summary.tsv"), sep="\t", quote=FALSE, na="")
cat("Done. Wrote gene lists to:", out_dir, "\n")
print(summary_dt)
