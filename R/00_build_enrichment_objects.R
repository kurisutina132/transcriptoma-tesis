# 00_build_enrichment_objects.R
# Construye TERM2GENE (GO/KEGG) + enrichment por contraste y guarda un RDS.
# Ejecutar desde: C:/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(GO.db)
})

project_root <- "C:/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma"
setwd(project_root)

# ---------- inputs ----------
map_file <- "deseq2_out/gene_id_mapping.tsv"
stopifnot(file.exists(map_file))
map <- readr::read_tsv(map_file, show_col_types = FALSE)

emapper_file <- list.files(
  "emapper_out",
  pattern = "emapper\\.annotations(\\.bak)?$|\\.annotations(\\.bak)?$",
  full.names = TRUE
)
if (length(emapper_file) == 0) stop("No se encontró archivo eggNOG annotations en emapper_out/")
emapper_file <- emapper_file[1]
emapper <- readr::read_tsv(emapper_file, comment = "##", show_col_types = FALSE)

deg_files <- list.files("deseq2_out", pattern = "^DEG_.*\\.tsv$", full.names = TRUE)
deg_files <- deg_files[!grepl("_sig_padj", basename(deg_files))]

contrast_name_from_path <- function(p) {
  x <- basename(p)
  x <- sub("^DEG_", "", x)
  x <- sub("\\.tsv$", "", x)
  x
}
names(deg_files) <- vapply(deg_files, contrast_name_from_path, character(1))

# Orden del reporte
report_order_estado <- c("T_vs_C")
report_order_exposicion <- c("As_vs_Cr", "Cr_vs_Control", "As_vs_Control")
report_contrasts <- c(report_order_estado, report_order_exposicion)

missing <- setdiff(report_contrasts, names(deg_files))
if (length(missing) > 0) stop("Faltan contrastes DEG_*.tsv: ", paste(missing, collapse = ", "))

# ---------- genes base ----------
genes_base <- emapper %>%
  dplyr::rename(gene_id = `#query`) %>%
  dplyr::left_join(map %>% dplyr::select(gene_id, pretty_id), by = "gene_id") %>%
  dplyr::mutate(pretty_id = dplyr::coalesce(pretty_id, gene_id))

# ---------- TERM2GENE ----------
go_term2gene <- genes_base %>%
  dplyr::select(pretty_id, GOs) %>%
  dplyr::filter(!is.na(GOs), GOs != "-") %>%
  tidyr::separate_rows(GOs, sep = ",") %>%
  dplyr::mutate(GOs = stringr::str_trim(GOs)) %>%
  dplyr::filter(stringr::str_detect(GOs, "^GO:\\d+")) %>%
  dplyr::transmute(term = GOs, gene = pretty_id) %>%
  dplyr::distinct()

ko_term2gene <- genes_base %>%
  dplyr::select(pretty_id, KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko), KEGG_ko != "-") %>%
  tidyr::separate_rows(KEGG_ko, sep = ",") %>%
  dplyr::mutate(KEGG_ko = stringr::str_trim(KEGG_ko)) %>%
  dplyr::filter(stringr::str_detect(KEGG_ko, "^K\\d+")) %>%
  dplyr::transmute(term = KEGG_ko, gene = pretty_id) %>%
  dplyr::distinct()

keggpath_term2gene <- genes_base %>%
  dplyr::select(pretty_id, KEGG_Pathway) %>%
  dplyr::filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
  tidyr::separate_rows(KEGG_Pathway, sep = ",") %>%
  dplyr::mutate(KEGG_Pathway = stringr::str_trim(KEGG_Pathway)) %>%
  dplyr::transmute(term = KEGG_Pathway, gene = pretty_id) %>%
  dplyr::distinct()

# Universos
universe_genes <- sort(unique(genes_base$pretty_id))
kegg_ko_universe <- sort(unique(ko_term2gene$gene))
kegg_path_universe <- sort(unique(keggpath_term2gene$gene))

# ---------- GO split + names ----------
go_terms <- sort(unique(go_term2gene$term))

ont_vec <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = go_terms,
  columns = c("ONTOLOGY", "TERM"),
  keytype = "GOID"
)

go_term2name <- ont_vec %>%
  dplyr::transmute(term = GOID, description = TERM) %>%
  dplyr::distinct()

go_term2ont <- ont_vec %>%
  dplyr::transmute(term = GOID, ontology = ONTOLOGY) %>%
  dplyr::distinct()

go_bp_term2gene <- go_term2gene %>%
  dplyr::inner_join(go_term2ont %>% dplyr::filter(ontology == "BP"), by = "term") %>%
  dplyr::select(term, gene)

go_mf_term2gene <- go_term2gene %>%
  dplyr::inner_join(go_term2ont %>% dplyr::filter(ontology == "MF"), by = "term") %>%
  dplyr::select(term, gene)

go_cc_term2gene <- go_term2gene %>%
  dplyr::inner_join(go_term2ont %>% dplyr::filter(ontology == "CC"), by = "term") %>%
  dplyr::select(term, gene)

# ---------- helpers ----------
read_deg_with_pretty <- function(deg_path, map) {
  readr::read_tsv(deg_path, show_col_types = FALSE) %>%
    dplyr::left_join(map %>% dplyr::select(gene_id, pretty_id), by = "gene_id") %>%
    dplyr::mutate(pretty_id = dplyr::coalesce(pretty_id, gene_id))
}

run_enricher <- function(gene_list, term2gene, term2name = NULL, universe, p_adjust_method = "BH") {
  gene_list <- unique(stats::na.omit(gene_list))
  if (length(gene_list) == 0) return(NULL)
  
  clusterProfiler::enricher(
    gene = gene_list,
    universe = universe,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = p_adjust_method
  )
}

# ---------- compute per contrast ----------
results_by_contrast <- list()

for (cn in report_contrasts) {
  deg <- read_deg_with_pretty(deg_files[[cn]], map) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))
  
  sig_genes <- deg %>%
    dplyr::filter(!is.na(padj) & padj < 0.05) %>%
    dplyr::pull(pretty_id) %>%
    unique()
  
  enr_bp <- run_enricher(sig_genes, go_bp_term2gene, go_term2name, universe = universe_genes)
  enr_mf <- run_enricher(sig_genes, go_mf_term2gene, go_term2name, universe = universe_genes)
  enr_cc <- run_enricher(sig_genes, go_cc_term2gene, go_term2name, universe = universe_genes)
  
  sig_kegg_path <- intersect(sig_genes, kegg_path_universe)
  sig_kegg_ko <- intersect(sig_genes, kegg_ko_universe)
  
  enr_kegg_path <- run_enricher(sig_kegg_path, keggpath_term2gene, term2name = NULL, universe = kegg_path_universe)
  enr_kegg_ko <- run_enricher(sig_kegg_ko, ko_term2gene, term2name = NULL, universe = kegg_ko_universe)
  
  n_path <- if (is.null(enr_kegg_path)) 0 else nrow(as.data.frame(enr_kegg_path))
  n_ko <- if (is.null(enr_kegg_ko)) 0 else nrow(as.data.frame(enr_kegg_ko))
  kegg_mode <- if (n_path >= n_ko) "kegg_path" else "kegg_ko"
  
  # normaliza: si no hay resultados, guarda NULL (para que UI no rompa)
  if (!is.null(enr_bp) && nrow(as.data.frame(enr_bp)) == 0) enr_bp <- NULL
  if (!is.null(enr_mf) && nrow(as.data.frame(enr_mf)) == 0) enr_mf <- NULL
  if (!is.null(enr_cc) && nrow(as.data.frame(enr_cc)) == 0) enr_cc <- NULL
  if (!is.null(enr_kegg_path) && nrow(as.data.frame(enr_kegg_path)) == 0) enr_kegg_path <- NULL
  if (!is.null(enr_kegg_ko) && nrow(as.data.frame(enr_kegg_ko)) == 0) enr_kegg_ko <- NULL
  
  results_by_contrast[[cn]] <- list(
    deg = deg,
    sig_genes = sig_genes,
    go_bp = enr_bp, go_mf = enr_mf, go_cc = enr_cc,
    kegg_path = enr_kegg_path, kegg_ko = enr_kegg_ko,
    kegg_mode = kegg_mode
  )
}

# ---------- save ----------
dir.create("r_objects", showWarnings = FALSE)
saveRDS(
  list(
    created_at = Sys.time(),
    report_contrasts = report_contrasts,
    results_by_contrast = results_by_contrast,
    # opcional: guardar term2gene/universos para reproducibilidad
    universe_genes = universe_genes,
    kegg_ko_universe = kegg_ko_universe,
    kegg_path_universe = kegg_path_universe
  ),
  file = "r_objects/enrichment_results_by_contrast.rds"
)

cat("OK: guardado en r_objects/enrichment_results_by_contrast.rds\n")
source("R/00_build_enrichment_objects.R")
