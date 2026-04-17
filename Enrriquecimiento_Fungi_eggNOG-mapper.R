#construir TERM2GENE (GO y KEGG) desde eggNOG

# Requiere: emapper ya leído (con columnas: #query, GOs, KEGG_ko, KEGG_Pathway)
# y map con gene_id -> pretty_id

# Tabla genes base
genes_base <- emapper %>%
  dplyr::rename(gene_id = `#query`) %>%
  dplyr::left_join(map %>% dplyr::select(gene_id, pretty_id), by = "gene_id") %>%
  dplyr::mutate(pretty_id = dplyr::coalesce(pretty_id, gene_id))

# -------- GO TERM2GENE: term = GO:xxxx ; gene = pretty_id ----------
go_term2gene <- genes_base %>%
  dplyr::select(pretty_id, GOs) %>%
  dplyr::filter(!is.na(GOs), GOs != "-") %>%
  tidyr::separate_rows(GOs, sep = ",") %>%
  dplyr::mutate(GOs = stringr::str_trim(GOs)) %>%
  dplyr::filter(stringr::str_detect(GOs, "^GO:\\d+")) %>%
  dplyr::transmute(term = GOs, gene = pretty_id) %>%
  dplyr::distinct()

# -------- KEGG KO TERM2GENE: term = Kxxxxx ; gene = pretty_id ----------
ko_term2gene <- genes_base %>%
  dplyr::select(pretty_id, KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko), KEGG_ko != "-") %>%
  tidyr::separate_rows(KEGG_ko, sep = ",") %>%
  dplyr::mutate(KEGG_ko = stringr::str_trim(KEGG_ko)) %>%
  dplyr::filter(stringr::str_detect(KEGG_ko, "^K\\d+")) %>%
  dplyr::transmute(term = KEGG_ko, gene = pretty_id) %>%
  dplyr::distinct()

# -------- KEGG Pathway TERM2GENE: term = mapXXXXX / koXXXXX / path string (depende del emapper)
keggpath_term2gene <- genes_base %>%
  dplyr::select(pretty_id, KEGG_Pathway) %>%
  dplyr::filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
  tidyr::separate_rows(KEGG_Pathway, sep = ",") %>%
  dplyr::mutate(KEGG_Pathway = stringr::str_trim(KEGG_Pathway)) %>%
  dplyr::transmute(term = KEGG_Pathway, gene = pretty_id) %>%
  dplyr::distinct()

# Universo de genes para enrich (todos los genes anotados/expresados)
universe_genes <- unique(genes_base$pretty_id)


#Separar GO en BP/MF/CC (recomendado con GO.db)
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db")

# Construir un mapping GO term -> Ontology (BP/MF/CC)
go_terms <- sort(unique(go_term2gene$term))

# GO.db devuelve ontology por GO id (BP/MF/CC)
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

# TERM2GENE por ontología
go_bp_term2gene <- go_term2gene %>% dplyr::inner_join(go_term2ont %>% dplyr::filter(ontology == "BP"), by = "term") %>% dplyr::select(term, gene)
go_mf_term2gene <- go_term2gene %>% dplyr::inner_join(go_term2ont %>% dplyr::filter(ontology == "MF"), by = "term") %>% dplyr::select(term, gene)
go_cc_term2gene <- go_term2gene %>% dplyr::inner_join(go_term2ont %>% dplyr::filter(ontology == "CC"), by = "term") %>% dplyr::select(term, gene)

#Enrichment con clusterProfiler::enricher() + tabla dinámica + dotplot

if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")

run_enricher <- function(gene_list, term2gene, term2name = NULL, universe = universe_genes,
                         p_adj_cut = 0.05) {
  if (length(gene_list) == 0) return(NULL)
  
  enr <- clusterProfiler::enricher(
    gene = gene_list,
    universe = universe,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH"
  )
  
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  
  # filtrar por padj si quieres (para mostrar)
  enr_df <- as.data.frame(enr)
  enr_df <- enr_df %>% dplyr::arrange(p.adjust, pvalue)
  
  enr
}

enrich_to_dt <- function(enr, caption) {
  if (is.null(enr)) return(NULL)
  df <- as.data.frame(enr) %>%
    dplyr::transmute(
      Term = ID,
      Description = Description,
      GeneRatio = GeneRatio,
      BgRatio = BgRatio,
      pvalue = pvalue,
      padjust = p.adjust,
      count = Count
    )
  DT::datatable(
    df,
    options = list(pageLength = 15, scrollX = TRUE),
    caption = caption
  )
}

enrich_dotplot <- function(enr, title, show_n = 15) {
  if (is.null(enr)) return(NULL)
  clusterProfiler::dotplot(enr, showCategory = show_n) + ggplot2::ggtitle(title)
}


#Iterar automáticamente por contrastes (DE + Enrichment)

# Localiza contrastes
deg_files <- list.files("deseq2_out", pattern = "^DEG_.*\\.tsv$", full.names = TRUE)

contrast_name_from_path <- function(p) {
  x <- basename(p)
  x <- sub("^DEG_", "", x)
  x <- sub("\\.tsv$", "", x)
  x
}

# lee DEG y agrega pretty_id
read_deg_with_pretty <- function(deg_path) {
  readr::read_tsv(deg_path, show_col_types = FALSE) %>%
    dplyr::left_join(map %>% dplyr::select(gene_id, pretty_id), by = "gene_id") %>%
    dplyr::mutate(pretty_id = dplyr::coalesce(pretty_id, gene_id))
}

# Ejecutar todo por contraste
results_by_contrast <- list()

# Requiere: clusterProfiler instalado/cargado
# if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")

run_enricher <- function(gene_list,
                         term2gene,
                         term2name = NULL,
                         universe = universe_genes,
                         p_adjust_method = "BH") {
  gene_list <- unique(stats::na.omit(gene_list))
  if (length(gene_list) == 0) return(NULL)
  
  # clusterProfiler espera columnas exactamente: term, gene
  if (!all(c("term", "gene") %in% names(term2gene))) {
    stop("term2gene debe tener columnas: term, gene. Columnas actuales: ", paste(names(term2gene), collapse = ", "))
  }
  if (!is.null(term2name) && !all(c("term", "description") %in% names(term2name))) {
    stop("term2name debe tener columnas: term, description. Columnas actuales: ", paste(names(term2name), collapse = ", "))
  }
  
  enr <- clusterProfiler::enricher(
    gene = gene_list,
    universe = universe,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = p_adjust_method
  )
  
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  enr
}

enrich_to_table <- function(enr) {
  if (is.null(enr)) return(dplyr::tibble())
  as.data.frame(enr) %>%
    dplyr::transmute(
      Term = ID,
      Description = Description,
      GeneRatio = GeneRatio,
      BgRatio = BgRatio,
      pvalue = pvalue,
      padjust = p.adjust,
      count = Count
    ) %>%
    dplyr::arrange(padjust, pvalue, dplyr::desc(count))
}

enrich_to_dt <- function(enr, caption) {
  df <- enrich_to_table(enr)
  if (nrow(df) == 0) {
    return(DT::datatable(df, caption = caption))
  }
  DT::datatable(
    df,
    options = list(pageLength = 15, scrollX = TRUE),
    caption = caption
  )
}

enrich_dotplot <- function(enr, title, show_n = 15) {
  if (is.null(enr)) return(NULL)
  clusterProfiler::dotplot(enr, showCategory = show_n) + ggplot2::ggtitle(title)
}

for (f in deg_files) {
  cn <- contrast_name_from_path(f)
  
  deg <- read_deg_with_pretty(f) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))
  
  # lista de genes significativos (puedes añadir |log2FC| si quieres)
  sig_genes <- deg %>%
    dplyr::filter(!is.na(padj) & padj < 0.05) %>%
    dplyr::pull(pretty_id) %>%
    unique()
  
  # GO BP/MF/CC
  enr_bp <- run_enricher(sig_genes, go_bp_term2gene, go_term2name, universe_genes)
  enr_mf <- run_enricher(sig_genes, go_mf_term2gene, go_term2name, universe_genes)
  enr_cc <- run_enricher(sig_genes, go_cc_term2gene, go_term2name, universe_genes)
  
  # KEGG: puedes elegir KO o Pathway (yo te sugiero Pathway si viene bien anotado)
  # Si KEGG_Pathway está muy vacío, usa KO.
  enr_kegg_path <- run_enricher(sig_genes, keggpath_term2gene, term2name = NULL, universe = universe_genes)
  enr_kegg_ko <- run_enricher(sig_genes, ko_term2gene, term2name = NULL, universe = universe_genes)
  
  results_by_contrast[[cn]] <- list(
    deg = deg,
    sig_genes = sig_genes,
    go_bp = enr_bp,
    go_mf = enr_mf,
    go_cc = enr_cc,
    kegg_path = enr_kegg_path,
    kegg_ko = enr_kegg_ko
  )
}

names(results_by_contrast)[1]
enrich_to_table(results_by_contrast[[1]]$go_bp) %>% head(10)

list.files("deseq2_out", pattern="^DEG_.*\\.tsv$")

#1) Crea universos específicos para KEGG (KO / Pathway)

kegg_ko_universe <- sort(unique(ko_term2gene$gene))
kegg_path_universe <- sort(unique(keggpath_term2gene$gene))

length(kegg_ko_universe)
length(kegg_path_universe)

#Modifica tu loop: KEGG con universo adecuado + fallback automático

results_by_contrast <- list()

for (f in deg_files) {
  cn <- contrast_name_from_path(f)
  
  deg <- read_deg_with_pretty(f) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))
  
  sig_genes <- deg %>%
    dplyr::filter(!is.na(padj) & padj < 0.05) %>%
    dplyr::pull(pretty_id) %>%
    unique()
  
  # GO
  enr_bp <- run_enricher(sig_genes, go_bp_term2gene, go_term2name, universe_genes)
  enr_mf <- run_enricher(sig_genes, go_mf_term2gene, go_term2name, universe_genes)
  enr_cc <- run_enricher(sig_genes, go_cc_term2gene, go_term2name, universe_genes)
  
  # KEGG: usar universo restringido a genes con anotación KEGG
  # (esto evita los mensajes y hace el test estadístico más correcto)
  sig_kegg_path <- intersect(sig_genes, kegg_path_universe)
  sig_kegg_ko <- intersect(sig_genes, kegg_ko_universe)
  
  enr_kegg_path <- run_enricher(sig_kegg_path, keggpath_term2gene, term2name = NULL, universe = kegg_path_universe)
  enr_kegg_ko <- run_enricher(sig_kegg_ko, ko_term2gene, term2name = NULL, universe = kegg_ko_universe)
  
  # escoger el "mejor" KEGG (el que tenga más términos enriquecidos)
  n_path <- if (is.null(enr_kegg_path)) 0 else nrow(as.data.frame(enr_kegg_path))
  n_ko <- if (is.null(enr_kegg_ko)) 0 else nrow(as.data.frame(enr_kegg_ko))
  kegg_mode <- if (n_path >= n_ko) "kegg_path" else "kegg_ko"
  
  results_by_contrast[[cn]] <- list(
    deg = deg,
    sig_genes = sig_genes,
    go_bp = enr_bp,
    go_mf = enr_mf,
    go_cc = enr_cc,
    kegg_path = enr_kegg_path,
    kegg_ko = enr_kegg_ko,
    kegg_mode = kegg_mode
  )
}


#Final Report” (tablas dinámicas + dotplots) por contraste

render_final_report <- function(cn, show_n = 15) {
  obj <- results_by_contrast[[cn]]
  if (is.null(obj)) stop("No existe el contraste: ", cn)
  
  cat("## Differential Expression Analysis\n")
  print(DT::datatable(
    obj$deg %>% dplyr::select(gene_id, pretty_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj),
    options = list(pageLength = 15, scrollX = TRUE),
    caption = paste0("DESeq2 Results: ", cn)
  ))
  
  cat("\n## Enrichment Analysis\n")
  
  cat("\n### Gene Ontology (GO): Biological Process (BP)\n")
  print(enrich_to_dt(obj$go_bp, paste0("GO BP: ", cn)))
  print(enrich_dotplot(obj$go_bp, paste0("Dotplot GO BP: ", cn), show_n = show_n))
  
  cat("\n### Gene Ontology (GO): Molecular Function (MF)\n")
  print(enrich_to_dt(obj$go_mf, paste0("GO MF: ", cn)))
  print(enrich_dotplot(obj$go_mf, paste0("Dotplot GO MF: ", cn), show_n = show_n))
  
  cat("\n### Gene Ontology (GO): Cellular Component (CC)\n")
  print(enrich_to_dt(obj$go_cc, paste0("GO CC: ", cn)))
  print(enrich_dotplot(obj$go_cc, paste0("Dotplot GO CC: ", cn), show_n = show_n))
  
  cat("\n### KEGG Pathway Ontology\n")
  enr_kegg <- obj[[obj$kegg_mode]]
  kegg_label <- if (obj$kegg_mode == "kegg_path") "KEGG_Pathway (eggNOG)" else "KEGG_ko (eggNOG)"
  print(enrich_to_dt(enr_kegg, paste0("KEGG: ", cn, " — ", kegg_label)))
  print(enrich_dotplot(enr_kegg, paste0("Dotplot KEGG: ", cn, " — ", kegg_label), show_n = show_n))
}


#Seleccionar contrastes (sin los _sig_...)

deg_files <- list.files("deseq2_out", pattern = "^DEG_.*\\.tsv$", full.names = TRUE)

# excluir los ya-filtrados
deg_files <- deg_files[!grepl("_sig_padj", basename(deg_files))]

contrast_name_from_path <- function(p) {
  x <- basename(p)
  x <- sub("^DEG_", "", x)
  x <- sub("\\.tsv$", "", x)
  x
}

names(deg_files) <- vapply(deg_files, contrast_name_from_path, character(1))
deg_files

#Orden del “Final Report” (Estado y Exposición
#report_order_estado <- c("T_vs_C")
report_order_exposicion <- c("As_vs_Cr", "Cr_vs_Control", "As_vs_Control")

# Validación
missing <- setdiff(c(report_order_estado, report_order_exposicion), names(deg_files))
if (length(missing) > 0) stop("Faltan estos contrastes DEG_*.tsv: ", paste(missing, collapse = ", "))


#Correr DE + Enrichment para SOLO esos contrastes

report_contrasts <- c(report_order_estado, report_order_exposicion)

results_by_contrast <- list()

for (cn in report_contrasts) {
  f <- deg_files[[cn]]
  
  deg <- read_deg_with_pretty(f) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))
  
  sig_genes <- deg %>%
    dplyr::filter(!is.na(padj) & padj < 0.05) %>%
    dplyr::pull(pretty_id) %>%
    unique()
  
  # GO
  enr_bp <- run_enricher(sig_genes, go_bp_term2gene, go_term2name, universe_genes)
  enr_mf <- run_enricher(sig_genes, go_mf_term2gene, go_term2name, universe_genes)
  enr_cc <- run_enricher(sig_genes, go_cc_term2gene, go_term2name, universe_genes)
  
  # KEGG con universo restringido
  sig_kegg_path <- intersect(sig_genes, kegg_path_universe)
  sig_kegg_ko <- intersect(sig_genes, kegg_ko_universe)
  
  enr_kegg_path <- run_enricher(sig_kegg_path, keggpath_term2gene, term2name = NULL, universe = kegg_path_universe)
  enr_kegg_ko <- run_enricher(sig_kegg_ko, ko_term2gene, term2name = NULL, universe = kegg_ko_universe)
  
  n_path <- if (is.null(enr_kegg_path)) 0 else nrow(as.data.frame(enr_kegg_path))
  n_ko <- if (is.null(enr_kegg_ko)) 0 else nrow(as.data.frame(enr_kegg_ko))
  kegg_mode <- if (n_path >= n_ko) "kegg_path" else "kegg_ko"
  
  results_by_contrast[[cn]] <- list(
    deg = deg,
    sig_genes = sig_genes,
    go_bp = enr_bp,
    go_mf = enr_mf,
    go_cc = enr_cc,
    kegg_path = enr_kegg_path,
    kegg_ko = enr_kegg_ko,
    kegg_mode = kegg_mode
  )
}


#Render “Final Report” (secciones + tablas dinámicas + dotplots)

render_final_report <- function(cn, show_n = 15) {
  obj <- results_by_contrast[[cn]]
  if (is.null(obj)) stop("No existe el contraste: ", cn)
  
  cat("\n\n### Differential Expression Analysis — DESeq2 Results\n")
  print(DT::datatable(
    obj$deg %>% dplyr::select(gene_id, pretty_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj),
    options = list(pageLength = 15, scrollX = TRUE),
    caption = paste0("DESeq2 Results: ", cn)
  ))
  
  cat("\n\n### Enrichment Analysis — Gene Ontology (GO): Biological Process (BP)\n")
  print(enrich_to_dt(obj$go_bp, paste0("GO BP: ", cn)))
  print(enrich_dotplot(obj$go_bp, paste0("Dotplot GO BP: ", cn), show_n))
  
  cat("\n\n### Enrichment Analysis — Gene Ontology (GO): Molecular Function (MF)\n")
  print(enrich_to_dt(obj$go_mf, paste0("GO MF: ", cn)))
  print(enrich_dotplot(obj$go_mf, paste0("Dotplot GO MF: ", cn), show_n))
  
  cat("\n\n### Enrichment Analysis — Gene Ontology (GO): Cellular Component (CC)\n")
  print(enrich_to_dt(obj$go_cc, paste0("GO CC: ", cn)))
  print(enrich_dotplot(obj$go_cc, paste0("Dotplot GO CC: ", cn), show_n))
  
  cat("\n\n### Enrichment Analysis — KEGG Pathway Ontology\n")
  enr_kegg <- obj[[obj$kegg_mode]]
  kegg_label <- if (obj$kegg_mode == "kegg_path") "KEGG_Pathway (eggNOG)" else "KEGG_ko (eggNOG)"
  print(enrich_to_dt(enr_kegg, paste0("KEGG: ", cn, " — ", kegg_label)))
  print(enrich_dotplot(enr_kegg, paste0("Dotplot KEGG: ", cn, " — ", kegg_label), show_n))
}

# --------- Ensamble del reporte completo (en tu Rmd) ----------
cat("# Enrichment Analysis — Final Report\n")

cat("\n## Overall Results — Estado (T/C)\n")
render_final_report("T_vs_C", show_n = 15)

cat("\n## Overall Results — Exposición\n")

cat("\n### As/Cr\n")
render_final_report("As_vs_Cr", show_n = 15)

cat("\n### Control/Cr\n")
render_final_report("Cr_vs_Control", show_n = 15)

cat("\n### Control/As\n")
render_final_report("As_vs_Control", show_n = 15)


## Exporta SOLO dotplots (GO BP/MF/CC y KEGG) por contraste a PNG y PDF
# Requiere:
# - results_by_contrast ya calculado
# - enrich_dotplot() definido
# - ggplot2 cargado
# - results_by_contrast[[cn]]$kegg_mode definido ("kegg_path" o "kegg_ko")

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)

show_n <- 15
width_in <- 8
height_in <- 6
dpi <- 300

sanitize_cn <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

export_one <- function(plot_obj, file_stub) {
  if (is.null(plot_obj)) return(invisible(FALSE))
  ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                  plot_obj, width = width_in, height = height_in, dpi = dpi)
  ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".pdf")),
                  plot_obj, width = width_in, height = height_in)
  invisible(TRUE)
}

for (cn in names(results_by_contrast)) {
  obj <- results_by_contrast[[cn]]
  cn2 <- sanitize_cn(cn)
  
  # GO BP/MF/CC
  p_bp <- enrich_dotplot(obj$go_bp, paste0("GO BP — ", cn), show_n = show_n)
  p_mf <- enrich_dotplot(obj$go_mf, paste0("GO MF — ", cn), show_n = show_n)
  p_cc <- enrich_dotplot(obj$go_cc, paste0("GO CC — ", cn), show_n = show_n)
  
  export_one(p_bp, paste0("dotplot_GO_BP_", cn2))
  export_one(p_mf, paste0("dotplot_GO_MF_", cn2))
  export_one(p_cc, paste0("dotplot_GO_CC_", cn2))
  
  # KEGG (usa el modo elegido automáticamente)
  enr_kegg <- obj[[obj$kegg_mode]]
  kegg_label <- if (obj$kegg_mode == "kegg_path") "KEGG_Pathway" else "KEGG_ko"
  p_kegg <- enrich_dotplot(enr_kegg, paste0(kegg_label, " — ", cn), show_n = show_n)
  
  export_one(p_kegg, paste0("dotplot_KEGG_", cn2))
}
