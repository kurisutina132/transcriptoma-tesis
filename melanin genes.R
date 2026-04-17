#Cargar eggNOG 

# encuentra .emapper.annotations o .bak
emapper_file <- list.files(
  "emapper_out",
  pattern = "emapper\\.annotations(\\.bak)?$|\\.annotations(\\.bak)?$",
  full.names = TRUE
)

if (length(emapper_file) == 0) {
  stop("No se encontró archivo de anotaciones eggNOG en emapper_out/.")
}
if (length(emapper_file) > 1) {
  message("Se encontraron múltiples archivos. Usando: ", emapper_file[1])
}
emapper_file <- emapper_file[1]

emapper <- readr::read_tsv(
  emapper_file,
  comment = "##",
  show_col_types = FALSE
)

# Master: por pretty_id (gene-IR...), listo para join con DEG
gene_annot_master <- emapper %>%
  dplyr::rename(gene_id = `#query`) %>%
  dplyr::left_join(map %>% dplyr::select(gene_id, pretty_id), by = "gene_id") %>%
  dplyr::mutate(
    pretty_id = dplyr::coalesce(pretty_id, gene_id),
    Description = dplyr::na_if(Description, "-"),
    Preferred_name = dplyr::na_if(Preferred_name, "-"),
    PFAMs = dplyr::na_if(PFAMs, "-"),
    GOs = dplyr::na_if(GOs, "-"),
    KEGG_ko = dplyr::na_if(KEGG_ko, "-"),
    KEGG_Pathway = dplyr::na_if(KEGG_Pathway, "-"),
    KEGG_Module = dplyr::na_if(KEGG_Module, "-")
  )

# chequeo rápido
message("Anotaciones cargadas: ", nrow(gene_annot_master), " filas; con pretty_id no-NA: ", sum(!is.na(gene_annot_master$pretty_id)))


#leer DEG y anexar anotación eggNOG

read_deg_annot <- function(deg_path) {
  readr::read_tsv(deg_path, show_col_types = FALSE) %>%
    dplyr::left_join(
      gene_annot_master %>%
        dplyr::select(pretty_id, Description, Preferred_name, PFAMs, GOs,
                      KEGG_ko, KEGG_Pathway, KEGG_Module, COG_category),
      by = "pretty_id"
    )
}


#metal tolerance + melanina (keywords + PFAM)

panel_patterns <- list(
  melanin = c(
    "melanin", "pigment", "pigmentation",
    # DHN-melanin keywords
    "polyketide", "PKS", "scytalone", "tetrahydroxynaphthalene", "naphthalene",
    # DOPA/oxidases
    "laccase", "multicopper", "cu-oxidase", "tyrosinase",
    # pyomelanin/tyrosine catabolism
    "homogentisate", "4-hydroxyphenylpyruvate", "hppd"
  ),
  transport = c(
    "ABC", "ABC_tran", "MFS", "MFS_1", "efflux", "transporter", "permease",
    "P-type ATPase", "ATPase", "heavy metal", "copper", "zinc", "iron", "manganese"
  ),
  oxidative = c(
    "glutathione", "transferase", "GST", "glutaredoxin", "thioredoxin",
    "peroxiredoxin", "catalase", "superoxide dismutase", "SOD",
    "oxidoreductase", "NADPH"
  ),
  vacuole_proteostasis = c(
    "vacuol", "V-ATPase", "autophag", "Hsp", "chaperone",
    "ubiquitin", "proteasome", "metallothionein"
  )
)

flag_panel_hits <- function(df) {
  text <- paste(df$Description, df$Preferred_name, df$PFAMs, df$GOs,
                df$KEGG_ko, df$KEGG_Pathway, df$KEGG_Module, sep = " | ")
  text <- tolower(ifelse(is.na(text), "", text))
  
  for (nm in names(panel_patterns)) {
    pats <- tolower(panel_patterns[[nm]])
    df[[nm]] <- Reduce(`|`, lapply(pats, function(p) stringr::str_detect(text, fixed(p))))
  }
  df$panel_any <- Reduce(`|`, df[names(panel_patterns)])
  df
}



#As vs Control y Cr vs Control (tablas listas)

read_deg_annot <- function(deg_path) {
  deg <- readr::read_tsv(deg_path, show_col_types = FALSE)
  
  if (!"gene_id" %in% names(deg)) {
    stop("El DEG no tiene columna gene_id: ", deg_path, "\nColumnas: ", paste(names(deg), collapse = ", "))
  }
  
  deg %>%
    # 1) mapear a pretty_id (gene-IR...)
    dplyr::left_join(map %>% dplyr::select(gene_id, pretty_id), by = "gene_id") %>%
    dplyr::mutate(pretty_id = dplyr::coalesce(pretty_id, gene_id)) %>%
    # 2) anexar anotación eggNOG (por gene_id; es la llave más directa)
    dplyr::left_join(
      emapper %>%
        dplyr::rename(gene_id = `#query`) %>%
        dplyr::select(gene_id, Description, Preferred_name, PFAMs, GOs,
                      KEGG_ko, KEGG_Pathway, KEGG_Module, COG_category, eggNOG_OGs, max_annot_lvl),
      by = "gene_id"
    ) %>%
    dplyr::mutate(
      Description = dplyr::na_if(Description, "-"),
      Preferred_name = dplyr::na_if(Preferred_name, "-"),
      PFAMs = dplyr::na_if(PFAMs, "-"),
      GOs = dplyr::na_if(GOs, "-"),
      KEGG_ko = dplyr::na_if(KEGG_ko, "-"),
      KEGG_Pathway = dplyr::na_if(KEGG_Pathway, "-"),
      KEGG_Module = dplyr::na_if(KEGG_Module, "-")
    )
}



deg_as <- read_deg_annot("deseq2_out/DEG_As_vs_Control.tsv") %>%
  flag_panel_hits() %>%
  dplyr::mutate(sig = !is.na(padj) & padj < 0.05)

deg_cr <- read_deg_annot("deseq2_out/DEG_Cr_vs_Control.tsv") %>%
  flag_panel_hits() %>%
  dplyr::mutate(sig = !is.na(padj) & padj < 0.05)

# Tabla: genes del panel y significativos
as_panel_sig <- deg_as %>%
  dplyr::filter(panel_any, sig) %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))

cr_panel_sig <- deg_cr %>%
  dplyr::filter(panel_any, sig) %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))

DT::datatable(
  as_panel_sig %>%
    dplyr::select(pretty_id, log2FoldChange, padj, Description, Preferred_name, PFAMs, GOs,
                  melanin, transport, oxidative, vacuole_proteostasis),
  options = list(pageLength = 15, scrollX = TRUE),
  caption = "As vs Control: genes significativos (padj<0.05) asociados a melanina / transporte / estrés oxidativo"
)

DT::datatable(
  cr_panel_sig %>%
    dplyr::select(pretty_id, log2FoldChange, padj, Description, Preferred_name, PFAMs, GOs,
                  melanin, transport, oxidative, vacuole_proteostasis),
  options = list(pageLength = 15, scrollX = TRUE),
  caption = "Cr vs Control: genes significativos (padj<0.05) asociados a melanina / transporte / estrés oxidativo"
)

#iterate_contrast_and_module_summary

# --- 0) Utilidades ---
# Detecta todos los DEG_*vs*.tsv disponibles
deg_files <- list.files("deseq2_out", pattern = "^DEG_.*\\.tsv$", full.names = TRUE)

if (length(deg_files) == 0) stop("No encontré archivos dese q2_out/DEG_*.tsv")

# Nombre humano del contraste a partir del filename
contrast_name_from_path <- function(p) {
  x <- basename(p)
  x <- sub("^DEG_", "", x)
  x <- sub("\\.tsv$", "", x)
  x
}

# Lee DEG + une map + une emapper (usa gene_id como llave)
read_deg_annot <- function(deg_path) {
  deg <- readr::read_tsv(deg_path, show_col_types = FALSE)
  if (!"gene_id" %in% names(deg)) {
    stop("El DEG no tiene columna gene_id: ", deg_path, "\nColumnas: ", paste(names(deg), collapse = ", "))
  }
  
  deg %>%
    dplyr::left_join(map %>% dplyr::select(gene_id, pretty_id), by = "gene_id") %>%
    dplyr::mutate(pretty_id = dplyr::coalesce(pretty_id, gene_id)) %>%
    dplyr::left_join(
      emapper %>%
        dplyr::rename(gene_id = `#query`) %>%
        dplyr::select(gene_id, Description, Preferred_name, PFAMs, GOs,
                      KEGG_ko, KEGG_Pathway, KEGG_Module,
                      COG_category, eggNOG_OGs, max_annot_lvl),
      by = "gene_id"
    ) %>%
    dplyr::mutate(
      Description = dplyr::na_if(Description, "-"),
      Preferred_name = dplyr::na_if(Preferred_name, "-"),
      PFAMs = dplyr::na_if(PFAMs, "-"),
      GOs = dplyr::na_if(GOs, "-"),
      KEGG_ko = dplyr::na_if(KEGG_ko, "-"),
      KEGG_Pathway = dplyr::na_if(KEGG_Pathway, "-"),
      KEGG_Module = dplyr::na_if(KEGG_Module, "-")
    )
}

# Módulos dirigidos
panel_patterns <- list(
  melanin = c(
    "melanin", "pigment", "pigmentation",
    "polyketide", "PKS", "scytalone", "tetrahydroxynaphthalene", "naphthalene",
    "laccase", "multicopper", "cu-oxidase", "tyrosinase",
    "homogentisate", "4-hydroxyphenylpyruvate", "hppd"
  ),
  transport = c(
    "ABC", "ABC_tran", "MFS", "MFS_1", "efflux", "transporter", "permease",
    "P-type ATPase", "heavy metal", "copper", "zinc", "iron", "manganese"
  ),
  oxidative = c(
    "glutathione", "transferase", "GST", "glutaredoxin", "thioredoxin",
    "peroxiredoxin", "catalase", "superoxide dismutase", "SOD",
    "oxidoreductase", "NADPH"
  ),
  vacuole_proteostasis = c(
    "vacuol", "V-ATPase", "autophag", "Hsp", "chaperone",
    "ubiquitin", "proteasome", "metallothionein"
  )
)

flag_panel_hits <- function(df) {
  text <- paste(df$Description, df$Preferred_name, df$PFAMs, df$GOs,
                df$KEGG_ko, df$KEGG_Pathway, df$KEGG_Module, sep = " | ")
  text <- tolower(ifelse(is.na(text), "", text))
  
  for (nm in names(panel_patterns)) {
    pats <- tolower(panel_patterns[[nm]])
    df[[nm]] <- Reduce(`|`, lapply(pats, function(p) stringr::str_detect(text, fixed(p))))
  }
  df$panel_any <- Reduce(`|`, df[names(panel_patterns)])
  df
}

# Conteos por módulo (significativos)
module_counts <- function(df, padj_cut = 0.05) {
  df %>%
    dplyr::mutate(sig = !is.na(padj) & padj < padj_cut) %>%
    dplyr::summarise(
      n_sig = sum(sig, na.rm = TRUE),
      melanin = sum(sig & melanin, na.rm = TRUE),
      transport = sum(sig & transport, na.rm = TRUE),
      oxidative = sum(sig & oxidative, na.rm = TRUE),
      vacuole_proteostasis = sum(sig & vacuole_proteostasis, na.rm = TRUE)
    )
}

# --- 1) Generar tabla por contraste y resumen comparativo ---
deg_tables <- list()
summary_tbl <- list()

for (f in deg_files) {
  cname <- contrast_name_from_path(f)
  
  df <- read_deg_annot(f) %>%
    flag_panel_hits() %>%
    dplyr::mutate(contrast = cname)
  
  deg_tables[[cname]] <- df
  
  summary_tbl[[cname]] <- module_counts(df) %>%
    dplyr::mutate(contrast = cname) %>%
    dplyr::select(contrast, dplyr::everything())
}

summary_tbl <- dplyr::bind_rows(summary_tbl) %>%
  dplyr::arrange(contrast)

# --- 2) Output: tablas (DT) por contraste ---
# En Rmd puedes imprimir esto en una sección:
render_contrast_table <- function(cname, padj_cut = 0.05) {
  df <- deg_tables[[cname]] %>%
    dplyr::mutate(sig = !is.na(padj) & padj < padj_cut) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))
  
  DT::datatable(
    df %>%
      dplyr::select(gene_id, pretty_id, baseMean, log2FoldChange, padj,
                    Description, Preferred_name, PFAMs, GOs,
                    melanin, transport, oxidative, vacuole_proteostasis),
    options = list(pageLength = 15, scrollX = TRUE),
    caption = paste0("Tabla anotada por contraste: ", cname, " (módulos marcados por keywords/PFAM/GO)")
  )
}

# Ejemplo de uso en el Rmd:
# render_contrast_table("As_vs_Control")
# render_contrast_table("Cr_vs_Control")

# --- 3) Resumen comparativo (As vs C vs Cr vs C) ---
# Esto te da una tabla compacta:
DT::datatable(
  summary_tbl,
  options = list(pageLength = 50, scrollX = TRUE),
  caption = "Resumen comparativo: conteos de DEG significativos (padj<0.05) por módulo"
)

# Si quieres filtrar solo los contrastes clave:
summary_key <- summary_tbl %>%
  dplyr::filter(contrast %in% c("As_vs_Control", "Cr_vs_Control", "As_vs_Cr"))

summary_key


#Barras: número de DEG up/down por contraste (visión global)

plot_deg_counts <- function(deg_tables, padj_cut = 0.05, lfc_cut = 0) {
  counts <- lapply(names(deg_tables), function(cn) {
    df <- deg_tables[[cn]]
    sig <- !is.na(df$padj) & df$padj < padj_cut
    up <- sig & df$log2FoldChange >  lfc_cut
    down <- sig & df$log2FoldChange < -lfc_cut
    data.frame(
      contrast = cn,
      direction = c("Up", "Down"),
      n = c(sum(up, na.rm = TRUE), sum(down, na.rm = TRUE))
    )
  }) |> dplyr::bind_rows()
  
  ggplot2::ggplot(counts, ggplot2::aes(x = contrast, y = n, fill = direction)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Número de DEG (padj<0.05)", fill = NULL,
                  title = "DEG Up/Down por contraste") +
    ggplot2::theme_bw()
}

plot_deg_counts(deg_tables)


#“Heatmap” (o tile plot) del resumen por módulo y contraste (tu objetivo principal)
#(melanina/transport/oxidative/vacuole)

plot_module_tile <- function(summary_tbl) {
  long <- summary_tbl %>%
    tidyr::pivot_longer(
      cols = c(melanin, transport, oxidative, vacuole_proteostasis),
      names_to = "module",
      values_to = "n_sig_module"
    )
  
  ggplot2::ggplot(long, ggplot2::aes(x = module, y = contrast, fill = n_sig_module)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
    ggplot2::labs(x = NULL, y = NULL, fill = "# DEG", title = "DEG significativos por módulo y contraste") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

plot_module_tile(summary_tbl)


#Volcano comparativo (uno por contraste) o uno interactivo

plot_volcano_panel <- function(df, contrast_title, padj_cut = 0.05) {
  df <- df %>%
    dplyr::mutate(
      sig = !is.na(padj) & padj < padj_cut,
      y = -log10(padj),
      panel = panel_any
    )
  
  ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = y)) +
    ggplot2::geom_point(color = "grey75", alpha = 0.5, size = 1) +
    ggplot2::geom_point(
      data = dplyr::filter(df, sig & panel),
      ggplot2::aes(color = melanin | transport | oxidative | vacuole_proteostasis),
      alpha = 0.9, size = 1.3
    ) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = 2) +
    ggplot2::labs(title = paste0("Volcano: ", contrast_title),
                  x = "log2FC", y = "-log10(padj)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

plot_volcano_panel(deg_tables[["As_vs_Control"]], "As vs Control")
plot_volcano_panel(deg_tables[["Cr_vs_Control"]], "Cr vs Control")


#Venn/UpSet entre contrastes (intersección de DEGs)
#1A) Preparar sets (DEG significativos) + tabla de intersección

padj_cut <- 0.05

deg_as_sig <- deg_tables[["As_vs_Control"]] %>%
  dplyr::filter(!is.na(padj) & padj < padj_cut)

deg_cr_sig <- deg_tables[["Cr_vs_Control"]] %>%
  dplyr::filter(!is.na(padj) & padj < padj_cut)

set_as <- unique(deg_as_sig$pretty_id)
set_cr <- unique(deg_cr_sig$pretty_id)

both <- intersect(set_as, set_cr)
as_only <- setdiff(set_as, set_cr)
cr_only <- setdiff(set_cr, set_as)

overlap_counts <- tibble::tibble(
  group = c("As_only", "Cr_only", "Both"),
  n = c(length(as_only), length(cr_only), length(both))
)

overlap_counts


install.packages("ggvenn")

venn_list <- list(
  "As vs Control" = set_as,
  "Cr vs Control" = set_cr
)

p_venn <- ggvenn::ggvenn(
  venn_list,
  fill_color = c("#1b9e77", "#d95f02"),
  stroke_size = 0.6,
  set_name_size = 4
) +
  ggplot2::ggtitle("Intersección de DEG (padj<0.05): As vs Control vs Cr vs Control")

p_venn


install.packages("ComplexUpset")

# construir tabla binaria (genes x sets)
all_genes <- sort(unique(c(set_as, set_cr)))

upset_df <- tibble::tibble(
  pretty_id = all_genes,
  As_vs_Control = pretty_id %in% set_as,
  Cr_vs_Control = pretty_id %in% set_cr
)

p_upset <- ComplexUpset::upset(
  upset_df,
  intersect = c("As_vs_Control", "Cr_vs_Control"),
  name = "DEG (padj<0.05)"
) +
  ggplot2::ggtitle("UpSet: intersección de DEG (As vs Control, Cr vs Control)")

p_upset


#As vs Control y Cr vs Control
plot_deg_counts_key <- function(deg_tables, padj_cut = 0.05, lfc_cut = 0,
                                contrasts = c("As_vs_Control", "Cr_vs_Control")) {
  counts <- lapply(contrasts, function(cn) {
    df <- deg_tables[[cn]]
    sig <- !is.na(df$padj) & df$padj < padj_cut
    up <- sig & df$log2FoldChange >  lfc_cut
    down <- sig & df$log2FoldChange < -lfc_cut
    data.frame(
      contrast = cn,
      direction = c("Up", "Down"),
      n = c(sum(up, na.rm = TRUE), sum(down, na.rm = TRUE))
    )
  }) |> dplyr::bind_rows()
  
  counts$contrast <- factor(counts$contrast, levels = contrasts)
  counts$direction <- factor(counts$direction, levels = c("Up", "Down"))
  
  p <- ggplot2::ggplot(counts, ggplot2::aes(x = contrast, y = n, fill = direction)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(Up = "#2b8cbe", Down = "#e34a33")) +
    ggplot2::labs(
      x = NULL, y = "Número de DEG (padj<0.05)",
      fill = NULL,
      title = "DEG por contraste",
      subtitle = "Conteos de genes diferencialmente expresados (Up/Down)"
    ) +
    ggplot2::theme_bw(base_size = 12)
  
  p
}

p_counts <- plot_deg_counts_key(deg_tables)
p_counts


plot_module_tile_key <- function(summary_tbl, contrasts = c("As_vs_Control", "Cr_vs_Control")) {
  long <- summary_tbl %>%
    dplyr::filter(contrast %in% contrasts) %>%
    tidyr::pivot_longer(
      cols = c(melanin, transport, oxidative, vacuole_proteostasis),
      names_to = "module",
      values_to = "n_sig_module"
    )
  
  long$contrast <- factor(long$contrast, levels = contrasts)
  long$module <- factor(long$module, levels = c("melanin", "transport", "oxidative", "vacuole_proteostasis"))
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = module, y = contrast, fill = n_sig_module)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = n_sig_module), color = "black", size = 4) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
    ggplot2::labs(
      x = NULL, y = NULL, fill = "# DEG",
      title = "DEG significativos por módulo",
      subtitle = "Conteo (padj<0.05) por contraste"
    ) +
    ggplot2::theme_bw(base_size = 12)
  
  p
}

p_mod <- plot_module_tile_key(summary_tbl)
p_mod

