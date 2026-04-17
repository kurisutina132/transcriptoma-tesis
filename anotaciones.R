
pacman::p_load(
  tidyverse,      # incluye ggplot2 y otros
  rio,            # importar/exportar
  here,           # localizador de ficheros
  stringr,        # trabajar con caracteres   
  scales,         # transformar números
  ggrepel,        # etiquetas colocadas inteligentemente
  gghighlight,    # resaltar una parte del gráfico
  RColorBrewer    # escalas de color
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)

KEGG_asvsco <- read.csv2("C:/Users/kuris/OneDrive/Escritorio/trancriptoma/KEGG_asvsco.csv")
View(KEGG_asvsco)
head(KEGG_asvsco)
ls()

                 
                    
as <- data.frame(KEEG_as)                    
as
dotplot(as)
library(ggplot2)


ggplot(KEEG_as, aes(x = GeneRatio, y = Description)) +
  geom_point(position = position_jitter(width = 0.1)) + # Add slight jitter for overlapping points
  labs(title = "Dotplot of KEEG_as", x = "KEEG_as values") +
  theme_classic()

KEGG_asvsco <- read.csv2("C:/Users/kuris/OneDrive/Escritorio/trancriptoma/KEGG_asvsco.csv")
head(KEGG_asvsco)
ggplot(data = KEGG_asvsco, aes(x = Description, y = count)) +
  geom_bar(stat = "identity", aes(fill = p.adjust)) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ))

a=ggplot(KEGG_asvsco, aes(x=rich_factor,y=Description)) +
  geom_point(aes(color=p.adjust, size=count)) +
  scale_color_gradientn(colours = rainbow(5)) +
  labs(
    x='Rich Factor', y=NULL,
    color='qvalue',size='Gene\nnumber'
  ) +
  theme(
    axis.title = element_text(face='bold'),
    axis.text = element_text(face='bold'),
    axis.text.y = element_text(size=15)
  )
a

  labs(y = "Description")+
  theme(axis.title.y = element_text(size = 20,
                                    color = "black",
                                    face = "bold"))


KEGG_cr <- read.csv2("C:/Users/kuris/OneDrive/Escritorio/trancriptoma/KEGG_cr.csv")
head(KEGG_cr)

b=ggplot(KEGG_cr, aes(x=rich_factor,y=Description)) +
  geom_point(aes(color=p.adjust, size=count)) +
  scale_color_gradientn(colours = rainbow(5)) +
  labs(
    x='Rich Factor', y=NULL,
    color='qvalue',size='Gene\nnumber'
  ) +
  theme(
    axis.title = element_text(face='bold'),
    #axis.text = element_text(face='bold'),
    axis.text.y = element_text(size=12)
  )
b

treatments <- read.csv2("C:/Users/kuris/OneDrive/Escritorio/trancriptoma/treatments.csv")
head(treatments)
c=ggplot(treatments, aes(x=rich_factor,y=Description)) +
  geom_point(aes(color=p.adjust, size=count)) +
  scale_color_gradientn(colours = rainbow(5)) +
  labs(
    x='Rich Factor', y=NULL,
    color='P-value',size='Gene\nnumber'
  ) +
  theme(
    axis.title = element_text(face='bold'),
    axis.text = element_text(face='bold'),
    axis.text.y = element_text(size=15)
  )

  
c
library(gridExtra)
grid.arrange(a,b,c, ncol=3)
grid.arrange(a,b, ncol=2)

###################################╔
GO_arsenic <- read.csv2("C:/Users/kuris/OneDrive/Escritorio/trancriptoma/go_arsenico.csv")
head(GO_arsenic)
str(GO_arsenic$rich_factor)
GO_arsenic$rich_factor <- as.numeric(as.character(GO_arsenic$rich_factor))
GO_arsenic$percentages <- as.numeric(as.character(GO_arsenic$percentages))


l=ggplot(GO_arsenic, aes(x = Description, y = percentages , fill = GO)) +
  geom_bar(position="stack", stat="identity",show.legend = FALSE) +
  facet_wrap(~GO, scales = "free")+
  coord_flip()+
  scale_color_gradientn(colours = rainbow(5)) +
  theme(
    axis.title = element_text(face='bold'),
    axis.text = element_text(face='bold'),
    axis.text.y = element_text(size=15))+
  scale_fill_manual(values = c("#CCEDB1", "#41B7C4"))
    # O ncol = 4
l
l+
  ylab("Percentages") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold"))
  

m=ggplot(GO_arsenic, aes(x=rich_factor,y=Description)) +
  geom_point(aes(color=p.adjust, size=Count)) +
  scale_color_gradientn(colours = rainbow(5)) +
       labs(
    x='Rich Factor', y=Description,
    color='qvalue',size='Gene\nnumber'
  ) +
  theme(
    axis.title = element_text(face='bold'),
    axis.text = element_text(face='bold'),
    axis.text.y = element_text(size=15)
  )+
  facet_wrap(~GO, strip.position = "bottom")+
  theme(axis.text = element_text(size = 14))
m
m+
  labs(y = "Description")+
  theme(axis.title.y = element_text(size = 20,
                                    color = "black",
                                    face = "bold"))

grid.arrange(l,m, a)

###########################################################################
set.seed(961)
GO_cromo <- read.csv2("C:/Users/kuris/OneDrive/Escritorio/trancriptoma/GO_cromo_2.csv")
head(GO_cromo)
View(GO_cromo)
p=ggplot(GO_cromo, aes(x = Description, y = percentages , fill = GO)) +
  geom_bar(position="stack", stat="identity",show.legend = FALSE) +
  facet_wrap(~GO, scales = "free")+
  coord_flip()+
  scale_color_gradientn(colours = rainbow(5)) +
  theme(axis.text = element_text(size = 12))+
  scale_fill_manual(values = c("#CCEDB1", "#41B7C4", "#FCFED4"))# O ncol = 4
p
p+
  ylab("Percentages") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold"))

o=ggplot(GO_cromo, aes(x=rich_factor,y=Description)) +
  geom_point(aes(color=p.adjust, size=Count)) +
  scale_color_gradientn(colours = rainbow(5)) +
  labs(
    x='Rich Factor', y=NULL,
    color='qvalue',size='Gene\nnumber'
  ) +
  theme(
    axis.title = element_text(face='bold'),
    axis.text = element_text(face='bold'),
    axis.text.y = element_text(size=15)
  )+
  facet_wrap(~GO, strip.position = "bottom")+
  theme(axis.text = element_text(size = 14))
o
o+
  labs(y = "Description")+
  theme(axis.title.y = element_text(size = 20,
                                    color = "black",
                                    face = "bold"))
grid.arrange(p,o,b, ncol=2, nrow =2, widths=c(4, 1.4), heights=c(1.4, 4))

grid.arrange(p, o, b,  ncol = 2, 
             layout_matrix = cbind(c(1,1,1), c(2,3,4)))

library(patchwork)
(p)/(o+b)

o + b
########################
GO_arsenic <- read.csv2("C:/Users/kuris/OneDrive/Escritorio/trancriptoma/GO_arsenic.csv")
head(GO_arsenic)

# Redondear valores de percentages y rich_factor
GO_arsenic <- GO_arsenic %>%
  mutate(
    percentages = round(percentages, 1),  # Redondear percentages a 1 decimal
    rich_factor = round(rich_factor, 1)   # Redondear rich_factor a 1 decimal
  )

# Gráfico de barras
l = ggplot(GO_arsenic, aes(x = Description, y = percentages, fill = GO)) +
  geom_bar(position="stack", stat="identity", show.legend = FALSE) +
  facet_wrap(~GO, scales = "free") +
  coord_flip() +
  scale_color_gradientn(colours = rainbow(5)) +
  scale_fill_manual(values = c("#CCEDB1", "#41B7C4"))  # O ncol = 4
l

# Gráfico de dispersión (rich factor vs Description)
m = ggplot(GO_arsenic, aes(x = rich_factor, y = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(colours = rainbow(5)) +
  scale_x_continuous(
    breaks = c(1.5, 2.0, 2.5, 3.0),
    labels = scales::number_format(accuracy = 0.1)  # Para formato bonito, opcional
  ) +
  labs(
    x = 'Rich Factor', y = NULL,
    color = 'qvalue', size = 'Gene\nnumber'
  ) +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text = element_text(face = 'bold'),
    axis.text.y = element_text(size = 15)
  ) +
  facet_wrap(~GO, strip.position = "bottom") +
  theme(axis.text = element_text(size = 14))

m
