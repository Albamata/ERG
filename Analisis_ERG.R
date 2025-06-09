# Script del TFG: ESTUDIO DEL MICROBIOMA DEL ESTUARIO DEL GUADALQUIVIR
# Autora: Alba Mata González
# Curso Académico: 2024/25

# Este script realiza un análisis taxonómico y funcional de los datos obtenidos 
# del estuario para evaluar la comunidad microbiana de este.

#Archivos requeridos
#OTU table
# OTUs.txt

#Taxonomy table
# taxonomy.txt

#Metadata table
#metadata.txt

#+++++++++++++++++++++++++++++++++++
# 1. CONFIGURACIÓN DEL ENTORNO-----
#+++++++++++++++++++++++++++++++++++

#Paquetes requeridos

library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(permute)
library(lattice)
library(vegan)
library(readr)
library(pheatmap)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(microbiome)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(FactoMineR)
library(factoextra)
library(viridis)
library(viridisLite)
library(readxl)
library(scales)
library(ggcorrplot)

# Directorio de trabajo
setwd("C:/Users/Usuario/OneDrive/Documentos/TFG/Estuario_Guadalquivir/Resultados_rio/resultados_Rio_11_11_24_MOTUS")

#+++++++++++++++++++++++++++++++++++
# 2. IMPORTACIÓN DE DATOS----
#+++++++++++++++++++++++++++++++++++

# Cargar tabla OTU (mOTUS)
otu <- read.table("MetaQVIR_motus_filtered.txt",
                  sep = "\t",
                  header = TRUE)

# ID OTU deben ser los nombres de las filas
rownames(otu) <- otu[,1]
otu[,1] <- NULL
#str(otu)  Para comprobar

# Guarda los rownames
rownames_otu <- rownames(otu)

# Convierte las columnas a numéricas
otu[] <- lapply(otu,
                function(x) suppressWarnings(as.numeric(as.character(x))))  

# Vuelve a asignar los rownames
rownames(otu) <- rownames_otu

# Eliminar X del nombre de las columnas
colnames(otu) <- gsub("^X",
                      "",
                      colnames(otu))

# Cargar tabla de taxonomía
tax <- read.table("MetaQVIR_motus_GTDB_taxonomy.tsv",
                  sep = "\t",
                  header = TRUE)
rownames(tax) <- tax[,1] 
tax[,1] <- NULL 
tax[tax == ""] <- NA

# Se debe convertir en matriz
taxmat <- as.matrix(tax)

# Cargar metadatos
metadata = read.table("analisis_metada.txt",
                      sep = "\t",
                      header = TRUE)
rownames(metadata) <- metadata[,1]
metadata[,1] <- NULL
metadata[metadata == ""] <- NA

#+++++++++++++++++++++++++++++++++++
# 3. CREACIÓN DE OBJETO PHYLOSEQ----
#+++++++++++++++++++++++++++++++++++

# Convirtiendo datos en un objeto phyloseq
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
sampledata = sample_data(metadata)

# Generando objeto phyloseq
physeq = phyloseq(OTU,
                  TAX,sampledata)

# Objetos usados en los distintos análisis
physeq1 = merge_phyloseq(physeq,
                         sampledata)

physeqreal <- prune_taxa(taxa_sums(physeq1) > 0,
                         physeq1)

#+++++++++++++++++++++++++++++++++++
# 4. PERSONALIZACIÓN GRÁFICOS----
#+++++++++++++++++++++++++++++++++++

#Para que todas las gráficas sean iguales
theme_set(theme_minimal() +
            theme(
              #text = element_text(family = "Helvetica", size = 12),
              axis.text.x = element_text(hjust = 1, color = "black", size = 12),
              axis.text.y = element_text(color = "black", size = 12),
              axis.title.x = element_text(color = "black", size = 12),
              axis.title.y = element_text(color = "black", size = 12),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 12)
            ))

#Para resetearlo
theme_set(theme())

#+++++++++++++++++++++++++++++++++++
# 5. ANÁLISIS TAXONÓMICO----
#+++++++++++++++++++++++++++++++++++

# Eliminamos taxones que no estan bien anotados
physeq_filtered <- subset_taxa(physeq1,
                               (Kingdom %in% c("d__Archaea",
                                               "d__Bacteria")))

# Agrupamos el nivel taxonomico y sumamos los taxones de nivel 
#inferior para evitar separaciones en el plot
physeq_grouped <- tax_glom(physeq_filtered,
                           taxrank = "Kingdom")

# Quitamos prefijos
tax_table(physeq_grouped)[, 
                          "Kingdom"] <- gsub("d__",
                                             "",
                                             tax_table(physeq_grouped)[,
                                                                       "Kingdom"])

## 5.1 ABUNDANCIA ABSOLUTA A NIVEL DE DOMINIO----

# Visualizamos la abundancia absoluta a nivel taxonómico de Dominio
plot_bar(physeq_grouped,
         fill = "Kingdom" )+
  labs(x = "Estaciones de muestreo", # Título ejes
       y = "Abundancia absoluta",
       fill = "Dominio") +
  scale_fill_manual(values = c("Bacteria" = "steelblue", # Acordar colores
                               "Archaea" = "#dadc8e")) +
  scale_y_continuous(breaks = seq(0,
                                  8000,
                                  by = 1000)) + # Escala eje y
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5)) # Localización y rotación elementos eje x

# Para conocer abundancias tanto de bacterias como arqueas:

# Extraemos la tabla de OTUs como matriz
otupor <- as(otu_table(physeq_grouped), "matrix")

# Sumamos abundancias totales para cada OTU (es decir, por fila)
abundancias_totales <- rowSums(otupor)

#print(abundancias_totales)

#      Bacterias     Arqueas
#       42574          588

# Suma total de todos los OTUs
suma_total <- sum(abundancias_totales)

# Porcentaje de cada dominio
dominios <- abundancias_totales / suma_total * 100

#print(dominios)

#  (%)    Bacterias     Arqueas
#         98.637691    1.362309

# Suma por estación de muestreo:

# sample_sums(physeqreal)
#A1   B1   B2   B3   C1   C2   C3   D1   D2   D3 V_A1 V_B2 V_C2 V_D2 
#4907 8091 4615 2919 4305 3277 1721 4703 3584 3396  334  176  691 1126 

## 5.2 ABUNDANCIA A NIVEL DE FILO----

# 5.2.1 ABUNDANCIA ABSOLUTA

# Agrupamos a nivel de filo
absolute_glom <- tax_glom(physeq = bac_arch_physeq, taxrank = "Phylum")

# Extraemos la tabla de taxonomía como data.frame
taxa_table <- as.data.frame(tax_table(absolute_glom))

# Modificamos la columna Phylum para eliminar "p__"
taxa_table$Phylum <- gsub("p__", "", taxa_table$Phylum)

# Convertimos la tabla modificada de nuevo a matriz y reasignarla
tax_table(absolute_glom) <- as.matrix(taxa_table)

# Convertimos a data.frame con psmelt()
absolute_df <- psmelt(absolute_glom)

# Ver estructura del resultado
#str(absolute_df)

#obtenemos los 10 más abundantes
top_10_filos <- phylum_abundance[order(phylum_abundance$Abundance,
                                       decreasing = TRUE),
][1:10, ]

#print(top_10_filos)
#               Phylum        Abundance
#            Proteobacteria     22834
#           Actinobacteriota     12711
#               Bacteroidota      2479
# Incongruent [d__Bacteria]      2479
#            Marinisomatota       982
#          Thermoplasmatota       544
#             Chloroflexota       427
#         Verrucomicrobiota       381
#           Patescibacteria       114
#        Desulfobacterota_D       100

# Escogemos los filos con mayor abundancia
phylum_deseado <- c("Actinobacteriota", "Bacteroidota", "Chloroflexota", 
                    "Incongruent Bacteria", "Marinisomatota", 
                    "Proteobacteria", "Thermoplasmatota", "Verrucomicrobiota")

# Asignamos colores
color_abundancias <- c("#8A2BE2", "#FFD700", "#32CD32", "mediumorchid4",
                       "#999999", "#00BFA6", "#FF1493", "#FF8C00", "#0000FF")

# Arreglamos nombres
absolute_df <- absolute_df %>%
  mutate(Phylum = gsub("\\[d__Bacteria\\]", 
                       "Bacteria",
                       Phylum)) 

# Agrupamos filos <1% abundancia
absolute_filtrado <- absolute_df %>%
  mutate(Phylum = ifelse(Phylum %in% phylum_deseado,
                         Phylum,
                         "Filo < 1% abundancia"))


abundancia_absoluta <- ggplot(absolute_filtrado,
                              aes(x = Sample,
                                  y = Abundance,
                                  fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_abundancias) +
  scale_y_continuous(breaks = seq(0, 8000, by = 1000)) +
  labs(x = "Estaciones de muestreo",
       y = "Abundancia Absoluta",
       fill = "Filo") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(vjust = 1,
                                   hjust=0.25))


# 5.2.2 ABUNDANCIA RELATIVA

percentages <- transform_sample_counts(bac_arch_physeq,
                                       function(x) x*100 / sum(x) )

# Agrupamos a nivel de filo
percentages_glom <- tax_glom(percentages,
                             taxrank = 'Phylum')
#View(percentages_glom@tax_table@.Data)

percentages_df <- psmelt(percentages_glom)
percentages_df$Phylum <- as.factor(percentages_df$Phylum)
phylum_colors_rel<- colorRampPalette(brewer.pal(12,
                                                "Paired")) (length(levels(percentages_df$Phylum)))

# Agrupamos los filos con <1% de abundancia relativa
percentages_df$Phylum <- as.character(percentages_df$Phylum) # Return the Phylum column to be of type character
percentages_df$Phylum[percentages_df$Abundance < 1] <- "Filo < 1% abundancia"

percentages_df_filtered <- percentages_df  

# Nos aseguramos de cómo se han agrupado
unique(percentages_df_filtered$Phylum)

percentages_df_filtered <- percentages_df_filtered%>%
  mutate(Phylum = gsub("p__", "", Phylum), 
         Phylum = gsub("\\[d__Bacteria\\]", # Arreglamos nomenclatura
                       "Bacteria", Phylum))

abundancia_relativa <- ggplot(data=percentages_df_filtered,
                        aes(x=Sample,
                            y=Abundance,
                            fill=Phylum))+ 
  geom_bar(aes(),
           stat="identity",
           position="stack")+
  labs(x = "Estaciones de muestreo",
       y ="Abundancia Relativa",
       fill = "Filo")+
  scale_fill_manual(values = color_abundancias)+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust=0.25),
        axis.line = element_line(colour = "black")) 

# Unimos ambas gráficas
(abundancia_absoluta | abundancia_relativa) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

#+++++++++++++++++++++++++++++++++++
# 6. MICROBIOMA BASAL----
#+++++++++++++++++++++++++++++++++++

# Cargamos objeto
coremicro <- transform(physeq1, "compositional")

# Para hacer este apartado, se hizo una previa búsqueda de los taxones
# y cuales podían componer el microbioma basal. Pero eso se puede ajustar
# según criterio propio con detection y prevalence

# Sólo se comentará el primer taxón (filo)

## 6.1. Filo----
coremicro <- aggregate_taxa(coremicro, 
                            level = "Phylum") # Agrupamos a nivel de filo

# Creamos vector lógico para filtrar y eliminar los no válidos
valid_taxa <- !grepl("Incongruent..d__Bacteria",  
                     tax_table(coremicro)[, "Phylum"])

# Eliminamos taxones no válidos
coremicro <- prune_taxa(valid_taxa,
                        coremicro)

phy_names <- gsub("p__",
                  "",
                  tax_table(coremicro)[, "Phylum"]) # Eliminamos prefijo

# Asignamos nombres únicos y compatibles con R a los taxones (filos)
taxa_names(coremicro) <- make.names(phy_names,
                                    unique = TRUE)

# Filtramos el objeto coremicro para identificar el "microbioma basal"
coremicro <- core(coremicro,
                  detection = 0.0001, # Umbral mínimo de abundancia relativa (0.0001 = 0.01%)
                  prevalence = 0.5) # Proporción mínima de muestras en las que debe aparecer (50%)

# Composición relativa, necesario para hacer el core microbiano
coremicro <- transform(coremicro,
                       "compositional") 

# Visualización
microbiome_filo <- plot_core(coremicro,
                             plot.type = "heatmap",
                             colours = c("#FFF380", "steelblue"), # Colores en gradiente
                             prevalences = seq(0.05, 1, 0.05), # Rango de prevalencias a mostrar
                             detections = 10^seq(log10(1e-4),
                                                 log10(0.1), 
                                                 length = 10)) + # Umbrales de detección
  scale_x_discrete(name = "Umbral de detección",
                   labels = function(x) scientific(as.numeric(x),
                                                   digits = 2))+ # Notación científica con dos dígitos
  scale_y_discrete(breaks = c("Bacteroidota", "Actinobacteriota", # Se muestran filos previamente visualizados 
                              "Proteobacteria"),                  # y confirmados como el microbioma basal
                   limits = c("Bacteroidota", "Actinobacteriota", 
                              "Proteobacteria")) +
  scale_fill_gradient(name = "Prevalencia (%)", # Gradiente de colores para representar la prevalencia (porcentaje)
                      low = "#FFF380", high = "steelblue",
                      breaks = seq(0, 1, by = 0.2),
                      labels = scales::percent_format(accuracy = 1)) 


## 6.2. Clase----
coremicro <- aggregate_taxa(coremicro,
                            level = "Class")
valid_taxa <- !grepl("Incongruent..d__Bacteria", 
                     tax_table(coremicro)[, "Class"])
coremicro <- prune_taxa(valid_taxa, coremicro)
cla_names <- gsub("c__", "", tax_table(coremicro)[, "Class"])
taxa_names(coremicro) <- make.names(cla_names, unique = TRUE)

coremicro <- core(coremicro, detection = 0.0001, prevalence = 0.5)
coremicro <- transform(coremicro, "compositional")

microbiome_clase <- plot_core(coremicro,
                              plot.type = "heatmap",
                              colours = c("#FFF380", "steelblue"),
                              prevalences = seq(0.05, 1, 0.05),
                              detections = 10^seq(log10(1e-4),
                                                  log10(0.1), 
                                                  length = 10)) +
  scale_x_discrete(name = "Umbral de detección",
                   labels = function(x) scientific(as.numeric(x), digits = 2))+
  scale_y_discrete(breaks = c("Acidimicrobiia", "Actinomycetia", 
                              "Alphaproteobacteria", "Gammaproteobacteria"),
                   limits = c("Acidimicrobiia", "Actinomycetia", 
                              "Alphaproteobacteria", "Gammaproteobacteria")) +
  scale_fill_gradient(name = "Prevalencia (%)",
                      low = "#FFF380", high = "steelblue",
                      breaks = seq(0, 1, by = 0.2),
                      labels = scales::percent_format(accuracy = 1)) 

## 6.3. Orden----
coremicro <- aggregate_taxa(coremicro, level = "Order")
valid_taxa <- !grepl("Incongruent..d__Bacteria",
                     tax_table(coremicro)[, "Order"])
coremicro <- prune_taxa(valid_taxa, coremicro)
or_names <- gsub("o__", "", tax_table(coremicro)[, "Order"])
taxa_names(coremicro) <- make.names(or_names, unique = TRUE)

coremicro <- core(coremicro, detection = 0.001, prevalence = 0.7)
coremicro <- transform(coremicro, "compositional")

microbiome_orden <- plot_core(coremicro,
                              plot.type = "heatmap",
                              colours = c("#FFF380", "steelblue"),
                              prevalences = seq(0.05, 1, 0.05),
                              detections = 10^seq(log10(1e-4),
                                                  log10(0.1), 
                                                  length = 10)) +
  scale_x_discrete(name = "Umbral de detección",
                   labels = function(x) scientific(as.numeric(x), digits = 2))+
  scale_y_discrete(breaks = c("Burkholderiales", "Pelagibacterales", 
                              "Nanopelagicales", "Acidimicrobiales"),
                   limits = c("Burkholderiales", "Pelagibacterales", 
                              "Nanopelagicales", "Acidimicrobiales")) +
  scale_fill_gradient(name = "Prevalencia (%)",
                      low = "#FFF380", high = "steelblue",
                      breaks = seq(0, 1, by = 0.2),
                      labels = scales::percent_format(accuracy = 1)) 

## 6.3. Familia----
coremicro <- aggregate_taxa(coremicro, level = "Family")
valid_taxa <- !grepl("Incongruent..d__Bacteria", 
                     tax_table(coremicro)[, "Family"])
coremicro <- prune_taxa(valid_taxa, coremicro)
fam_names <- gsub("f__", "", tax_table(coremicro)[, "Family"])
taxa_names(coremicro) <- make.names(fam_names, unique = TRUE)

coremicro <- core(coremicro, detection = 0.001, prevalence = 0.7)
coremicro <- transform(coremicro, "compositional")

microbiome_familia <- plot_core(coremicro,
                                plot.type = "heatmap",
                                colours = c("#FFF380", "steelblue"),
                                prevalences = seq(0.05, 1, 0.05),
                                detections = 10^seq(log10(1e-4),
                                                    log10(0.1), 
                                                    length = 10)) +
  scale_x_discrete(name = "Umbral de detección",
                   labels = function(x) scientific(as.numeric(x), digits = 2))+
  scale_y_discrete(breaks = c("Burkholderiaceae", "Pelagibacteraceae", 
                              "Methylophilaceae"),
                   limits = c("Burkholderiaceae", "Pelagibacteraceae", 
                              "Methylophilaceae")) +
  scale_fill_gradient(name = "Prevalencia (%)",
                      low = "#FFF380", high = "steelblue",
                      breaks = seq(0, 1, by = 0.2),
                      labels = scales::percent_format(accuracy = 1)) 
## 6.4. Género----
coremicro <- aggregate_taxa(coremicro, level = "Genus")
gen_names <- gsub("g__", "", tax_table(coremicro)[, "Genus"])
taxa_names(coremicro) <- make.names(gen_names, unique = TRUE)

coremicro <- core(coremicro, detection = 0.0001, prevalence = 0.7)
coremicro <- transform(coremicro, "compositional")

microbiome_genus <- plot_core(coremicro,
                              plot.type = "heatmap",
                              colours = c("#FFF380", "steelblue"),
                              prevalences = seq(0.05, 1, 0.05),
                              detections = 10^seq(log10(1e-4),
                                                  log10(0.1), 
                                                  length = 10)) +
  scale_x_discrete(name = "Umbral de detección",
                   labels = function(x) scientific(as.numeric(x), digits = 2))+
  scale_y_discrete(breaks = c("Polynucleobacter", "Fonsibacter"),
                   limits = c("Polynucleobacter", "Fonsibacter")) +
  scale_fill_gradient(name = "Prevalencia (%)",
                      low = "#FFF380", high = "steelblue",
                      breaks = seq(0, 1, by = 0.2),
                      labels = scales::percent_format(accuracy = 1)) 


#Agrupamos todos los niveles taxonómicos en la misma figura
(microbiome_filo / microbiome_clase / microbiome_orden / microbiome_familia / microbiome_genus) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

#+++++++++++++++++++++++++++++++++++
# 7. PROPORCIÓN DE VESICULACIÓN----
#+++++++++++++++++++++++++++++++++++

## 7.1. Elección de las familias con mayor tasa de vesiculación

# Agrupamos a nivel de familia
physeqfam <- tax_glom(physeq1,
                      taxrank = "Family")

#Filtramos por las estaciones "A1", "B2", "C2" y "D2"
zonas <- c("A1", "B2", "C2", "D2",
           "V_A1", "V_B2", "V_C2", "V_D2")
physeqfamfil <- subset_samples(physeqfam,
                               sample_names(physeqfam) %in% zonas)

# Extraemos tabla OTU y taxonomía
otu <- as.data.frame(otu_table(physeqfamfil))
if (!taxa_are_rows(physeqfamfil)) {
  otu <- t(otu)
}

# Añadimos columna con los nombres de familia
tax <- as.data.frame(tax_table(physeqfamfil))
otu$Family <- tax$Family

# Convertimos a formato largo
otu_long <- otu %>%
  pivot_longer(cols = -Family, names_to = "Sample", 
               values_to = "Abundance")

# Extraemos zona base (ej., A1) y tipo (ej., V para validación)
otu_long <- otu_long %>%
  mutate(Base = str_replace(Sample, "V_", ""),
         Tipo = ifelse(str_starts(Sample, "V_"), "Vesicular", "Celular"))

# Creamos tabla comparativa: una columna para cada tipo por familia y zona base
otu_wide <- otu_long %>%
  select(Family, Base, Tipo, Abundance) %>%
  pivot_wider(names_from = Tipo,
              values_from = Abundance, 
              values_fill = 0)

# Ordenamos columnas para claridad
otu_wide <- otu_wide %>%
  relocate(Celular, Vesicular, .after = Base) %>%
  filter(!(Celular == 0 & Vesicular == 0))

# Para guardar tabla
write.csv2(otu_wide, 
           file = "comparacion_familias_celular_vesicular.csv", 
           row.names = FALSE)
# Para ver cómo es graficamente, en Excel ordenamos por abundancias quedándonos con 
# las 10 familias con más de 20 conteos tanto en la fracción celular como vesicular

# Para visualizar todas las familias
otu_long_plot <- otu_wide %>%
  pivot_longer(cols = c(Celular, Vesicular),
               names_to = "Fraccion", 
               values_to = "Abundancia")

ggplot(otu_long_plot, aes(x = Family, y = Abundancia, fill = Fraccion)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Base, scales = "free_x") +
  labs(title = "Comparación de Abundancias por Familia",
       x = "Familia", y = "Abundancia") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## 7.2. Heatmap familias para fracción celular y vesicular----

# Obtenemos abundancias relativa
physeq_fam_rel <- transform_sample_counts(physeqfam,
                                          function(x) x / sum(x))

# Quitamos el prefijo "f__" de los nombres de familia
tax_table(physeq_fam_rel)[, "Family"] <- gsub("^f__", "", tax_table(physeq_fam_rel)[, "Family"])

# Filtramos por las 10 familias de interés
familia_interes <- c("S36-B12", "Methylophilaceae", "Nanopelagicaceae", 
                     "Pelagibacteraceae", "Ilumatobacteraceae", "Burkholderiaceae",
                     "Flavobacteriaceae", "HIMB59", "Limnocylindraceae",
                     "AcAMD-5")
physeq_fam_rel <- subset_taxa(physeq_fam_rel, Family %in% familia_interes)

# Filtramos por zonas
physeq_fam_rel <- subset_samples(physeq_fam_rel, sample_names(physeq_fam_rel) %in% zonas)

# Extraemos tabla de abundancias relativas como data frame
otu_relab_df <- as.data.frame(otu_table(physeq_fam_rel))

# Agregamos nombres taxonómicos 
taxa_df <- as.data.frame(tax_table(physeq_fam_rel))

# Eliminamos columnas de género y especie
taxa_df <- taxa_df[, !(colnames(taxa_df) %in% c("Species", "Genus"))]

# Asignamos la nueva tax_table al objeto phyloseq
tax_table(physeq_fam_rel) <- tax_table(as.matrix(taxa_df))

# Transponemos el data frame para que las filas sean muestras y columnas sean familias
heatmap_data <- t(otu_relab_df)

# Asignamos nombres más legibles a las columnas (familias)
colnames(heatmap_data) <- tax_table(physeq_fam_rel)[, "Family"]

# Asignamos colores y breaks a la gráfica

valor_blanco <- 0 # Para saber qué valores no tengo 
val_abund_min <- 0.0002037905 # Valor más bajo con abundancia
val_abund_max <- 0.4 # Valor más alto esperado
n_colors <- 100 # Afinamos número de pasos en el gradiente

# Creamos el breaks desde 0 hasta valor máximo
breaks <- c(valor_blanco, seq(val_abund_min, 
                              val_abund_max, 
                              length.out = n_colors + 1))

# Creamos paleta de colores
colors <- c("white",
            viridis(n_colors, # Viridis es una paleta usada 
                    begin = 0.2, # en el ámbito científico
                    end = 1,
                    direction = -1)) 

# Generamos heatmap de las familias para la fracción celular y vesicular
pheatmap(heatmap_data,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         fontsize_row = 10,
         fontsize_col = 10)

# Guardamos matriz
heatmap_datadf <- as.data.frame(heatmap_data)
write.xlsx(heatmap_datadf,
           file = "resultado_matriz.xlsx",
           rowNames = TRUE)

## 7.3. Heatmap tasa de vesiculación----

# Cargamos archivo de Excel
proporcion <- read_excel("C:/Users/Usuario/OneDrive/Documentos/TFG/Estuario_Guadalquivir/Resultados_rio/resultados_Rio_11_11_24_MOTUS/proporciones_vesiculas.xlsx")

# Aplica as.numeric a cada celda manteniendo la estructura
proporcion <- as.data.frame(lapply(proporcion, as.numeric))
proporcion <- proporcion[,-1]

# Copiamos el dataframe para hacer la transformación logarítmica
log_proporcionVS <- proporcion

# Sustituimos 0 por 1 para que no de error
log_proporcionVS[log_proporcionVS == 0] <- 1

# Aplicamos Log2
log_proporcionVS <- log2(log_proporcionVS)

# Rango de colores: queremos que blanco sea en 0
range_low <- -4
range_high <- 4
n_colores <- 200  # total de transiciones de color
n_half <- n_colores / 2  # mitad para negativo, mitad para positivo

# Parte negativa: amarillo (#ffda03) a blanco (#ffffff)
col_neg <- colorRampPalette(c("#ffda03", "#ffffff"))(n_half)

# Parte positiva: blanco (#ffffff) a verde oscuro (#006856)
col_pos <- colorRampPalette(c("#ffffff", "#00846b"))(n_half)

# Unir los colores
colores <- c(col_neg, col_pos)

# Los breaks deben tener longitud igual a length(colores) + 1
breakes <- seq(range_low, range_high, length.out = length(colores) + 1)

# Graficar
pheatmap(log_proporcionVS,
         cluster_rows = FALSE, # Para que no se clusterice de forma jerárquica
         cluster_cols = FALSE,
         color = colores, # Personalizado
         breaks = breakes, # Personalizado
         fontsize_row = 10, # Ajusta tamaño de la fuente
         fontsize_col = 10,
         legend_breaks = seq(range_low, range_high, by = 2), # Leyenda configurada para
         legend_labels = seq(range_low, range_high, by = 2)) # mostrar cortes cada dos uds de log2

#+++++++++++++++++++++++++++++++++++
# 8. DIVERSIDAD ALFA----
#+++++++++++++++++++++++++++++++++++

bacteria_archaea_physeq <- subset_taxa(physeqreal,
                                       Kingdom %in% c("d__Bacteria",
                                                      "d__Archaea"))

plot_richness(bacteria_archaea_physeq, 
              measures = c("Chao1", # Calcula índice de Chao1
                           "Shannon"), # Calcula índice de Shannon
              color = "Fraction") +
  scale_color_manual(values = c("cells" = "firebrick",
                                "vesicles" = "forestgreen"))+
  labs(x = "Zonas",
       y = "Diversidad alpha")+
  geom_point(size = 4) # Ajusta el tamaño de los puntos

#Para obtener el dato numérico
richness_values <- estimate_richness(bacteria_archaea_physeq,
                                     measures = c("Shannon",
                                                  "Chao1"))

#print(richness_values$Shannon)
#[1] 4.327588 4.550463 4.666110 5.051419 4.201012 4.554720 4.323622 4.259939 4.057454
#[10] 4.224056 3.035157 3.579817 3.693863 3.528284

#print(richness_values$Chao1)
#[1] 458.12500 670.78947 572.14706 530.93671 329.83721 336.16418 278.00000 333.38636
#[9] 297.09756 264.81395  83.21429  90.10000 132.27273 116.66667

#+++++++++++++++++++++++++++++++++++
# 9. DIVERSIDAD BETA----
#+++++++++++++++++++++++++++++++++++

meta_ordb <- ordinate(physeq = percentages,
                      method = "PCoA", 
                      distance = "bray")

# Comprobar matriz
#head(meta_ordb)

#Para cambiar el nombre de las zonas
sample_data_df <- data.frame(sample_data(physeq_phylum))

sample_data_df <- sample_data_df%>%
  mutate(Zonas = case_when(str_detect(Zona, "A") ~ "MA",
                           str_detect(Zona, "B") ~ "ZB",
                           str_detect(Zona, "C") ~ "ZM",
                           str_detect(Zona, "D") ~ "ZA",
                           TRUE ~ NA_character_  ))
sample_data(physeq_phylum) <- sample_data(sample_data_df)

# Elegimos color y forma de los datos
colores <- c("MA" = "#DF536B",
             "ZB"="#26828E",
             "ZM"="#6ADE4B",
             "ZA"="#FDE725")
forma <- c("vesicles" = 17, "cells" = 19)

# Graficamos diversidad beta
pb <- plot_ordination(physeq_phylum,
                      ordination = meta_ordb,
                      color = "Zonas",
                      shape = "Fraction") +
  scale_color_manual(values = colores,
                     breaks = c("MA", "ZB", "ZM", "ZA")) + 
  scale_shape_manual(values = forma,
                     labels = c("cells" = "Celular",
                                "vesicles" = "Vesicular")) +
  labs(shape = "Fracción") +
  geom_point(size = 4) +
  theme_classic() 

#print(pb) 

#+++++++++++++++++++++++++++++++++++
# 10. PERMANOVA----
#+++++++++++++++++++++++++++++++++++

# Combinamos los datos taxonómicos y los metadatos
physeq_permat = merge_phyloseq(physeq,  
                               sampledata)

# Eliminamos los que están mal anotados
physeq_permat <- physeq_permat %>%
  subset_taxa(!(Phylum %in% c("Incongruent [d__Bacteria]", 
                              "Not_annotated"))) 


# Comprobamos si hay algún dato como NA
#any(is.na(otu_table(physeq_permat)))  
#FALSE

# Probamos datos, creamos dataframe con la suma de los conteos
permat_df <- data.frame(sum = sample_sums(physeq_permat)) #sumo spp de cada zona

# Visualizamos la distribución de lecturas
ggplot(permat_df, aes(x = sum))+
  geom_histogram(color = "black",
                 fill = "indianred",
                 binwidth = 2500) +
  xlab("Read counts") +
  theme_minimal()

# Estadísticas básicas sobre la abundancia total del muestreo
smin <- min(sample_sums(physeq_permat)) # Mínima
#print(smin) 
#134

smean <- mean(sample_sums(physeq_permat)) # Media
#print(smean) 
#2907.786

smax <- max(sample_sums(physeq_permat)) # Máximo
#print(smax) 
#6961

# Lo hacemos a nivel de filo
physeq_permat_phylum <- physeq_permat %>%
  tax_glom(taxrank = "Phylum") %>% # Agrupamos
  transform_sample_counts(function(x) {x * 100 /sum(x)} ) # Abundancia relativa

# Calculamos datos con la distancia de Bray-Curtis
dist_bray <- distance(physeq = physeq_permat_phylum,
                      method = "bray")

physeq_permat_df <- physeq_permat_phylum %>%
  psmelt() %>% # Convertimos en dataframe
  filter(Abundance > 2) %>% # Filtramos taxones poco abundantes
  arrange(Phylum) # Ordenamos de menor a mayor

# Definimos estructura de permutación: libre entre grupos de la variable "Zona"
perm <- how(within = Within(type = "none"),
            plots = with(metadata,
                         Plots(strata = Zona,
                               type = "free")))

# Agrupamos niveles de Zona para análisis (AB vs CD)
sample_data <- sample_data(physeq_permat_phylum)
ZonaAg <- ifelse(grepl("^A|^B",
                       sample_data$Zona),
                 "AB", "CD")
# Realizamos PERMANOVA con distancia Bray-Curtis y 999 permutaciones
adonis2(dist_bray ~ ZonaAg,
        data = metadata,
        permutations = 999,
        method = "bray")

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dist_bray ~ ZonaAg, data = metadata, permutations = 999, method = "bray")
#          Df SumOfSqs  R2      F     Pr(>F)    
#Model     1  0.46689 0.57237 16.061  0.001 ***
# Residual 12  0.34883 0.42763                  
#Total    13  0.81571 1.00000                  
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#+++++++++++++++++++++++++++++++++++
# 11. CORRELACIÓN----
#+++++++++++++++++++++++++++++++++++

# Cargamos el texto donde aparecen los nutrientes, 
# propiedades fisicoquímicas y componentes biológicos
correlaciones <- read.table("correlaciones_tfg.txt",
                            sep = "\t",
                            header = TRUE)


correlaciones <- correlaciones %>%
  select(-X.1, -X.2, -Zone) %>% # Salen columnas que no interesan, borrar
  mutate(EDAR = ifelse(Station %in% c("2", "5", "6", "7", "9", "10"), 0, # Añadimos nueva columna: cerca de EDAR = 1
                       ifelse(Station %in% c("11", "13", "14", "15"), 1, NA))) %>%              #  lejos de EDAR = 0
  rename(CDOM = CDOMM, # Renombramos variables para que no de M = media
         Turb = TurbM,
         pH = pHM,
         Chl = ChlM,
         O.Sat = O.SatM,
         Sal = SalM,
         Temp = TempM)

# Una vez añadida la columna EDAR, eliminamos la columna Station
correlaciones <- correlaciones %>%
  select(- Station)

#view(correlaciones)

#Las correlaciones deben tener valores numéricos
correlaciones[-1] <- lapply(correlaciones[-1], function(x) {
  if (is.character(x)) {
    # Sustituir coma por punto
    x_clean <- gsub(",", ".", x)
    # Intentar convertir a numérico
    if (all(!is.na(as.numeric(x_clean)))) {
      as.numeric(x_clean)  # Si es posible convertir, lo convierte
    } else {
      x  # Si no, deja el texto original
    }
  } else {
    as.numeric(x)
  }
})

# Para coger solo las columnas numéricas
correlaciones_num <- correlaciones [,-1]

# Matriz p valores, pruebas añadiendo Fracción Vesicular
matriz_p <- cor_pmat(correlaciones_num)

cor(select(correlaciones, -c(X, Fraction))) %>% 
  ggcorrplot(., 
             method = "circle",
             type = "upper",
             p.mat = cor_pmat(correlaciones[,-1]),
             insig = "blank",
             hc.order = TRUE,
             show.diag = F) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.1))

# Lo mismo pero sin vesículas
correlacionesSV <- correlaciones[1:(nrow(correlaciones) - 4), ]

# Matriz correlación
matriz_corSV <- cor(select(correlacionesSV, -c(X, Fraction)), 
                    method = "spearman")
# Matriz p valores
matriz_pSV <- cor_pmat(select(correlacionesSV, -c(X, Fraction)))

# Convertir a formato largo
cor_df <- as.data.frame(as.table(matriz_corSV))
pval_df <- as.data.frame(as.table(matriz_pSV))
names(pval_df)[3] <- "p_value"

# Lo combinamos en un sólo dataframe
full_df <- left_join(cor_df, pval_df, by = c("Var1", "Var2"))
names(full_df)[3] <- "cor"

# Nos aseguramos que Var1 y Var2 tengan el mismo orden de factores
orden_vars <- unique(full_df$Var1)
full_df$Var1 <- factor(full_df$Var1, levels = orden_vars)
full_df$Var2 <- factor(full_df$Var2, levels = orden_vars)

# Incluimos triángulo inferior **con diagonal**
full_df <- full_df %>% 
  filter(as.numeric(Var1) <= as.numeric(Var2)) %>%
  mutate(is_diagonal = Var1 == Var2)

# Graficar
correlacion <- ggplot(full_df, aes(x = Var1, y = Var2)) +
  geom_point(data = full_df %>% filter(!is_diagonal), # Dibuja los valores de la correlación
             aes(x = Var1,
                 y = Var2,
                 size = ifelse(p_value <= 0.05, # Tamaño círculo según p valor
                               p_value, NA),
                 fill = cor), 
             shape = 21,
             color = "black",
             na.rm = TRUE) + # Omite NA si los hubiese
  geom_point(data = full_df %>% filter(is_diagonal), # Grafica la diagonal, da referencia
             aes(x = Var1,
                 y = Var2,
                 size = ifelse(p_value <= 0.05,
                               p_value, NA),
                 fill = cor),
             shape = 21,
             color = "red",
             stroke = 1.5, # Puntos de la diagonal con borde rojo y más grueso
             na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", # Asignamos colores,
                       mid = "white", # colorea según correlación
                       high = "red",
                       midpoint = 0,
                       breaks = c(-1, -0.5, 0, 0.5, 1), # Acotamos y señalamos
                       limits = c(-1, 1)) +
  scale_size_continuous(range = c(2, 10), # El tamaño círculo más grande cuanto más significativo
                        trans = "reverse",
                        guide = guide_legend(title = "Significancia")) +
  scale_y_discrete(limits = orden_vars,
                   position = "right") +
  scale_x_discrete(limits = rev(orden_vars),
                   position = "top") + # Matriz de tipo triangular
  theme_minimal() + # Estilo de gráfica
  coord_fixed() +
  labs(size = "Significancia", fill = "Correlación")  + 
  guides(fill = guide_colorbar(barwidth = 15, barheight = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title.position = "top")

#+++++++++++++++++++++++++++++++++++
# 12. GÉNEROS ASOCIADOS CON PATOGÉNESIS Y AGUAS RESIDUALES ----
#+++++++++++++++++++++++++++++++++++

# Ponemos abundancias relativas
relativa <- transform_sample_counts(physeq1,
                                    function(x) x*100 / sum(x))

# Agrupamos a nivel de géneros
relativa <- tax_glom(relativa,
                     taxrank = 'Genus')

# Elimina géneros de los que no tengamos abundancia
relativa <- prune_taxa(taxa_sums(relativa) >0, relativa)

# Géneros de interés (se ha hecho una previa búsqueda para confirmar su patogeneidad)
generos_conjunto <- c("Sulfuritalea", "Runella","Malikia",
                      "Hydrogenophaga", "Brachymonas", "Alicycliphilus",
                      "Azospira", "Aquirickettsiella", "Acidovorax",
                      "Stenotrophomonas", "Sphingomonas", "Pseudomonas",
                      "Novosphingobium", "Cutibacterium", "Comamonas", "Brevundimonas",
                      "Aliarcobacter")

# Extraemos OTU y taxonomía
otu_df <- as.data.frame(otu_table(relativa))
tax_df <- as.data.frame(tax_table(relativa))
tax_df$OTU <- rownames(tax_df)
otu_df$OTU <- rownames(otu_df)

# Reorganizamos en formato largo y combinamos con taxonomía
otu_long <- otu_df %>%
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance") %>%
  left_join(tax_df, by = "OTU")

# Ponemos la tabla taxonómica tambien como nombre base
otu_long <- otu_long %>%
  mutate(
    Genus = str_trim(Genus),
    Genus_base = str_replace(Genus, "^g__([A-Za-z]+).*", "\\1")
  )

otu_long$Fraccion <- ifelse(grepl("^V_", otu_long$Sample), "Vesicular", "Celular")

# Sumamos las abundancias relativas de los géneros de interés por muestra
abundancia_rel_por_muestra <- otu_long %>%
  filter(Genus_base %in% generos_conjunto) %>% 
  group_by(Sample) %>%
  summarise(Abundancia_total = sum(Abundance)) %>%
  arrange(desc(Abundancia_total))  # Opcional: ordenar de mayor a menor

# Mostrar el resultado
abundancia_rel_por_muestra

#Para asegurarnos
otu_filtrado <- otu_long %>%
  filter(Genus_base %in% generos_conjunto)

# Gráfica con todos los generos juntos

otu_conjunto <- otu_long %>%
  filter(Genus_base %in% generos_conjunto) %>%
  mutate(Fraccion = ifelse(grepl("^V_", Sample), "Vesicular", "Celular"),
         Genus_base = factor(Genus_base, levels = generos_conjunto))

# Para saber abundancias mínima, máxima y media (+ desviación típica)
media_conjunto <- otu_conjunto %>%
  group_by(Fraccion, Genus_base) %>%
  summarise(Min = min(Abundance),
            Max = max(Abundance),
            Media = mean(Abundance),
            SD = sd(Abundance),
            .groups = "drop")

#Unimos media con abundancias por muestra
conjunto_plot <- otu_conjunto %>%
  left_join(media_conjunto, by = c("Fraccion", "Genus_base")) %>%
  mutate(Desviacion = (Abundance - Media) / SD)

#Nos aseguramos que mantenga ese orden
conjunto_plot <- conjunto_plot %>%
  mutate(Genus_base = factor(Genus_base, levels = generos_conjunto))

# Filtramos por fracción
celular_plot_conjunto <- conjunto_plot %>% filter(Fraccion == "Celular")
vesicular_plot_conjunto <- conjunto_plot %>% filter(Fraccion == "Vesicular")

# Gráfico para fracción Celular
grafica_celular_conjunto <- ggplot(celular_plot_conjunto, aes(x = Sample, y = Genus_base)) +
  geom_point(aes(size = Abundance, color = Desviacion)) +
  scale_color_gradient2(
    low = "#31688E",
    mid = "gray90",
    high = "#aa3a3a",
    midpoint = 0,
    limits = c(-max(abs(celular_plot_conjunto$Desviacion)), max(abs(celular_plot_conjunto$Desviacion))),
    name = "Desviación estándar") +
  scale_size_continuous(range = c(1, 10)) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  labs(x = "Zonas",
       y = "Género",
       size = "Abundancia relativa (%)") + 
  theme(axis.text.x = element_text(hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Gráfico para fracción Vesicular
grafica_vesicular_conjunto <- ggplot(vesicular_plot_conjunto, aes(x = Sample, y = Genus_base)) +
  geom_point(aes(size = Abundance, color = Desviacion)) +
  scale_color_gradient2(
    low = "#31688E",
    mid = "gray90",
    high = "#aa3a3a",
    midpoint = 0,
    limits = c(-max(abs(vesicular_plot_conjunto$Desviacion)), max(abs(vesicular_plot_conjunto$Desviacion))),
    name = "Desviación estándar") +
  scale_size_continuous(range = c(1, 10)) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  labs(x = "Zonas",
       y = "Género",
       size = "Abundancia relativa (%)") +
  theme(axis.text.x = element_text(hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

#+++++++++++++++++++++++++++++++++++
# 13. FAPROTAX----
#+++++++++++++++++++++++++++++++++++

# Cargar archivos (versión de Windows)
# Windows
 tab1 <- read.delim("C:/Users/34644/OneDrive/Documentos/UNIVERSIDAD/4 CARRERA/TFG/Estuario_Guadalquivir/Resultados_rio/resultados_Rio_11_11_24_MOTUS/MetaQVIR_motus_filtered.txt", header = TRUE, sep = "\t", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
 tab2 <- read.delim("C:/Users/34644/OneDrive/Documentos/UNIVERSIDAD/4 CARRERA/TFG/Estuario_Guadalquivir/Resultados_rio/resultados_Rio_11_11_24_MOTUS/MetaQVIR_motus_GTDB_taxonomy.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Unimos las tablas por columna común (asumimos que se llama "X")
tabla_unida <- merge(tab1, tab2, by = "X", all.x = TRUE)

# Guardamos la tabla unida
write.table(tabla_unida, 
            "~/Documentos/TFG/Estuario_Guadalquivir/Resultados_rio/resultados_Rio_11_11_24_MOTUS/otu_table.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Establecemos nuevo directorio de trabajo
setwd("~/Documentos/TFG/Estuario_Guadalquivir/Resultados_rio/resultados_Rio_11_11_24_MOTUS/FAPROTAX")

# Cargamos los datos funcionales
functional_data <- read.table("functional_table.tsv", 
                              header = TRUE, 
                              sep = "",  # Detecta espacios
                              quote = "", 
                              fill = TRUE, 
                              comment.char = "", 
                              check.names = FALSE)

# Exploración inicial
dim(functional_data)
head(functional_data)
colnames(functional_data)

# Eliminamos la primera columna y convertimos a numérico
functional_data_num <- apply(functional_data[,-1], 2, as.numeric)
rownames(functional_data_num) <- functional_data[[1]]

# Ahora a matriz y limpiamos
functional_matrix <- as.matrix(functional_data_num)
functional_matrix[is.na(functional_matrix)] <- 0
functional_matrix <- functional_matrix[rowSums(functional_matrix) != 0, ]

# Visualización 
pheatmap(functional_matrix, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Mapa de calor de funciones ecológicas")

# Exportamos a TIFF con alta resolución
tiff("mapa_calor.tiff", width = 4000, height = 4000, res = 300)
pheatmap(functional_matrix, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Mapa de calor de funciones ecológicas")
dev.off()

# Índices de diversidad
shannon_index <- diversity(functional_matrix, index = "shannon")
simpson_index <- diversity(functional_matrix, index = "simpson")

# Creamos dataframe con índices
diversity_indices <- data.frame(
  Sample = rownames(functional_matrix),
  Shannon = shannon_index,
  Simpson = simpson_index)

# Shannon
ggplot(diversity_indices, aes(x = Sample, y = Shannon)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Índice de Diversidad de Shannon", 
       x = "Muestras", 
       y = "Índice de Shannon") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Simpson
ggplot(diversity_indices, aes(x = Sample, y = Simpson)) +
  geom_bar(stat = "identity", fill = "salmon", color = "black") +
  theme_minimal() +
  labs(title = "Índice de Diversidad de Simpson", 
       x = "Muestras", 
       y = "Índice de Simpson") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

