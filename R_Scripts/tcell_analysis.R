setwd("C:/Users/kaaru/Box/Bigley Lab/Kaarunya Nachimuthu/Nachimuthu Projects/Nachimuthu MRV R848 Spleen scRNAseq vdj/")
options(future.globals.maxSize = 100000 * 1024^2)

library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(harmony)
library(patchwork)
library(enrichR)

mrv <- readRDS("./mrv_harmonized.rds")

DefaultAssay(mrv) <- "RNA"

data.markers <- FindAllMarkers(mrv, logfc.threshold = 0.5, only.pos = TRUE)
write.csv(data.markers, file = "CSV/AllMarkers.csv")

tcell <- subset(mrv, idents = c(2,7,10,16,32), subset = Cd19 <0.5 & Cd79a < 0.5 & Cd79b < 0.5 & Ms4a1 < 0.5 & Cd3e > 0)

DefaultAssay(tcell) <- "RNA"
tcell <- FindVariableFeatures(tcell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tcell)
tcell <- ScaleData(tcell, features = all.genes)
tcell <- RunPCA(tcell, features = VariableFeatures(tcell))
ElbowPlot(tcell, ndims = 40)

tcell <- FindNeighbors(tcell, dims = 1:10)  # Tcell-specific SNN
tcell <- FindClusters(tcell, resolution = 0.4)  # Uses new SNN
tcell <- RunUMAP(tcell, dims = 1:10)

DimPlot(tcell, label = TRUE, pt.size = 0.8) + labs(x = "UMAP 1", y = "UMAP 2")
VlnPlot(tcell, features = "Cd3e")
VlnPlot(tcell, features = "Cd4")
VlnPlot(tcell, features = "Cd8a")
VlnPlot(tcell, features = "Cd8b1")

tmark <- FindAllMarkers(tcell, logfc.threshold = 0.5, only.pos = TRUE)

#markers for annotation
zero <- FindMarkers(tcell, ident.1 = "0", logfc.threshold = 0.5, only.pos = TRUE)
zero$gene <- rownames(zero)
zero$cluster <- "0"

one <- FindMarkers(tcell, ident.1 = "1", logfc.threshold = 0.5, only.pos = TRUE)
one$gene <- rownames(one)
one$cluster <- "1"

two <- FindMarkers(tcell, ident.1 = "2", logfc.threshold = 0.5, only.pos = TRUE)
two$gene <- rownames(two)
two$cluster <- "2"

three <- FindMarkers(tcell, ident.1 = "3", logfc.threshold = 0.5, only.pos = TRUE)
three$gene <- rownames(three)
three$cluster <- "3"

four <- FindMarkers(tcell, ident.1 = "4", logfc.threshold = 0.5, only.pos = TRUE)
four$gene <- rownames(four)
four$cluster <- "4"

five <- FindMarkers(tcell, ident.1 = "5", logfc.threshold = 0.5, only.pos = TRUE)
five$gene <- rownames(five)
five$cluster <- "5"

six <- FindMarkers(tcell, ident.1 = "6", logfc.threshold = 0.5, only.pos = TRUE)
six$gene <- rownames(six)
six$cluster <- "6"

seven <- FindMarkers(tcell, ident.1 = "7", logfc.threshold = 0.5, only.pos = TRUE)
seven$gene <- rownames(seven)
seven$cluster <- "7"

eight <- FindMarkers(tcell, ident.1 = "8", logfc.threshold = 0.5, only.pos = TRUE)
eight$gene <- rownames(eight)
eight$cluster <- "8"

nine <- FindMarkers(tcell, ident.1 = "9", logfc.threshold = 0.5, only.pos = TRUE)
nine$gene <- rownames(nine)
nine$cluster <- "9"

FeaturePlot(tcell, features = c("Cd3e","Cd4","Cd8a","Cd8b1"))
FeaturePlot(tcell, features = c("Ncr1","Klrk1","Trdc","Tcrg-C1"))

markt <- rbind(zero, one, two, three, four, five, six, seven, eight, nine)
write.csv(markt, "cluster_markers.csv")

# CD4 - 0, 3, 4, 5, 6, 8
# CD8 - 1, 2, 
# GDT - 9
# NKT - 7 (contains Cd3e, Ncr1, Klrk1, Klrb1c)

cluster_labels <- c("0" = "CD4", "1" = "CD8", "2" = "CD8", "3" = "CD4", "4" = "CD4", 
                    "5" = "CD4", "6" = "CD4", "7" = "NKT","8" = "CD4", "9" = "GDT")

tcell <- RenameIdents(tcell, cluster_labels)

DimPlot(tcell, label = TRUE, pt.size = 0.5) + labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(tcell, label = TRUE, split.by = "Condition", pt.size = 1) + labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(tcell, label = TRUE, split.by = "Affstat", pt.size = 1) + labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(tcell, label = TRUE, split.by = "Sample", pt.size = 1) + labs(x = "UMAP 1", y = "UMAP 2")

tcell$annotations <- Idents(tcell)
saveRDS(tcell, "tcells_updated.rds")

# Pie chart
cluster_counts <- table(tcell$annotations)
cluster_df <- as.data.frame(cluster_counts)
colnames(cluster_df) <- c("Cluster", "Count")
cluster_df$Percentage <- (cluster_df$Count / sum(cluster_df$Count)) * 100

cluster_df

cluster_colors <- c(
  "CD8" = "#1f78b4",
  "CD4" = "#33a02c",
  "NKT" = "#e31a1c",
  "GDT" = "#ff7f00"
)

label_df <- cluster_df %>%
  arrange(desc(Cluster)) %>%
  mutate(
    ymax = cumsum(Percentage),
    ymin = ymax - Percentage,
    label_y = (ymax + ymin) / 2,
    label_text = paste0(round(Percentage, 1), "%")
  ) %>%
  ungroup()

ggplot(cluster_df, aes(x = 1, y = Percentage, fill = Cluster)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "black"
  ) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Cell composition",
    fill = "T cell type"
  ) +
  geom_label_repel(
    data = label_df,
    aes(
      x = 1.35,              # outside the pie
      y = label_y,
      label = label_text,
      fill = Cluster
    ),
    inherit.aes = FALSE,
    size = 2.8,
    show.legend = FALSE,
    label.padding = unit(0.15, "lines"),
    label.size = 0.2,
    segment.color = "grey50",
    segment.size = 0.3,
    force = 2,
    max.overlaps = Inf
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12)
  )


# Pie chart split by condition
cluster_cond_df <- tcell@meta.data %>%
  dplyr::count(Sample, annotations) %>%
  group_by(Sample) %>%
  mutate(
    Percentage = n / sum(n) * 100
  ) %>%
  ungroup()

cluster_cond_df

cluster_colors <- c(
  "CD8" = "#1f78b4",
  "CD4" = "#33a02c",
  "NKT" = "#e31a1c",
  "GDT" = "#ff7f00"
)


label_df <- cluster_cond_df %>%
  group_by(Sample) %>%
  arrange(desc(annotations)) %>%
  mutate(
    ymax = cumsum(Percentage),
    ymin = ymax - Percentage,
    label_y = (ymax + ymin) / 2,
    label_text = paste0(round(Percentage, 1), "%")
  ) %>%
  ungroup()

ggplot(cluster_cond_df, aes(x = 1, y = Percentage, fill = annotations)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "black"
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample) +
  theme_void() +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Cell composition split by sample",
    fill = "T cell type"
  ) +
  geom_label_repel(
    data = label_df,
    aes(
      x = 1.35,              # outside the pie
      y = label_y,
      label = label_text,
      fill = annotations
    ),
    inherit.aes = FALSE,
    size = 2.8,
    show.legend = FALSE,
    label.padding = unit(0.15, "lines"),
    label.size = 0.2,
    segment.color = "grey50",
    segment.size = 0.3,
    force = 2,
    max.overlaps = Inf
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12)
  )


#################################################################################################
# Pathway analysis of all clusters
cd4_de <- FindMarkers(tcell, ident.1 = "CD4", min.pct = 0.25, only.pos = TRUE)
cd8_de <- FindMarkers(tcell, ident.1 = "CD8", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(cd4_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 20, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Pathway analysis of CD4+ T cells")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(cd8_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 20, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Pathway analysis of CD8+ T cells")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}


#############################################################################################
# Pathway analysis of CD4 subset
cd4 <- subset(tcell, idents = "CD4")
Idents(cd4) <- "Sample"

#mock_Cont
mockcont_de <- FindMarkers(cd4, ident.1 = "mock_Cont", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mockcont_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 15, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  Cd4 (Mock Cont)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}


# MRV_Cont
mrvcont_de <- FindMarkers(cd4, ident.1 = "MRV_Cont", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mrvcont_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 15, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  Cd4 (MRV Cont)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}

#mock_R848
mockr848_de <- FindMarkers(cd4, ident.1 = "mock_R848", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mockr848_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 20, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  Cd4 (Mock R848)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}

#MRV_R848
mrvr848_de <- FindMarkers(cd4, ident.1 = "MRV_R848", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mrvr848_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 20, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  Cd4 (MRV R848)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}

#######################################################################################################
# Pathway analysis of CD8 subset
cd8 <- subset(tcell, idents = "CD8")
Idents(cd8) <- "Sample"

#mock_Cont
mockcont_de <- FindMarkers(cd8, ident.1 = "mock_Cont", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mockcont_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 15, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  CD8 (Mock Cont)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}


# MRV_Cont
mrvcont_de <- FindMarkers(cd8, ident.1 = "MRV_Cont", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mrvcont_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 15, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  CD8 (MRV Cont)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}

#mock_R848
mockr848_de <- FindMarkers(cd8, ident.1 = "mock_R848", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mockr848_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 20, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  CD8 (Mock R848)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}

#MRV_R848
mrvr848_de <- FindMarkers(cd8, ident.1 = "MRV_R848", min.pct = 0.25, only.pos = TRUE)

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  
  # Define the MSigDB database
  dbs <- c("MSigDB_Hallmark_2020")
  
  # Perform enrichment analysis on differentially expressed genes
  enriched <- enrichr(rownames(mrvr848_de), dbs)
  
  # Plot top 20 enriched terms, ordered by P-value
  if (enrichRLive && !is.null(enriched[["MSigDB_Hallmark_2020"]])) {
    plotEnrich(enriched[["MSigDB_Hallmark_2020"]], 
               showTerms = 20, numChar = 50, y = "Count", 
               orderBy = "P.value", title = "Enrichment analysis of  CD8 (MRV R848)")
  } else {
    message("No enrichment results found for MSigDB_Hallmark_2020.")
  }
}


# GSEA Analysis
msig <- msigdbr(species = "Mus musculus", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

gene_list <- mockcontr848_de$avg_log2FC
names(gene_list) <- rownames(mockcontr848_de)

gene_list <- sort(gene_list, decreasing = TRUE)
fgsea_res <- fgsea(
  pathways = pathways,
  stats    = gene_list,
  nperm    = 10000
)

mockcont_cd8_df <- fgsea_res %>%
  arrange(desc(NES))

top_pos <- mockcont_cd8_df %>% slice_max(NES, n = 10)
top_neg <- mockcont_cd8_df %>% slice_min(NES, n = 10)

mockcont_cd8_df_gsea <- bind_rows(top_pos, top_neg)


mockcont_cd8_df_gsea <- mockcont_cd8_df_gsea %>%
  mutate(
    pathway = pathway %>%
      str_remove("^HALLMARK_") %>%   
      str_replace_all("_", " ") %>%  
      str_to_title()                 
  )

ggplot(mockcont_cd8_df_gsea,
       aes(x = NES,
           y = reorder(pathway, NES),
           fill = padj)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_gradientn(
    colours = c("red",'darkviolet', "blue"),   
    name = "P value"
  ) +
  labs(x = "Normalized Enrichment Score",
       y = NULL,
       fill = "Adj p-value",
       title = "CD8: mock_R848 vs MRV_R848") +
  theme_classic()

































































































































