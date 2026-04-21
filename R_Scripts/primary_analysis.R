setwd("C:/Users/kaaru/Box/Bigley Lab/Kaarunya Nachimuthu/Nachimuthu Projects/Nachimuthu MRV R848 Spleen scRNAseq vdj/")
options(future.globals.maxSize = 100000 * 1024^2)

library(enrichR)
library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(harmony)

mrv <- readRDS("./RDS Files/MRV R848 scRNAseq cellcycle_corrected.rds")

DefaultAssay(mrv) <- "SCT"

mito_genes <- grep(pattern = "^mt-", x = rownames(mrv), value = TRUE)
globin_genes <- grep(pattern = c("^Hb"), x = rownames(mrv), value = TRUE)
mouse_genes <- c("Mouse1", "Mouse2", "Mouse3")
genes_to_remove <- c(mito_genes, globin_genes, mouse_genes)

counts <- GetAssayData(mrv, assay = "SCT")
counts <- counts[-(which(rownames(counts) %in% c(genes_to_remove))),]
mrv <- subset(mrv, features = rownames(counts))

mrv <- RunHarmony(mrv, group.by.vars = "Sample", project.dim = FALSE)
mrv <- FindNeighbors(mrv, reduction = "harmony", dims = 1:30, graph.name = "harmony_snn")
mrv <- FindClusters(mrv, resolution = 0.8, algorithm = 2, graph = "harmony_snn")
mrv <- RunUMAP(mrv, reduction = "harmony", dims = 1:30)

DimPlot(mrv, group.by = "seurat_clusters", label = TRUE, pt.size = 0.7, label.size = 2.5, raster = FALSE)
DimPlot(mrv, group.by = "Sample", split.by = "Sample", label = FALSE, pt.size = 0.9) + NoLegend() + ggtitle('UMAP split by sample')
DimPlot(mrv, group.by = "seurat_clusters", split.by = "Condition", label = TRUE, pt.size = 0.9, label.size = 2.5)
ggsave("mrv_clusters.png", plot = p, width = 15, height = 15, dpi = 300)

#----------------------------------------------------------------------------------------------------------
genes_t <- c("Cd3e","Tcf7","Lef1","Sell","Cd4","Il2ra","Il7r","Cd8a","Cd8b1","Tbx21",
             "Ifng","Gata3","Runx2","Rorc","Ccr6","Anxa1","Foxp3","Il32","Ifit3","Mx1",
             "Gzmk","Ccl5","Klrb1","Pask","Prf1","Gzma","Tox2","Foxb","Cd7","Ncr",
             "Nkg7","Ctla4","Xcl1","Mki67")

all_genes <- rownames(mrv)
valid_genes <- genes_t[genes_t %in% all_genes]

for (gene in valid_genes) {
  featureplot <- FeaturePlot(mrv,features = gene, cols = c("lightgrey","blue"), pt.size = 0.9)
  file_name <- paste0(gene, "_feature_plot_Tcell.png")
  ggsave(file_name, plot = featureplot, width = 6, height = 5.5)
}

#----------------------------------------------------------------------------------------------------------
genes_b <- c("Cd19", "Cd22", "Cd79a", "Igha1", "Ly86", "Ms4a1", "Ighd", "Ighg1", "Mzb1", "Vpreb3",
             "Cd79b", "Ctla4", "Ptprj", "Bhlhe41", "Zbtb20", "Plac8", "Jchain", "Irf4", "Cd70", 
             "Cd38", "Cxcr3", "Mki67")

all_genes <- rownames(mrv)
valid_genes <- genes_b[genes_b %in% all_genes]

for (gene in valid_genes) {
  featureplot <- FeaturePlot(mrv,features = gene, cols = c("lightgrey","blue"), pt.size = 0.9)
  file_name <- paste0(gene, "_feature_plot_Bcell.png")
  ggsave(file_name, plot = featureplot, width = 6, height = 5.5)
}

#----------------------------------------------------------------------------------------------------------
mrv_markers <- FindAllMarkers(mrv, logfc.threshold = 0.5, only.pos = TRUE, assay = "RNA")

clusters <- as.character(0:45)
for (cluster in clusters) {
  deg <- FindMarkers(mrv, ident.1 = cluster,
                     only.pos = TRUE,
                     logfc.threshold = 0.5,
                     min.pct = 0.25,
                     assay = "RNA")
  
  write.csv(deg, paste0("cluster_", cluster, "_markers.csv"))
}































#----------------------------------------------------------------------------------------------------------
new.cluster.ids <- c('B cell','B cell','T cell','B cell','Proliferative cells','Erythrocyte',
                     'B cell','T cell','B cell','B cell','T cell','Proliferative cells','Macrophage',
                     'Erythrocyte','Erythrocyte','B cell','T cell','Erythrocyte','Erythrocyte','Dendritic cells',
                     'Monocyte','Neutrophils','Proliferative cells','Erythrocyte','NK Cells','B cell','Erythrocyte',
                     'Monocyte','Erythrocyte','Endothelial cells','Endothelial cells','B cell','T cell','Endothelial cells',
                     'Fibroblasts','B cell','NK Cells','Mast cells','B cell','Dendritic cells','Dendritic cells','Neutrophils',
                     'Dendritic cells','B cell','B cell','Dendritic cells')
names(new.cluster.ids) <- levels(mrv)
mrv <- RenameIdents(mrv, new.cluster.ids)

DimPlot(mrv, raster = FALSE, label = TRUE)
p <- DimPlot(mrv, raster = FALSE, label = TRUE, split.by = 'Sample')
ggsave(p, filename = 'mrv_clusters.png', height = 14, width = 18)






