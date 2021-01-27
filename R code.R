library(tidyverse)
library(Seurat)
library(harmony)
library(SCINA)


#Remove all metadata from samples (only one dataset is shown)
Lungs <- DietSeurat(Lungs, counts = TRUE, data = TRUE, assays = 'RNA')
Lungs@meta.data$nCount_SCT <- NULL
Lungs@meta.data$nFeature_SCT <- NULL
Lungs@meta.data$SCT_snn_res.0.8 <- NULL
Lungs@meta.data$seurat_clusters <- NULL
Lungs@meta.data$cell_labels <- NULL
Lungs$organ <- 'Lungs'
write_rds(Lungs, 'Lungs.rds')

#Merge datasets
Endothelial <- merge(Adipose, y= c(Colon, Liver, Lungs, Mesenteric, Nigra1, Nigra2, Nigra5, Nigra6, Nigra7, Pancreas, Prostate, Prostate1, Prostate2, Testis1, Testis2, Testis3, Testis4, Testis5))

#Standard Seurat workflow
Idents(Endothelial) <- 'organ'
Endothelial <- NormalizeData(Endothelial)
Endothelial <- FindVariableFeatures(Endothelial, selection.method = "vst", nfeatures = 3000)
Endothelial <- ScaleData(Endothelial)
Endothelial <- RunPCA(Endothelial, features = VariableFeatures(object = Endothelial))

#Run jackstraw
Endothelial <- JackStraw(Endothelial, num.replicate = 200, verbose = TRUE)
Endothelial <- ScoreJackStraw(Endothelial, reduction = 'pca', dims = 1:50, do.plot = TRUE)
JackStrawPlot(Endothelial, dims = 1:10)

#Run Harmony for integration
Endothelial <- RunHarmony(Endothelial,
                          group.by.vars = 'orig.ident',
                          max.iter.harmony = 100,
                          theta = 0,
                          plot_convergence = TRUE)

Endothelial <- RunUMAP(Endothelial, reduction = 'harmony' , slot = 'counts' , dims = 1:30, verbose = FALSE)
Endothelial <- RunTSNE(Endothelial,reduction = 'harmony')
Endothelial <- FindNeighbors(Endothelial, reduction = 'harmony' , dims = 1:30, verbose = FALSE)
Endothelial <- FindClusters(Endothelial, resolution = 0.2, verbose = FALSE)
DimPlot(Endothelial, reduction = 'umap', label = TRUE)
VlnPlot(Endothelial, features = c('PECAM1', 'FLT1', 'CLDN5', 'VWF', 'CLEC14A', 'ICAM2'), ncol = 3)

#Run SCINA
Endothelial_Markers <- list(c('CLDN5', 'VWF', 'FLT1', 'PECAM1'))
names(Endothelial_Markers) <- 'Endothelial cells'

SCINA_results <- SCINA(Endothelial@assays$RNA@data,
                       Endothelial_Markers,
                       max_iter = 2000, 
                       convergence_n = 100, 
                       convergence_rate = 1, 
                       sensitivity_cutoff = 0.5, 
                       rm_overlap=FALSE, 
                       allow_unknown=TRUE)

Endothelial$cell_types <- SCINA_results$cell_labels

