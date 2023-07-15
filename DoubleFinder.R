
# filter out doublets: DoubleFinder

# load libraries

library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubleFinder)

# create counts matrix

cts <- ReadMtx(mtx ='matrix.mtx.gz',
               features ='features.tvs.gz',
               cells = 'barcodes.tvs.gz')

cts[1:10, 1:10]

# create Seurat object

pbmc.seurat <- CreateSeuratObject(counts = cts)
str(pbmc.seurat)

# QC and Filtering
# explore QC

pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')

pbmc.seurat.filtered <- PercentFeatureSet(pbmc.seurat, pattern = '^MT-')

pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                               nFeature_RNA > 500 &
                               mitoPercent < 10)
                               
pbmc.seurat
pbmc.seurat.filtered
 
# pre-process standard workflow

pbmc.seurat.filtered <- NormalizedData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
ElbowPlot(pbmc.seurat.filtered)

pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object= pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)


## pk Identification (no ground truth)

sweep.res.list_nsclc <- paramSweep_v3(nsclc.seurat.obj, PCs = 1:20, sct = FALSE)
sweep.stats_nslc <- summarizeSweep(sweep.res.list_nsclc, Gt = FALSE)
bcmvn_nsclc <- find.pk(sweep.stats_nsclc)
        

ggplot(bcmvn_nsclc, aes(pK, BCmetric, group = 1))+
     geom_point()+
     geom_line()        
        
pK <- bcmvn_nsclc %>%
              filter(BCmetric == max(BCmetric)) %>%
              select(pK)
                     
pK <- as.numeric(as.character(pK[[1]])) 

## Homotypic Doublet Proportion Estimate------------------------------------

annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data)) 
nExp_poi.aj <- round(nExp_poi*(1-homotypic.prop))

# run doubleFinder

pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered,
                                         PCs = 1:20,
                                         pN = 0.25,
                                         pK = pK,
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)

# visualize doublets

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classification_0.25_0")
                                         
# number of singlets and doublets

table(pbmc.seurat.filtered@meta.data$DF.classification_0.25_0)     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                                        
