# DoubletFinder-

DoubletFinder is an R package that predicts doublets in single-cell RNA sequencing data.

DoubletFinder Overview

DoubletFinder can be broken up into 4 steps:

(1) Generate artificial doublets from existing scRNA-seq data

(2) Pre-process merged real-artificial data

(3) Perform PCA and use the PC distance matrix to find each cell's proportion of artificial k nearest neighbors (pANN)

(4) Rank order and threshold pANN values according to the expected number of doublets

alternativetext

DoubletFinder takes the following arguments:

![image](https://github.com/sukirtipriya/DoubletFinder-/assets/88479900/4c1aec23-3cda-4f01-b06a-fa60d05c2208)


seu ~ This is a fully-processed Seurat object (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been run).

PCs ~ The number of statistically-significant principal components, specified as a range (e.g., PCs = 1:10)

pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 25%, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).

pK ~ This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values should be estimated using the strategy described below.

nExp ~ This defines the pANN threshold used to make final doublet/singlet predictions. This value can best be

estimated from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated proportion of homotypic doublets.
