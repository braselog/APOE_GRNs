library(ArchR)
library(Seurat)
library(qs)

micro <- qread("/home/brasel/SingleCellProjects/MyProjects/pySCENIC_APOE/microglia/data/ROSMAP_microglia_Kellis.qs")
atac <- readRDS("/home/brasel/SingleCellProjects/MyProjects/pySCENIC_APOE/microglia/data/microglia.snATAC.ArchRobj.rds")

atac <- addGeneIntegrationMatrix(atac, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = micro, addToArrow = F, groupRNA = "seurat_clusters", nameCell = "predictedCell_Un", nameGroup = "predictedGroup_Un", nameScore = "predictedScore_Un")
atac_mat <- getMatrixFromProject(atac, useMatrix = "GeneScoreMatrix")
