
set.seed(1234)
# Parameters
args <- commandArgs(T)
options(stringsAsFactors=F)
i <- as.integer(args[1])

# Load packages
library(Seurat)
library(DoubletFinder)
library(magrittr)
library(stringr)
library(purrr)
library(PRROC)

# Function
# 1.DelDoublet: Run DoubletFinder
DelDoublet <- function(SeuratObject){
  DoubletRate <- nrow(RealDataObject@meta.data) / 500 * 0.004
  sweep.res.list <- paramSweep_v3(SeuratObject, PCs = 1:16, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  homotypic.prop <- modelHomotypic(SeuratObject$seurat_clusters)   # celltype is the best choice
  nExp_poi <- round(DoubletRate * ncol(SeuratObject))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  SeuratObject <- doubletFinder_v3(SeuratObject, PCs = 1:16, pN = 0.25, pK = pK_bcmvn,
                                   nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  return(SeuratObject)
}

# 2.SeuratWorkFlow
SeuratWorkFlow <- function(RealDataObject){
  RealDataObject <- NormalizeData(RealDataObject, normalization.method = "LogNormalize", scale.factor = 1e4)
  RealDataObject <- FindVariableFeatures(RealDataObject, selection.method = "vst", nfeatures = 2000)
  RealDataObject <- ScaleData(RealDataObject)
  if (nrow(RealDataObject@meta.data) >= 50) {
    RealDataObject <- RunPCA(RealDataObject, verbose = TRUE)
    RealDataObject <- RunUMAP(RealDataObject, dims = 1:30)
    RealDataObject <- RunTSNE(RealDataObject, dims = 1:30)
    RealDataObject <- FindNeighbors(RealDataObject, dims = 1:30)
    RealDataObject <- FindClusters(RealDataObject, resolution = 1)
    return(RealDataObject)
  } else {
    npcs <- nrow(RealDataObject@meta.data) - 1
    RealDataObject <- RunPCA(RealDataObject, npcs = npcs, verbose = TRUE)
    RealDataObject <- RunUMAP(RealDataObject, dims = 1:npcs, n.neighbors = npcs)
    RealDataObject <- RunTSNE(RealDataObject, dims = 1:npcs, n.neighbors = npcs)
    RealDataObject <- FindNeighbors(RealDataObject, dims = 1:npcs)
    RealDataObject <- FindClusters(RealDataObject, resolution = 1)
    return(RealDataObject)
  }
}

# Real-Datasets ----
locs <- c('./Datasets/real_datasets/pbmc-ch.rds',
          './Datasets/real_datasets/cline-ch.rds',
          './Datasets/real_datasets/mkidney-ch.rds',
          './Datasets/real_datasets/hm-12k.rds',
          './Datasets/real_datasets/hm-6k.rds',
          './Datasets/real_datasets/pbmc-1A-dm.rds',
          './Datasets/real_datasets/pbmc-1B-dm.rds',
          './Datasets/real_datasets/pbmc-1C-dm.rds',
          './Datasets/real_datasets/pbmc-2ctrl-dm.rds',
          './Datasets/real_datasets/pbmc-2stim-dm.rds',
          './Datasets/real_datasets/J293t-dm.rds',
          './Datasets/real_datasets/pdx-MULTI.rds',
          './Datasets/real_datasets/HMEC-orig-MULTI.rds',
          './Datasets/real_datasets/HMEC-rep-MULTI.rds',
          './Datasets/real_datasets/HEK-HMEC-MULTI.rds',
          './Datasets/real_datasets/nuc-MULTI.rds')

DataID <- str_split_fixed(locs, "/", 4)[ ,4]
DataID <- gsub(".rds", "", DataID)

Results <- data.frame(DataID = DataID, AUPRC1st = NA, AUROC1st = NA,
                      AUPRC2ed = NA, AUROC2ed = NA,
                      AUPRC3rd = NA, AUROC3rd = NA,
                      AUPRC4th = NA, AUROC4th = NA,
                      AUPRC5th = NA, AUROC5th = NA,
                      AUPRC6th = NA, AUROC6th = NA,
                      AUPRC7th = NA, AUROC7th = NA,
                      AUPRC8th = NA, AUROC8th = NA,
                      AUPRC9th = NA, AUROC9th = NA,
                      AUPRC10th = NA, AUROC10th = NA,
                      AUPRCMultiCycle = NA, AUROCMultiCycle = NA)

# ----
setwd("/home/sheyong/Rprogram/SingleCellDoublet/")
RealData <- readRDS(locs[i])
setwd(paste0("./Results/DoubletFinder/", DataID[i]))

# First-time run ----
# Construct Seurat object
RealDataObject <- CreateSeuratObject(counts = RealData[[1]],
                                     project = "RealDataObject")
RealDataObject <- SeuratWorkFlow(RealDataObject)

# Calculate doublet score
RealDataObject <- DelDoublet(RealDataObject)
saveRDS(RealDataObject, paste0(DataID[i], "RealDataObject.rds"))










