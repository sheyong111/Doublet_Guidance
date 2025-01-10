setwd("/home/sheyong/Rprogram/SingleCellDoublet/")

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
    RealDataObject <- FindNeighbors(RealDataObject, dims = 1:30)
    RealDataObject <- FindClusters(RealDataObject, resolution = 1)
    return(RealDataObject)
  } else {
    npcs <- nrow(RealDataObject@meta.data) - 1
    RealDataObject <- RunPCA(RealDataObject, npcs = npcs, verbose = TRUE)
    RealDataObject <- RunUMAP(RealDataObject, dims = 1:npcs, n.neighbors = npcs)
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
RealData <- readRDS(locs[i])
dir.create(paste0("./Results/DoubletFinder/", DataID[i]))
setwd(paste0("./Results/DoubletFinder/", DataID[i]))

# First-time run ----
# Construct Seurat object
RealDataObject <- CreateSeuratObject(counts = RealData[[1]],
                                     project = "RealDataObject")
RealDataObject <- SeuratWorkFlow(RealDataObject)

# Calculate doublet score
RealDataObject <- DelDoublet(RealDataObject)

# Calculate AUPRC and AUROC
label <- ifelse(RealData[[2]] == 'doublet', 1, 0)
score <- RealDataObject@meta.data[ ,ncol(RealDataObject@meta.data)]
score <- ifelse(score == 'Doublet', 1, 0)
saveRDS(score, paste0(DataID[i], "_1st_.rds"))

fg <- score[label==1]
bg <- score[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 2] <- pr$auc.integral
Results[i, 3] <- roc$auc


# Second-time run ----
names(RealDataObject@meta.data)[RealDataObject@meta.data %>% ncol()] <- "Doubletstate1st"
RealDataObject2ed <- subset(RealDataObject, Doubletstate1st == "Singlet")
RealDataObject2ed <- SeuratWorkFlow(RealDataObject2ed)
RealDataObject2ed <- DelDoublet(RealDataObject2ed)

score2ed <- score
score2ed[match(rownames(RealDataObject2ed@meta.data)[RealDataObject2ed@meta.data[ ,ncol(RealDataObject2ed@meta.data)] == "Doublet"],
            colnames(RealData[[1]]))] <- 1
saveRDS(score2ed, paste0(DataID[i], "_2ed_.rds"))

fg <- score2ed[label==1]
bg <- score2ed[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 4] <- pr$auc.integral
Results[i, 5] <- roc$auc


# Third-time run ----
names(RealDataObject2ed@meta.data)[RealDataObject2ed@meta.data %>% ncol()] <- "Doubletstate2ed"
RealDataObject3rd <- subset(RealDataObject2ed, Doubletstate2ed == "Singlet")
RealDataObject3rd <- SeuratWorkFlow(RealDataObject3rd)
RealDataObject3rd <- DelDoublet(RealDataObject3rd)

score3ed <- score2ed
score3ed[match(rownames(RealDataObject3rd@meta.data)[RealDataObject3rd@meta.data[ ,ncol(RealDataObject3rd@meta.data)] == "Doublet"],
               colnames(RealData[[1]]))] <- 1
saveRDS(score3ed, paste0(DataID[i], "_3rd_.rds"))

fg <- score3ed[label==1]
bg <- score3ed[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 6] <- pr$auc.integral
Results[i, 7] <- roc$auc


# fourth-time run ----
names(RealDataObject3rd@meta.data)[RealDataObject3rd@meta.data %>% ncol()] <- "Doubletstate3rd"
RealDataObject4th <- subset(RealDataObject3rd, Doubletstate3rd == "Singlet")
RealDataObject4th <- SeuratWorkFlow(RealDataObject4th)
RealDataObject4th <- DelDoublet(RealDataObject4th)

score4th <- score3ed
score4th[match(rownames(RealDataObject4th@meta.data)[RealDataObject4th@meta.data[ ,ncol(RealDataObject4th@meta.data)] == "Doublet"],
               colnames(RealData[[1]]))] <- 1
saveRDS(score4th, paste0(DataID[i], "_4th_.rds"))

fg <- score4th[label==1]
bg <- score4th[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 8] <- pr$auc.integral
Results[i, 9] <- roc$auc


# fifth-time run ----
names(RealDataObject4th@meta.data)[RealDataObject4th@meta.data %>% ncol()] <- "Doubletstate4th"
RealDataObject5th <- subset(RealDataObject4th, Doubletstate4th == "Singlet")
RealDataObject5th <- SeuratWorkFlow(RealDataObject5th)
RealDataObject5th <- DelDoublet(RealDataObject5th)

score5th <- score4th
score5th[match(rownames(RealDataObject5th@meta.data)[RealDataObject5th@meta.data[ ,ncol(RealDataObject5th@meta.data)] == "Doublet"],
               colnames(RealData[[1]]))] <- 1
saveRDS(score5th, paste0(DataID[i], "_5th_.rds"))

fg <- score5th[label==1]
bg <- score5th[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 10] <- pr$auc.integral
Results[i, 11] <- roc$auc


# sixth-time run ----
names(RealDataObject5th@meta.data)[RealDataObject5th@meta.data %>% ncol()] <- "Doubletstate5th"
RealDataObject6th <- subset(RealDataObject5th, Doubletstate5th == "Singlet")
RealDataObject6th <- SeuratWorkFlow(RealDataObject6th)
RealDataObject6th <- DelDoublet(RealDataObject6th)

score6th <- score5th
score6th[match(rownames(RealDataObject6th@meta.data)[RealDataObject6th@meta.data[ ,ncol(RealDataObject6th@meta.data)] == "Doublet"],
               colnames(RealData[[1]]))] <- 1
saveRDS(score6th, paste0(DataID[i], "_6th_.rds"))

fg <- score6th[label==1]
bg <- score6th[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 12] <- pr$auc.integral
Results[i, 13] <- roc$auc


# seventh-time run ----
names(RealDataObject6th@meta.data)[RealDataObject6th@meta.data %>% ncol()] <- "Doubletstate6th"
RealDataObject7th <- subset(RealDataObject6th, Doubletstate6th == "Singlet")
RealDataObject7th <- SeuratWorkFlow(RealDataObject7th)
RealDataObject7th <- DelDoublet(RealDataObject7th)

score7th <- score6th
score7th[match(rownames(RealDataObject7th@meta.data)[RealDataObject7th@meta.data[ ,ncol(RealDataObject7th@meta.data)] == "Doublet"],
               colnames(RealData[[1]]))] <- 1
saveRDS(score7th, paste0(DataID[i], "_7th_.rds"))

fg <- score7th[label==1]
bg <- score7th[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 14] <- pr$auc.integral
Results[i, 15] <- roc$auc


# eightth-time run ----
names(RealDataObject7th@meta.data)[RealDataObject7th@meta.data %>% ncol()] <- "Doubletstate7th"
RealDataObject8th <- subset(RealDataObject7th, Doubletstate7th == "Singlet")
RealDataObject8th <- SeuratWorkFlow(RealDataObject8th)
RealDataObject8th <- DelDoublet(RealDataObject8th)

score8th <- score7th
score8th[match(rownames(RealDataObject8th@meta.data)[RealDataObject8th@meta.data[ ,ncol(RealDataObject8th@meta.data)] == "Doublet"],
               colnames(RealData[[1]]))] <- 1
saveRDS(score8th, paste0(DataID[i], "_8th_.rds"))

fg <- score8th[label==1]
bg <- score8th[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 16] <- pr$auc.integral
Results[i, 17] <- roc$auc


# nineth-time run ----
names(RealDataObject8th@meta.data)[RealDataObject8th@meta.data %>% ncol()] <- "Doubletstate8th"
RealDataObject9th <- subset(RealDataObject8th, Doubletstate8th == "Singlet")
RealDataObject9th <- SeuratWorkFlow(RealDataObject9th)
RealDataObject9th <- DelDoublet(RealDataObject9th)

score9th <- score8th
score9th[match(rownames(RealDataObject9th@meta.data)[RealDataObject9th@meta.data[ ,ncol(RealDataObject9th@meta.data)] == "Doublet"],
               colnames(RealData[[1]]))] <- 1
saveRDS(score9th, paste0(DataID[i], "_9th_.rds"))

fg <- score9th[label==1]
bg <- score9th[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 18] <- pr$auc.integral
Results[i, 19] <- roc$auc


# tenth-time run ----
names(RealDataObject9th@meta.data)[RealDataObject9th@meta.data %>% ncol()] <- "Doubletstate9th"
RealDataObject10th <- subset(RealDataObject9th, Doubletstate9th == "Singlet")
RealDataObject10th <- SeuratWorkFlow(RealDataObject10th)
RealDataObject10th <- DelDoublet(RealDataObject10th)

score10th <- score9th
score10th[match(rownames(RealDataObject10th@meta.data)[RealDataObject10th@meta.data[ ,ncol(RealDataObject10th@meta.data)] == "Doublet"],
                colnames(RealData[[1]]))] <- 1
saveRDS(score10th, paste0(DataID[i], "_10th_.rds"))

fg <- score10th[label==1]
bg <- score10th[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 20] <- pr$auc.integral
Results[i, 21] <- roc$auc



# Multi-cycle ----
SplitSeurat <- SplitObject(RealDataObject, split.by = "seurat_clusters")
SplitSeurat <- lapply(SplitSeurat, function(i) {
  tryCatch(
    {
      SeuratWorkFlow(i)
      DelDoublet(i)
    },
    error = function(e) {
      NA
    }
  )
})
SplitSeurat <- keep(SplitSeurat, function(x) inherits(x, "Seurat"))
DoubletID <- lapply(SplitSeurat, function(x){
  names(x@meta.data)[x@meta.data %>% ncol()] <- "DoubletstateMulti"
  ID <- rownames(x@meta.data)[x$DoubletstateMulti == "Doublet"]
  return(ID)
})
DoubletID <- Reduce(c, DoubletID)

score3st <- score
score3st[match(DoubletID, colnames(RealData[[1]]))] <- 1
saveRDS(score3st, paste0(DataID[i], "_Multi_.rds"))

fg <- score3st[label==1]
bg <- score3st[label==0]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
Results[i, 22] <- pr$auc.integral
Results[i, 23] <- roc$auc

saveRDS(Results, "Results.rds")










