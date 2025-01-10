# Parameters
args <- commandArgs(T)
options(stringsAsFactors=F)
# Parameter
# i: The number of datasets
i <- as.integer(args[1])

# Load packages
library(scds)
library(SingleCellExperiment)
library(Seurat)
library(magrittr)
library(stringr)
library(purrr)
library(PRROC)

# Set random seed
set.seed(19981031)

# Function ----
# RealDataObject: Matrix
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

# RealData: List from locs
# CombineFunctionResults: The results of CombineFunction

CombineFunction <- function(RealData){
  RealData <- RealData[[1]]
  example_sce <- SingleCellExperiment(assays = list(counts = RealData))
  # Run cxds, bcds and hybrid
  Results <- cxds_bcds_hybrid(example_sce)
  cxds_Results <- Results$cxds_score
  bcds_Results <- Results$bcds_score
  hybrid_Results <- colData(Results)$hybrid_score
  return(list(cxds_Results = cxds_Results, bcds_Results = bcds_Results, hybrid_Results = hybrid_Results))
}
Cal3Index <- function(CombineFunctionResults, RealData){
  Results <- data.frame(cxds = c(NA, NA), bcds = c(NA, NA), hybrid = c(NA, NA))
  score <- CombineFunctionResults$cxds_Results
  label <- ifelse(RealData[[2]] == 'doublet', 1, 0)
  fg <- score[label==1]
  bg <- score[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc
  Results[, 1] <- c(pr, roc)
  
  score <- CombineFunctionResults$bcds_Results
  fg <- score[label==1]
  bg <- score[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc
  Results[, 2] <- c(pr, roc)
  
  score <- CombineFunctionResults$hybrid_Results
  fg <- score[label==1]
  bg <- score[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc
  Results[, 3] <- c(pr, roc)
  rownames(Results) <- c("pr", "roc")
  return(Results)
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

cxds <- data.frame(DataID = DataID, AUPRC1st = NA, AUROC1st = NA,
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
bcds <- cxds
hybrid <- cxds
All_Results <- list(cxds = cxds, bcds = bcds, hybrid = hybrid)

# First-time run ----
# Construct SingleCellExperiment
setwd("/home/sheyong/Rprogram/SingleCellDoublet/")
dir.create(paste0("./Results/scds/scds_", DataID[i]))
RealData <- readRDS(locs[i])
setwd(paste0("./Results/scds/scds_", DataID[i]))
j <- 1 # Set the number of times
DoubletID1st <- CombineFunction(RealData)
names(DoubletID1st$bcds_Results) <- names(DoubletID1st$cxds_Results)
DoubletNum <- round((length(DoubletID1st$cxds_Results) / 500 * 0.004) * length(DoubletID1st$cxds_Results))
DoubletID1st_Results <- list(cxds_Results = rep(0, length(DoubletID1st$cxds_Results)),
                             bcds_Results = rep(0, length(DoubletID1st$cxds_Results)),
                             hybrid_Results = rep(0, length(DoubletID1st$cxds_Results)))
DoubletID1st_Results_ID <- lapply(DoubletID1st, function(x){order(x, decreasing = T)[1:DoubletNum]})
DoubletID1st_Results$cxds_Results[DoubletID1st_Results_ID$cxds_Results] <- 1
DoubletID1st_Results$bcds_Results[DoubletID1st_Results_ID$bcds_Results] <- 1
DoubletID1st_Results$hybrid_Results[DoubletID1st_Results_ID$hybrid_Results] <- 1
saveRDS(DoubletID1st_Results, paste0(DataID[i], "_DoubletID1st.rds"))
Index_Results <- Cal3Index(DoubletID1st_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]


# Second-time run ----
j <- 2 # Set the number of times
DoubletID2ed_Results <- DoubletID1st_Results

## cxds ----
DoubletNum <- round((length(DoubletID1st$cxds_Results) / 500 * 0.004) * length(DoubletID1st$cxds_Results))
cxds_DoubletID <- names(sort(DoubletID1st$cxds_Results, decreasing = T)[1:DoubletNum])
RealData2 <- list(RealData[[1]][ ,-match(cxds_DoubletID, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID, colnames(RealData[[1]]))])
DoubletID2ed <- CombineFunction(RealData2)
DoubletNum2ed <- round((length(DoubletID2ed$cxds_Results) / 500 * 0.004) * length(DoubletID2ed$cxds_Results))
DoubletID2ed_Results$cxds_Results[match(names(sort(DoubletID2ed$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID2 <- names(sort(DoubletID2ed$cxds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(cxds_DoubletID)

## bcds ----
DoubletNum <- round((length(DoubletID1st$bcds_Results) / 500 * 0.004) * length(DoubletID1st$bcds_Results))
bcds_DoubletID <- names(sort(DoubletID1st$bcds_Results, decreasing = T)[1:DoubletNum])
RealData2 <- list(RealData[[1]][ ,-match(bcds_DoubletID, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID, colnames(RealData[[1]]))])
DoubletID2ed <- CombineFunction(RealData2)
names(DoubletID2ed$bcds_Results) <- colnames(RealData2[[1]])
DoubletNum2ed <- round((length(DoubletID2ed$bcds_Results) / 500 * 0.004) * length(DoubletID2ed$bcds_Results))
DoubletID2ed_Results$bcds_Results[match(names(sort(DoubletID2ed$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID2 <- names(sort(DoubletID2ed$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID)

## hybrid ----
DoubletNum <- round((length(DoubletID1st$hybrid_Results) / 500 * 0.004) * length(DoubletID1st$hybrid_Results))
hybrid_DoubletID <- names(sort(DoubletID1st$hybrid_Results, decreasing = T)[1:DoubletNum])
RealData2 <- list(RealData[[1]][ ,-match(hybrid_DoubletID, colnames(RealData[[1]]))],
                  RealData[[2]][-match(hybrid_DoubletID, colnames(RealData[[1]]))])
DoubletID2ed <- CombineFunction(RealData2)
DoubletNum2ed <- round((length(DoubletID2ed$hybrid_Results) / 500 * 0.004) * length(DoubletID2ed$hybrid_Results))
DoubletID2ed_Results$hybrid_Results[match(names(sort(DoubletID2ed$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID2 <- names(sort(DoubletID2ed$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID)

# Combine 2 cycle results
saveRDS(DoubletID2ed_Results, paste0(DataID[i], "_DoubletID2ed.rds"))
Index_Results <- Cal3Index(DoubletID2ed_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]


# Third-time run ----
j <- 3 # Set the number of times
DoubletID3rd_Results <- DoubletID2ed_Results

## cxds ----
RealData3 <- list(RealData[[1]][ ,-match(cxds_DoubletID2, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID2, colnames(RealData[[1]]))])
DoubletID3rd <- CombineFunction(RealData2)
DoubletNum3rd <- round((length(DoubletID3rd$cxds_Results) / 500 * 0.004) * length(DoubletID3rd$cxds_Results))
DoubletID3rd_Results$cxds_Results[match(names(sort(DoubletID3rd$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID3 <- names(sort(DoubletID3rd$cxds_Results, decreasing = T)[1:DoubletNum3rd]) %>% union(cxds_DoubletID2)

## bcds ----
RealData2 <- list(RealData[[1]][ ,-match(bcds_DoubletID2, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID2, colnames(RealData[[1]]))])
DoubletID3rd <- CombineFunction(RealData2)
names(DoubletID3rd$bcds_Results) <- colnames(RealData2[[1]])
DoubletNum2ed <- round((length(DoubletID3rd$bcds_Results) / 500 * 0.004) * length(DoubletID3rd$bcds_Results))
DoubletID3rd_Results$bcds_Results[match(names(sort(DoubletID3rd$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID3 <- names(sort(DoubletID3rd$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID2)

## hybrid ----
RealData2$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID2, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID2, colnames(RealData[[1]]))])
DoubletID3rd <- CombineFunction(RealData2)
DoubletNum2ed <- round((length(DoubletID3rd$hybrid_Results) / 500 * 0.004) * length(DoubletID3rd$hybrid_Results))
DoubletID3rd_Results$hybrid_Results[match(names(sort(DoubletID3rd$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID3 <- names(sort(DoubletID3rd$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID2)

# Combine 2 cycle results
saveRDS(DoubletID3rd_Results, paste0(DataID[i], "_DoubletID3rd.rds"))
Index_Results <- Cal3Index(DoubletID3rd_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]



# Fouth-time run ----
j <- 4 # Set the number of times
DoubletID4th_Results <- DoubletID3rd_Results

## cxds ----
RealData3 <- list(RealData[[1]][ ,-match(cxds_DoubletID3, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID3, colnames(RealData[[1]]))])
DoubletID4th <- CombineFunction(RealData3)
DoubletNum3rd <- round((length(DoubletID4th$cxds_Results) / 500 * 0.004) * length(DoubletID4th$cxds_Results))
DoubletID4th_Results$cxds_Results[match(names(sort(DoubletID4th$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID4 <- names(sort(DoubletID4th$cxds_Results, decreasing = T)[1:DoubletNum3rd]) %>% union(cxds_DoubletID3)

## bcds ----
RealData3 <- list(RealData[[1]][ ,-match(bcds_DoubletID3, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID3, colnames(RealData[[1]]))])
DoubletID4th <- CombineFunction(RealData3)
names(DoubletID4th$bcds_Results) <- colnames(RealData3[[1]])
DoubletNum2ed <- round((length(DoubletID4th$bcds_Results) / 500 * 0.004) * length(DoubletID4th$bcds_Results))
DoubletID4th_Results$bcds_Results[match(names(sort(DoubletID4th$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID4 <- names(sort(DoubletID4th$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID3)

## hybrid ----
RealData3$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID3, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID3, colnames(RealData[[1]]))])
DoubletID4th <- CombineFunction(RealData3)
DoubletNum2ed <- round((length(DoubletID4th$hybrid_Results) / 500 * 0.004) * length(DoubletID4th$hybrid_Results))
DoubletID4th_Results$hybrid_Results[match(names(sort(DoubletID4th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID4 <- names(sort(DoubletID4th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID3)

# Combine 2 cycle results
saveRDS(DoubletID4th_Results, paste0(DataID[i], "_DoubletID4th.rds"))
Index_Results <- Cal3Index(DoubletID4th_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]



# Fifth-time run ----
j <- 5 # Set the number of times
DoubletID5th_Results <- DoubletID4th_Results

## cxds ----
RealData4 <- list(RealData[[1]][ ,-match(cxds_DoubletID4, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID4, colnames(RealData[[1]]))])
DoubletID5th <- CombineFunction(RealData4)
DoubletNum2ed <- round((length(DoubletID5th$cxds_Results) / 500 * 0.004) * length(DoubletID5th$cxds_Results))
DoubletID5th_Results$cxds_Results[match(names(sort(DoubletID5th$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID5 <- names(sort(DoubletID5th$cxds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(cxds_DoubletID4)

## bcds ----
RealData4 <- list(RealData[[1]][ ,-match(bcds_DoubletID4, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID4, colnames(RealData[[1]]))])
DoubletID5th <- CombineFunction(RealData4)
names(DoubletID5th$bcds_Results) <- colnames(RealData4[[1]])
DoubletNum2ed <- round((length(DoubletID5th$bcds_Results) / 500 * 0.004) * length(DoubletID5th$bcds_Results))
DoubletID5th_Results$bcds_Results[match(names(sort(DoubletID5th$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID5 <- names(sort(DoubletID5th$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID4)

## hybrid ----
RealData4$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID4, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID4, colnames(RealData[[1]]))])
DoubletID5th <- CombineFunction(RealData4)
DoubletNum2ed <- round((length(DoubletID5th$hybrid_Results) / 500 * 0.004) * length(DoubletID5th$hybrid_Results))
DoubletID5th_Results$hybrid_Results[match(names(sort(DoubletID5th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID5 <- names(sort(DoubletID5th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID4)

# Combine 2 cycle results
saveRDS(DoubletID5th_Results, paste0(DataID[i], "_DoubletID5th.rds"))
Index_Results <- Cal3Index(DoubletID5th_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]



# Sixth-time run ----
j <- 6 # Set the number of times
DoubletID6th_Results <- DoubletID5th_Results

## cxds ----
RealData5 <- list(RealData[[1]][ ,-match(cxds_DoubletID5, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID5, colnames(RealData[[1]]))])
DoubletID6th <- CombineFunction(RealData5)
DoubletNum2ed <- round((length(DoubletID6th$cxds_Results) / 500 * 0.004) * length(DoubletID6th$cxds_Results))
DoubletID6th_Results$cxds_Results[match(names(sort(DoubletID6th$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID6 <- names(sort(DoubletID6th$cxds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(cxds_DoubletID5)

## bcds ----
RealData5 <- list(RealData[[1]][ ,-match(bcds_DoubletID5, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID5, colnames(RealData[[1]]))])
DoubletID6th <- CombineFunction(RealData5)
names(DoubletID6th$bcds_Results) <- colnames(RealData5[[1]])
DoubletNum2ed <- round((length(DoubletID6th$bcds_Results) / 500 * 0.004) * length(DoubletID6th$bcds_Results))
DoubletID6th_Results$bcds_Results[match(names(sort(DoubletID6th$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID6 <- names(sort(DoubletID6th$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID5)

## hybrid ----
RealData5$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID5, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID5, colnames(RealData[[1]]))])
DoubletID6th <- CombineFunction(RealData5)
DoubletNum2ed <- round((length(DoubletID6th$hybrid_Results) / 500 * 0.004) * length(DoubletID6th$hybrid_Results))
DoubletID6th_Results$hybrid_Results[match(names(sort(DoubletID6th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID6 <- names(sort(DoubletID6th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID5)

# Combine 2 cycle results
saveRDS(DoubletID6th_Results, paste0(DataID[i], "_DoubletID6th.rds"))
Index_Results <- Cal3Index(DoubletID6th_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]



# seventh-time run ----
j <- 7 # Set the number of times
DoubletID7th_Results <- DoubletID6th_Results
RealData6 <- list(RealData[[1]][ ,-match(cxds_DoubletID6, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID6, colnames(RealData[[1]]))])
DoubletID7th <- CombineFunction(RealData6)
DoubletNum2ed <- round((length(DoubletID7th$cxds_Results) / 500 * 0.004) * length(DoubletID7th$cxds_Results))
DoubletID7th_Results$cxds_Results[match(names(sort(DoubletID7th$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID7 <- names(sort(DoubletID7th$cxds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(cxds_DoubletID6)

## bcds ----
RealData6 <- list(RealData[[1]][ ,-match(bcds_DoubletID6, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID6, colnames(RealData[[1]]))])
DoubletID7th <- CombineFunction(RealData6)
names(DoubletID7th$bcds_Results) <- colnames(RealData6[[1]])
DoubletNum2ed <- round((length(DoubletID7th$bcds_Results) / 500 * 0.004) * length(DoubletID7th$bcds_Results))
DoubletID7th_Results$bcds_Results[match(names(sort(DoubletID7th$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID7 <- names(sort(DoubletID7th$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID6)

## hybrid ----
RealData6$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID6, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID6, colnames(RealData[[1]]))])
DoubletID7th <- CombineFunction(RealData6)
DoubletNum2ed <- round((length(DoubletID7th$hybrid_Results) / 500 * 0.004) * length(DoubletID7th$hybrid_Results))
DoubletID7th_Results$hybrid_Results[match(names(sort(DoubletID7th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID7 <- names(sort(DoubletID7th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID6)

# Combine 2 cycle results
saveRDS(DoubletID7th_Results, paste0(DataID[i], "_DoubletID7th.rds"))
Index_Results <- Cal3Index(DoubletID7th_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]



# eightth-time run ----
j <- 8 # Set the number of times
DoubletID8th_Results <- DoubletID7th_Results

## cxds ----
RealData7 <- list(RealData[[1]][ ,-match(cxds_DoubletID7, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID7, colnames(RealData[[1]]))])
DoubletID8th <- CombineFunction(RealData7)
DoubletNum2ed <- round((length(DoubletID8th$cxds_Results) / 500 * 0.004) * length(DoubletID8th$cxds_Results))
DoubletID8th_Results$cxds_Results[match(names(sort(DoubletID8th$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID8 <- names(sort(DoubletID8th$cxds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(cxds_DoubletID7)

## bcds ----
RealData7 <- list(RealData[[1]][ ,-match(bcds_DoubletID7, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID7, colnames(RealData[[1]]))])
DoubletID8th <- CombineFunction(RealData7)
names(DoubletID8th$bcds_Results) <- colnames(RealData7[[1]])
DoubletNum2ed <- round((length(DoubletID8th$bcds_Results) / 500 * 0.004) * length(DoubletID8th$bcds_Results))
DoubletID8th_Results$bcds_Results[match(names(sort(DoubletID8th$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID8 <- names(sort(DoubletID8th$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID7)

## hybrid ----
RealData7$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID7, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID7, colnames(RealData[[1]]))])
DoubletID8th <- CombineFunction(RealData7)
DoubletNum2ed <- round((length(DoubletID8th$hybrid_Results) / 500 * 0.004) * length(DoubletID8th$hybrid_Results))
DoubletID8th_Results$hybrid_Results[match(names(sort(DoubletID8th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID8 <- names(sort(DoubletID8th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID7)

# Combine 2 cycle results
saveRDS(DoubletID8th_Results, paste0(DataID[i], "_DoubletID8th.rds"))
Index_Results <- Cal3Index(DoubletID8th_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]



# nineth-time run ----
j <- 9 # Set the number of times
DoubletID9th_Results <- DoubletID8th_Results

## cxds ----
RealData8 <- list(RealData[[1]][ ,-match(cxds_DoubletID8, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID8, colnames(RealData[[1]]))])
DoubletID9th <- CombineFunction(RealData8)
DoubletNum2ed <- round((length(DoubletID9th$cxds_Results) / 500 * 0.004) * length(DoubletID9th$cxds_Results))
DoubletID9th_Results$cxds_Results[match(names(sort(DoubletID9th$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID9 <- names(sort(DoubletID9th$cxds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(cxds_DoubletID8)

## bcds ----
RealData8 <- list(RealData[[1]][ ,-match(bcds_DoubletID8, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID8, colnames(RealData[[1]]))])
DoubletID9th <- CombineFunction(RealData8)
names(DoubletID9th$bcds_Results) <- colnames(RealData8[[1]])
DoubletNum2ed <- round((length(DoubletID9th$bcds_Results) / 500 * 0.004) * length(DoubletID9th$bcds_Results))
DoubletID9th_Results$bcds_Results[match(names(sort(DoubletID9th$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID9 <- names(sort(DoubletID9th$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID8)

## hybrid ----
RealData8$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID8, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID8, colnames(RealData[[1]]))])
DoubletID9th <- CombineFunction(RealData8)
DoubletNum2ed <- round((length(DoubletID9th$hybrid_Results) / 500 * 0.004) * length(DoubletID9th$hybrid_Results))
DoubletID9th_Results$hybrid_Results[match(names(sort(DoubletID9th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID9 <- names(sort(DoubletID9th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID8)

# Combine 2 cycle results
saveRDS(DoubletID9th_Results, paste0(DataID[i], "_DoubletID9th.rds"))
Index_Results <- Cal3Index(DoubletID9th_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]



# tenth-time run ----
j <- 10 # Set the number of times
DoubletID10th_Results <- DoubletID9th_Results

## cxds ----
RealData9 <- list(RealData[[1]][ ,-match(cxds_DoubletID8, colnames(RealData[[1]]))],
                  RealData[[2]][-match(cxds_DoubletID8, colnames(RealData[[1]]))])
DoubletID10th <- CombineFunction(RealData9)
DoubletNum2ed <- round((length(DoubletID10th$cxds_Results) / 500 * 0.004) * length(DoubletID10th$cxds_Results))
DoubletID10th_Results$cxds_Results[match(names(sort(DoubletID10th$cxds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
cxds_DoubletID9 <- names(sort(DoubletID10th$cxds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(cxds_DoubletID8)

## bcds ----
RealData9 <- list(RealData[[1]][ ,-match(bcds_DoubletID8, colnames(RealData[[1]]))],
                  RealData[[2]][-match(bcds_DoubletID8, colnames(RealData[[1]]))])
DoubletID10th <- CombineFunction(RealData9)
names(DoubletID10th$bcds_Results) <- colnames(RealData9[[1]])
DoubletNum2ed <- round((length(DoubletID10th$bcds_Results) / 500 * 0.004) * length(DoubletID10th$bcds_Results))
DoubletID10th_Results$bcds_Results[match(names(sort(DoubletID10th$bcds_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
bcds_DoubletID9 <- names(sort(DoubletID10th$bcds_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(bcds_DoubletID8)

## hybrid ----
RealData9$bcds <- list(RealData[[1]][ ,-match(hybrid_DoubletID8, colnames(RealData[[1]]))],
                       RealData[[2]][-match(hybrid_DoubletID8, colnames(RealData[[1]]))])
DoubletID10th <- CombineFunction(RealData9)
DoubletNum2ed <- round((length(DoubletID10th$hybrid_Results) / 500 * 0.004) * length(DoubletID10th$hybrid_Results))
DoubletID10th_Results$hybrid_Results[match(names(sort(DoubletID10th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]), colnames(RealData[[1]]))] <- 1
hybrid_DoubletID9 <- names(sort(DoubletID10th$hybrid_Results, decreasing = T)[1:DoubletNum2ed]) %>% union(hybrid_DoubletID8)

# Combine 2 cycle results
saveRDS(DoubletID10th_Results, paste0(DataID[i], "_DoubletID10th.rds"))
Index_Results <- Cal3Index(DoubletID10th_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]




# Multi-cycle ----
j <- 11 # Set the number of times
RealDataObject <- CreateSeuratObject(counts = RealData[[1]],
                                     project = "RealDataObject")
RealDataObject <- SeuratWorkFlow(RealDataObject)
SplitSeurat <- SplitObject(RealDataObject, split.by = "seurat_clusters")
SplitSeurat <- lapply(SplitSeurat, function(x){x@assays$RNA@counts})
CombineFunctionMulti <- function(RealData){
  example_sce <- SingleCellExperiment(assays = list(counts = RealData))
  # Run cxds, bcds and hybrid
  Results <- cxds_bcds_hybrid(example_sce)
  cxds_Results <- Results$cxds_score
  bcds_Results <- Results$bcds_score
  hybrid_Results <- colData(Results)$hybrid_score
  return(list(cxds_Results = cxds_Results, bcds_Results = bcds_Results, hybrid_Results = hybrid_Results))
}
DoubletIDMulti <- lapply(SplitSeurat, CombineFunctionMulti)
merged_list <- unlist(DoubletIDMulti, recursive = FALSE)
Barcode <- lapply(SplitSeurat, function(x){colnames(x)})

# Add barcode to bcds_Results
for (x in 1:length(grep("bcds_Results", names(merged_list)))) {
  names(merged_list[[grep("bcds_Results", names(merged_list))[x]]]) <- Barcode[[x]]
}

# Calculate Doublet numbers in each cluster
DoubletID <- unlist(lapply(merged_list, function(x){round((length(x) / 500 * 0.004) * length(x))}))

for (y in 1:length(DoubletID)) {
  merged_list[[y]] <- sort(merged_list[[y]], decreasing = T)[1:DoubletID[y]]
}
DoubletIDMulti <- list(cxds_Results = Reduce(c, merged_list[grep("cxds_Results", names(merged_list))]),
                       bcds_Results = Reduce(c, merged_list[grep("bcds_Results", names(merged_list))]),
                       hybrid_Results = Reduce(c, merged_list[grep("hybrid_Results", names(merged_list))]))
DoubletIDMulti_Results <- DoubletID1st_Results
for (z in 1:length(DoubletIDMulti)) {
  DoubletIDMulti_Results[[1]][match(names(DoubletIDMulti[[1]]), colnames(RealData[[1]]))] <- 1
  DoubletIDMulti_Results[[2]][match(names(DoubletIDMulti[[2]]), colnames(RealData[[1]]))] <- 1
  DoubletIDMulti_Results[[3]][match(names(DoubletIDMulti[[3]]), colnames(RealData[[1]]))] <- 1
}

# Combine Multi cycle results
saveRDS(DoubletIDMulti_Results, paste0(DataID[i], "_DoubletIDMulti.rds"))
Index_Results <- Cal3Index(DoubletIDMulti_Results, RealData)
All_Results$cxds[i, c(2*j ,2*j+1)] <- Index_Results[ ,1]
All_Results$bcds[i, c(2*j ,2*j+1)] <- Index_Results[ ,2]
All_Results$hybrid[i, c(2*j ,2*j+1)] <- Index_Results[ ,3]

saveRDS(All_Results, paste0(DataID[i], "_All_Results.rds"))

