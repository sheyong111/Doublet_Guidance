# 2024/12/21

# Load packages
library(reticulate)
library(Matrix)
library(Seurat)
library(dplyr)
library(stringr)
library(PRROC)
library(pbapply)
# source('utility.R')
set.seed(2020)


# read python module
# Sys.setenv(RETICULATE_PYTHON = "~/soft/conda/miniconda3/envs/Doublet/bin/python")
reticulate::py_config() #检查配置
# use_condaenv("~/soft/conda/miniconda3/envs/Doublet/bin/python")
use_condaenv(condaenv = "Doublet", required = T,conda = "/home/sheyong/soft/conda/miniconda3/bin/conda")
py_available()

scr <- import('scrublet')  # 这个包加载的时候时常会报错， 重启一下R然后重新跑
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')



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
# list to save doublet scores
score.list <- list()

# i <- 4
for (i in 1:length(locs)) {   # i=14时发生报错
  
  # first run
  setwd("/home/sheyong/Rprogram/SingleCellDoublet/")
  dir.create(paste0("./Results/scrublet/", DataID[i]))
  RealData <- readRDS(locs[i])
  setwd(paste0("./Results/scrublet/", DataID[i]))
  label <- ifelse(RealData[[2]] == 'doublet', 1, 0)
  count <- RealData[[1]]
  DoubletRate <- ncol(count) / 500 * 0.004
  
  # Run Scrublet
  result <- scr$Scrublet(counts_matrix = t(count),
                         expected_doublet_rate = DoubletRate,
                         random_state = 10L)
  results <- result$scrub_doublets(min_counts = 1,
                                   min_cells = 1,
                                   min_gene_variability_pctl = 85,
                                   n_prin_comps = 30L)
  
  table(as.numeric(results[[2]]))
  score <- as.numeric(results[[2]])
  fg <- score[label==1]
  bg <- score[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  Results[i, 2] <- pr$auc.integral
  Results[i, 3] <- roc$auc
  saveRDS(score, paste0(DataID[i], "_1st.rds"))
  
  
  # Two cycle
  count <- count[ ,-which(score == 1)]
  DoubletRate <- ncol(count) / 500 * 0.004
  result <- scr$Scrublet(counts_matrix = t(count),
                         expected_doublet_rate = DoubletRate,
                         random_state = 10L)
  results <- result$scrub_doublets(min_counts = 1,
                                   min_cells = 1,
                                   min_gene_variability_pctl = 85,
                                   n_prin_comps = 30L)
  table(as.numeric(results[[2]]))
  
  score2 <- score
  score2[match(colnames(count)[which(as.numeric(results[[2]]) == 1)], colnames(RealData[[1]]))] <- 1
  table(score2)
  
  fg <- score2[label==1]
  bg <- score2[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  Results[i, 4] <- pr$auc.integral
  Results[i, 5] <- roc$auc
  saveRDS(score2, paste0(DataID[i], "_2nd.rds"))
  
  
  # Three cycle
  count <- count[ ,-which(score2 == 1)]
  DoubletRate <- ncol(count) / 500 * 0.004
  result <- scr$Scrublet(counts_matrix = t(count),
                         expected_doublet_rate = DoubletRate,
                         random_state = 10L)
  results <- result$scrub_doublets(min_counts = 1,
                                   min_cells = 1,
                                   min_gene_variability_pctl = 85,
                                   n_prin_comps = 30L)
  table(as.numeric(results[[2]]))
  
  score3 <- score2
  score3[match(colnames(count)[which(as.numeric(results[[2]]) == 1)], colnames(RealData[[1]]))] <- 1
  table(score)
  table(score2)
  table(score3)
  
  fg <- score3[label==1]
  bg <- score3[label==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  Results[i, 6] <- pr$auc.integral
  Results[i, 7] <- roc$auc
  saveRDS(score3, paste0(DataID[i], "_3rd.rds"))
  

  # # Multi cycle
  # RealDataObject <- CreateSeuratObject(counts = RealData[[1]],
  #                                      project = "RealDataObject")
  # RealDataObject <- SeuratWorkFlow(RealDataObject)
  # SplitSeurat <- SplitObject(RealDataObject, split.by = "seurat_clusters")
  # SplitSeurat <- lapply(SplitSeurat, function(x){x@assays$RNA@counts})
  # SplitSeurat_Results <- lapply(SplitSeurat, function(i) {
  #   tryCatch(
  #     {
  #       scr$Scrublet(counts_matrix = t(x),
  #                    expected_doublet_rate = ncol(x) / 500 * 0.004,
  #                    random_state = 10L)
  #       result$scrub_doublets(min_counts = 1,
  #                             min_cells = 1,
  #                             min_gene_variability_pctl = 85,
  #                             n_prin_comps = 30L)
  #     },
  #     error = function(e) {
  #       NA
  #     }
  #   )
  # })
  # 
  # score_multi <- score
  # for (x in 1:length(SplitSeurat_Results)) {
  #   if(is.na(SplitSeurat_Results[[x]])){
  #     next
  #   }
  #   names(as.numeric(SplitSeurat_Results[[x]][[2]])) <- SplitSeurat[[x]]
  #   score_multi[match(names(as.numeric(SplitSeurat_Results[[x]][[2]]))[which(as.numeric(SplitSeurat_Results[[x]][[2]]) == 1)], RealData[[1]])] <- 1
  # }
  # 
  # fg <- score_multi[label==1]
  # bg <- score_multi[label==0]
  # pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  # roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  # Results[i, 22] <- pr$auc.integral
  # Results[i, 23] <- roc$auc
  
}


saveRDS(Results, "/home/sheyong/Rprogram/SingleCellDoublet/Results/scrublet/Results.rds")


Results_scrublet_wider <- readRDS("/home/sheyong/Rprogram/SingleCellDoublet/Results/scrublet/Results.rds")

Results_scrublet <- reshape2::melt(Results_scrublet_wider[-c(4:5,14),1:7])



AUROC_dif <- subset(Results_scrublet, variable %in% c("AUROC1st", "AUROC2ed"))
mean_values <- aggregate(value ~ variable, data = AUROC_dif, FUN = mean)
P1 <- ggplot(AUROC_dif, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot()+
  stat_summary(fun.y = "mean", geom = "point", color = "red") +
  stat_summary(fun.y = "mean", geom = "line", aes(group = 1), color = "red", linetype = "dashed")+
  stat_compare_means(method = "wilcox.test",  paired = T,   # wilcox.test    t.test
                     comparisons = list(c("AUROC1st","AUROC2ed") ))+
  geom_text(data = mean_values, aes(label = round(value, 3)), color = "white", vjust = 1.5)+
  scale_fill_manual(values = c("#BCB27E","#B6090A"))+
  ylim(0,0.8)+
  commontheme1(0)+theme(axis.text.x = element_text(angle = 75))+
  ylab("AUROC")+xlab("")


AUPRC_dif <- subset(Results_scrublet, variable %in% c("AUPRC1st", "AUPRC2ed"))
mean_values <- aggregate(value ~ variable, data = AUPRC_dif, FUN = mean)
P2 <- ggplot(AUPRC_dif, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot()+
  stat_summary(fun.y = "mean", geom = "point", color = "red") +
  stat_summary(fun.y = "mean", geom = "line", aes(group = 1), color = "red", linetype = "dashed")+
  stat_compare_means(method = "wilcox.test",  paired = T,   # wilcox.test    t.test
                     comparisons = list(c("AUPRC1st","AUPRC2ed") ))+
  geom_text(data = mean_values, aes(label = round(value, 3)), color = "white", vjust = 1.5)+
  scale_fill_manual(values = c("#BCB27E","#B6090A"))+
  ylim(0,0.8)+
  commontheme1(0)+theme(axis.text.x = element_text(angle = 75))+
  ylab("AUPRC")+xlab("")


P1+P2
ggsave("Results/plot/scrublet_AUROC_AUPRC.pdf", width = 7, height = 6)





