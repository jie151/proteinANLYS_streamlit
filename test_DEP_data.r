library(DEP)
library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(httr)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(NormalyzerDE)
library(SummarizedExperiment)
library(biomaRt)

# ---------範例資料----------
# data <- UbiLength
# data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
# data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
# columns <- grep("LFQ.", colnames(data_unique))
# exp_design <- UbiLength_ExpDesign

# ---------檔案----------
data <- read.csv(file="1.proteinGroups_Recovery (1).txt", header=TRUE, fileEncoding ="UTF-8", sep = '\t')
exp_design <- read.csv('1.exp_design.csv',header=TRUE ,fileEncoding ="UTF-8")

cat("file's row * column =", dim(data), "\n")
colname_proteinIDs <- "Protein.IDs"
colname_geneNames <- "Gene.names"
cat('* Are there any duplicated gene names? ', data[, colname_geneNames] %>% duplicated() %>% any(), "\n")

data_unique <- make_unique(data, colname_geneNames, colname_proteinIDs, delim = ";")
columns <- grep("Reporter.intensity.corrected.", colnames(data_unique)) # get LFQ column numbers
exp_design$label = gsub(" ", ".", exp_design$label)

maxReplicate <- max(exp_design$replicate) # max replicate
exp_design_condition <- unique(exp_design$condition)

# ----------參數----------
nThr <- 1
normalizeOption <- "Log2"
control_r <- "T0"
contrast_r <- "T1"
alpha_r <- 0.05
lfc_r <- 2
contrastSample <- paste(contrast_r, "_vs_", control_r, sep = "")
# ----------參數----------
data_se <- make_se(data_unique, LFQ_columns, exp_design)
data_filt <- filter_missval(data_se, thr = nThr) #讓使用者選0~4(重複)
data_norm <- normalize_proteiNorm(data_filt, normalizeOption)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_diff <- test_diff(data_imp, type = "control", control = control_r)
dep <- add_rejections(data_diff, alpha = alpha_r, lfc = log2(lfc_r)) #alpha、lfc調整

plot_volcano(dep, contrast = contrastSample, label_size = 3, add_names = TRUE,adjusted = FALSE, plot = TRUE)
