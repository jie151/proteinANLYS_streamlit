library(DEP)
library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(NormalyzerDE)
library(SummarizedExperiment)
source("normalize_function.r")

print("--------------dep3_1_9.r start--------------")
# ----------參數----------
Args <- commandArgs(TRUE)
filename <- Args[1]
id       <- Args[2]
colname_geneNames  <- Args[3]
colname_proteinIDs <- Args[4]
normalizeOption    <- Args[5]
filter_option <- Args[6]
nThr          <- as.numeric(Args[7])
filter_min    <- as.numeric(Args[8])
control_r <- Args[9]
alpha_r <- as.numeric(Args[10])
lfc_r <- as.numeric(Args[11])
protein_type_r <- Args[12]


proteinData <- readLines(paste("./file/", id, "/proteinData.txt", sep=""))
print("--------------------")
print(protein_type_r)
print(proteinData)
print("--------------------")

# ----------參數----------
exp_design_filename <- paste("./file/", id, "/exp_design.csv", sep="")
pic_filename <- paste("./file/", id, "/image/", sep="")

data <- read.csv(file=filename, header=TRUE, fileEncoding ="UTF-8", sep = '\t')
cat("file's row * column =", dim(data), "\n")
if ( "Reverse" %in% names(data) ){ #如果有的話，要做篩選
    data <- filter(data, Reverse != "+")
}
if ( "Potential.contaminant" %in% names(data) ){
    data <- filter(data, Potential.contaminant != "+")
}


data_unique <- make_unique(data, colname_geneNames, colname_proteinIDs, delim = ";")
columns <- grep("Reporter.intensity.corrected.", colnames(data_unique)) # get LFQ column numbers
exp_design <- read.csv(exp_design_filename, header=TRUE ,fileEncoding ="UTF-8")
exp_design$label = gsub(" ", ".", exp_design$label)
exp_design$condition <- make.names(exp_design$condition)
exp_design_condition <- unique(exp_design$condition)
data_se <- make_se(data_unique, columns, exp_design)

data_filt <- filter_proteins(data_se, type = filter_option, thr = nThr, min = filter_min)
data_norm <- normalize_proteiNorm(data_filt, normalizeOption)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

data_diff <- test_diff(data_imp, type = "control", control = control_r)
dep <- add_rejections(data_diff, alpha = alpha_r, lfc = log2(lfc_r)) #alpha、lfc調整

pic1 <- plot_single(dep, proteins = unlist( proteinData) , type = protein_type_r)
save_plot(paste("./file/", id, "/image/plot1_9.png", sep=""), pic1)