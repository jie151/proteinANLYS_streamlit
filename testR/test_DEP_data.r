library(DEP)
library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(NormalyzerDE)
library(SummarizedExperiment)
# source("normalize_function.r")
source("C:/Users/jie/Desktop/pyEnv/web_proteinAnls/code/docker/normalize_function.r")
print("--------------pic1_1_to_1_5.r start--------------")
# ----------參數----------
Args <- commandArgs(TRUE)
filename <- Args[1]
id       <- Args[2]
species  <- Args[3]
colname_geneNames  <- Args[4]
colname_proteinIDs <- Args[5]
normalizeOption    <- Args[6]
filter_option <- Args[7]
nThr          <- as.numeric(Args[8])
filter_min    <- as.numeric(Args[9])
control_r <- Args[10]
contrast_r <- Args[11]
alpha_r <- as.numeric(Args[12])
lfc_r <- as.numeric(Args[13])
# ----------參數----------

# ----------參數----------
filename <- "./file/test/proteinGroups_HsinYuan_Rat.txt"
id       <- "test"
species  <- "rat"
colname_geneNames  <- "Gene.names"
colname_proteinIDs <- "Protein.IDs"
normalizeOption    <- "Log2"
nThr <- 1
filter_option <- "condition"
filter_min    <- 0.66
control_r <- "A"
contrast_r <- "B"
alpha_r <- 0.05
lfc_r <- 2
# ----------參數----------

contrastSample <- paste(contrast_r, "_vs_", control_r, sep = "")
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

# Are there any duplicated gene names?
cat('* Are there any duplicated gene names? ', data[, colname_geneNames] %>% duplicated() %>% any(), "\n")

if ( data[, colname_geneNames] %>% duplicated() %>% any() ){
    # Make a table of duplicated gene names
    write.table(data %>% group_by_(.dots = colname_geneNames) %>% summarize(frequency = n()) %>%
        arrange(desc(frequency)) %>% filter(frequency > 1), file = "./file/my_data1.txt", row.names =FALSE)

    cat("a table of duplicated gene names: (table ", dim(table), "\n")
    table <- read.csv(file="./file/my_data1.txt", header=TRUE, fileEncoding ="UTF-8", sep = ' ')
    print(head(table, 7))
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
data_results <- get_results(dep)
write.csv(data_results,"./file/dep_output.csv", row.names = FALSE, quote=F)
