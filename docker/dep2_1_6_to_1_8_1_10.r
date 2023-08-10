library(DEP)
library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(NormalyzerDE)
library(SummarizedExperiment)
library(venneuler)
library(eulerr)
source("normalize_function.r")

print("--------------dep2_1_6.r start--------------")
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
contrast_r <- Args[10]
alpha_r <- as.numeric(Args[11])
lfc_r <- as.numeric(Args[12])

cat("!!!!control: ", control_r, "\n")

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
write.csv(data_results, paste("./file/", id, "/dep_output.csv", sep=""), row.names = FALSE, quote=F)


# ----------畫圖----------
draw_plot1_6_to_1_8_1_10 <- function(data_se, data_filt, data_norm, data_imp, dep) {
    # plot1_6
    tryCatch(
        {
            pic1 <- plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)
            save_plot(paste("./file/", id, "/image/plot1_6_1.png", sep=""), pic1)
        },
        error = function(e) {
            message("Error!!!!! plot1_6_1.png ")
            print(e)
        }
    )
    tryCatch(
        {
            png(file=paste("./file/", id, "/image/plot1_6_2.png", sep=""))
            pic1 <- plot_cor(dep, significant = FALSE, lower = 0.9, upper = 1, pal = "Reds")
            dev.off()
        },
        error = function(e) {
            message("Error!!!!! plot1_6_2.png ")
            print(e)
        }
    )
    # plot1_7
    tryCatch(
        {
            png(file=paste("./file/", id, "/image/plot1_7_1.png", sep=""))
            pic1 <- plot_heatmap(dep, type = "centered", kmeans = TRUE,
                        k = 6, col_limit = 4, show_row_names = FALSE,
                        indicate = c("condition", "replicate"))
            dev.off()
        },
        error = function(e) {
            message("Error!!!!! plot1_7_1.png ")
            print(e)
        }
    )
    tryCatch(
        {
            png(file=paste("./file/", id, "/image/plot1_7_2.png", sep=""))
            pic1 <- plot_heatmap(dep, type = "contrast", kmeans = TRUE,
                        k = 6, col_limit = 10, show_row_names = FALSE)
            dev.off()
        },
        error = function(e) {
            message("Error!!!!! plot1_7_2.png ")
            print(e)
        }
    )
    # plot1_8
    tryCatch(
        {
            pic <- plot_volcano(dep, contrast = contrastSample, label_size = 3, add_names = TRUE,adjusted = FALSE, plot = TRUE)
            save_plot(paste("./file/", id, "/image/plot1_8.png", sep=""), pic)
        },
        error = function(e) {
            message("Error!!!!! plot1_7_2.png ")
            print(e)
        }
    )
    # plot1_10
    tryCatch(
        {
            png(file=paste("./file/", id, "/image/plot1_10.png", sep=""))
            pic1 <- plot_cond(dep)
            dev.off()

            # merged_list <- paste(cond_list_colnames, collapse = ",;")
            # merged_list2 <- paste(cond_list_proteins, collapse = ",;")

            # # 將新的list轉換為一個字串，並加入換行符號"\n"
            # output_string <- paste(merged_list, merged_list2, sep = "\n")

            # # 將該字串寫入txt檔案
            # file_path <- paste("./file/", id, "/venn_data.txt", sep="")
            # writeLines(output_string, con = file_path)
        },
        error = function(e) {
            message("Error!!!!! plot1_10.png ")
            print(e)
        }
    )

}

cond_list <- plot_cond(dep, plot=FALSE)
cond_list_colnames <- t(cond_list$counts['conditions'])
cond_list_proteins <- t(cond_list$counts['proteins'])
cond_list_length <- length(cond_list_colnames)

for (i in c(1: cond_list_length)){
    cond_list_colnames[i] = gsub(" ", "&", cond_list_colnames[i] )
}
colnames <- cond_list_colnames
data <- as.numeric(cond_list_proteins)
venn_data <- setNames(data, colnames)
print(venn_data)

vd <- euler(venn_data)

png(file=paste("./file/", id, "/image/plot1_10_2.png", sep=""))
plot(vd,
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    alpha = 0.45,
    labels = list(col = "red", font = 8),
    edges = list(col = "black", lex = 2),
    quantities = TRUE,
    legend = TRUE)

dev.off()

draw_plot1_6_to_1_8_1_10(data_se, data_filt, data_norm, data_imp, dep)