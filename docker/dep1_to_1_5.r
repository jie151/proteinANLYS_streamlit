library(DEP)
library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(NormalyzerDE)
library(SummarizedExperiment)
source("normalize_function.r")

print("--------------pic1_1_to_1_5.r start--------------")
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
write.csv(data_unique, paste("./file/", id, "/data_unique.csv", sep=""), row.names = FALSE, quote=F)

columns <- grep("Reporter.intensity.corrected.", colnames(data_unique)) # get LFQ column numbers
exp_design <- read.csv(exp_design_filename, header=TRUE ,fileEncoding ="UTF-8")
exp_design$label = gsub(" ", ".", exp_design$label)
exp_design$condition <- make.names(exp_design$condition)
exp_design_condition <- unique(exp_design$condition)
data_se <- make_se(data_unique, columns, exp_design)

data_filt <- filter_proteins(data_se, type = filter_option, thr = nThr, min = filter_min)
data_norm <- normalize_proteiNorm(data_filt, normalizeOption)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Let's have a look at the SummarizedExperiment object
cat("* the SummarizedExperiment object:  \n")
print(data_se)

# ----------畫圖----------
draw_plot1_1_to_1_5 <- function(data_se, data_filt, data_norm, data_imp) {
    # plot1_1
    tryCatch(
        {
            save_plot(paste(pic_filename, "plot1_1.png", sep=""), plot =  plot_frequency(data_se), , base_height = 4, base_width = 4.5)
        },
        error = function(e) {
            message("Error!!!!! plot1_1.png ")
            print(e)
        }
    )
    # plot1_2
    tryCatch(
        {
            save_plot(paste(pic_filename, "plot1_2.png", sep=""), plot = plot_numbers(data_filt), base_height = 4, base_width = 4.5)
        },
        error = function(e) {
            message("Error!!!!! plot1_2.png ")
            print(e)
        }
    )
    # plot1_3
    tryCatch(
        {
            save_plot(paste(pic_filename, "plot1_3.png", sep=""), plot = plot_coverage(data_filt), base_height = 4, base_width = 4)
        },
        error = function(e) {
            message("Error!!!!! plot1_3.png ")
            print(e)
        }
    )
    # plot1_4
    tryCatch(
        {
            pic1 <- plot_normalization(data_filt, data_norm)
            save_plot(paste(pic_filename, "plot1_4.png", sep=""), plot =  pic1, base_height = 4, base_width = 4)
        },
        error = function(e) {
            message("Error!!!!! plot1_4.png ")
            print(e)
        }
    )
    # plot1_5_1
    tryCatch(
        {
            png(file = paste(pic_filename,"plot1_5_1.png", sep=""))
            pic1 <- plot_missval(data_filt)
            dev.off()
        },
        error = function(e) {
            message("Error!!!!! plot1_5_1.png ")
            print(e)
        }
    )
    # plot1_5_2
    tryCatch(
        {
            png(file = paste(pic_filename,"plot1_5_2.png", sep=""))
            pic1 <- plot_detect(data_filt)
            dev.off()

            pic1 <- plot_imputation(data_norm, data_imp)
            save_plot( paste(pic_filename, "plot1_5_3.png", sep=""), plot =  pic1, base_width = 4, base_height = 4)
        },
        error = function(e) {
            message("Error!!!!! plot1_5.png ")
            print(e)
        }
    )
    # plot1_5_3
    tryCatch(
        {
            pic1 <- plot_imputation(data_norm, data_imp)
            save_plot( paste(pic_filename, "plot1_5_3.png", sep=""), plot =  pic1, base_width = 4, base_height = 4)
        },
        error = function(e) {
            message("Error!!!!! plot1_5.png ")
            print(e)
        }
    )
}

draw_plot1_1_to_1_5(data_se, data_filt, data_norm, data_imp)