library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(SummarizedExperiment)

# ----------參數----------
Args <- commandArgs(TRUE)
id       <- Args[1]
ratioName <- Args[2]
pvalue_2_8 <- as.numeric(Args[3])
# ----------參數----------

data <- read.csv(file=paste("./file/", id, "/dep_output_result.csv", sep=""), header=TRUE, fileEncoding ="UTF-8", sep = ',')
## feature 1: numeric vector
geneList = data[,ratioName]

## feature 2: named vector
names(geneList) = as.character(data[, 'human_entrez'])

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)
print(head(geneList))

draw_plot2_8 <- function(geneList, pvalue_2_8) {
    print("--------------r_plot_upsetplot_with_splider_2_8() start--------------")
    kk2 <- gseKEGG( geneList = geneList,
                    organism = 'hsa',
                    minGSSize = 120,
                    pvalueCutoff = pvalue_2_8,
                    verbose = FALSE)
    if (length(kk2[,2]) < 1) {
        print("no term enriched under specific pvalueCutoff")
    }else {
        pic16_2 <- upsetplot(kk2)
        save_plot(paste("./file/", id, "/image/plot2_8.png", sep=""), pic16_2, base_height = 10, base_aspect_ratio = 1.5)
   }
}

draw_plot2_8(geneList, pvalue_2_8)