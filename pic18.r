library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(DOSE)
library(enrichplot)


Args <- commandArgs(TRUE)
ratioName <- Args[1]


data <- read.csv(file="./file/dep_output_result.csv", header=TRUE, fileEncoding ="UTF-8", sep = ',')
## feature 1: numeric vector
geneList = data[,ratioName] #-log(pValue)

## feature 2: named vector
names(geneList) = as.character(data[, 'human_entrez'])

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)
print(head(geneList))

de <- names(geneList)[abs(log(geneList)) > 1]
edo <- enrichDGN(de)
edo2 <- gseDO(geneList, pvalueCutoff=1)


pic18_1 <- gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
pic18_2 <- gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])

pic18_3 <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])

pic18_4 <- gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])
pic18_5 <- gseaplot2(edo2, geneSetID = 1:3)
pic18_6 <- gseaplot2(edo2, geneSetID = 1:3, pvalue_table = TRUE,
    color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

pic18_7 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1) + xlab("Rank")
pic18_8 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1:2)# + xlab("Rank")

pic18_9 <- gsearank(edo2, 1, title = edo2[1, "Description"])

pic18_10 <- lapply(1:3, function(i){
    anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
    gsearank(edo2, i, edo2[i, 2]) + xlab("Rank")+
    annotate("text", 5000, edo2[i, "enrichmentScore"]* .75, label =lab, hjust=0, vjust=0)})


save_plot("./file/image/plot2_10_1.png", pic18_1, base_height = 10, base_aspect_ratio = 1)
save_plot("./file/image/plot2_10_2.png", pic18_2, base_height = 10, base_aspect_ratio = 1)



png("./file/image/plot2_10_3.png")
    pic18_3
dev.off()

png("./file/image/plot2_10_4.png")
    pic18_4
dev.off()

png("./file/image/plot2_10_5.png")
    pic18_5
dev.off()

png("./file/image/plot2_10_6.png")
    pic18_6
dev.off()

png("./file/image/plot2_10_7.png")
    pic18_7
dev.off()

png("./file/image/plot2_10_8.png")
    pic18_8
dev.off()

png("./file/image/plot2_10_9.png")
    pic18_9
dev.off()

png("./file/image/plot2_10_10.png")
    pic18_10
dev.off()