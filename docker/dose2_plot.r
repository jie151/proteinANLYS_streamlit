library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(SummarizedExperiment)
library(ReactomePA)
library(AnnotationHub)
library(msigdbr)

# ----------參數----------
Args <- commandArgs(TRUE)
id       <- Args[1]
ratioName <- Args[2]
de_up_down <- Args[3]
range <- as.numeric(Args[4])
enrichment_analysis_methods <- Args[5]
universal_enrichment_category <- Args[6]
universal_enrichment_subcategory <- Args[7]

# ----------參數----------

data <- read.csv(file=paste("./file/", id, "/dep_output_result.csv", sep=""), header=TRUE, fileEncoding ="UTF-8", sep = ',')
## feature 1: numeric vector
geneList = data[,ratioName]

## feature 2: named vector
names(geneList) = as.character(data[, 'human_entrez'])

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)
print(head(geneList))

r_geneList_de_up_down <- function(de_up_down, range, enrichment_analysis_methods, universal_enrichment_category, universal_enrichment_subcategory){
    print("--------------r_geneList_de_up_down(de_up_down_py, range_py) start--------------")
    cat("head(geneList)\n", head(geneList), "\n")
    print(head(geneList))

    print("r_geneList_de_up_down!!!!")
    if (de_up_down == "up") {
        print("up!")
        de <- names(geneList)[ geneList > range]
    }else if(de_up_down == "down") {
        print("down!")
        de <- names(geneList)[ geneList < range]
    }else{
        print("de!")
        de <- names(geneList)[ abs(geneList) > range]
    }
    cat("head(de)\n", head(de, 10), "\n")

    print(enrichment_analysis_methods)

    if (enrichment_analysis_methods == "DGN") {
        # origin
        print("DGN!")
        edo <- enrichDGN(de)
    }else if(enrichment_analysis_methods == "KEGG") {
        # 7 KEGG enrichment analysis
        print("KEGG!")
        edo <- enrichKEGG(gene         = de,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
    }else if(enrichment_analysis_methods == "WikiPathways") {
        # 8 WikiPathways analysis
        print("WikiPathways!")
        edo <- enrichWP(gene = de, organism = "Homo sapiens")
    }else if(enrichment_analysis_methods == "Reactome") {
        # 9 Reactome enrichment analysis
        print("Reactome!")
        edo <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)
    }else if(enrichment_analysis_methods == "Disease") {
        # 10 Disease enrichment analysis
        print("Disease!")
        edo <- enrichDO(gene          = de,
                        ont           = "DO",
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        universe      = names(geneList),
                        minGSSize     = 5,
                        maxGSSize     = 500,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)
    }else{
        # 12 Universal enrichment analysis
        print("Universal!")

        C3_t2g <- msigdbr(species = "Homo sapiens", category = universal_enrichment_category, subcategory = universal_enrichment_subcategory) %>%
            dplyr::select(gs_name, entrez_gene)
        edo <- enricher(gene=de, TERM2GENE=C3_t2g)
    }
    return (edo)
}

draw_plot2_1_to_2_7_2_9 <- function(edo, edo2, edox, geneList) {
    # plot2_1
    pic_name <- paste("./file/", id, "/image/plot2_1.png", sep="")
    tryCatch(
        {
            pic1 <- barplot(edo,x = "Count", color="p.adjust", showCategory=20)+ xlab("Count")
            save_plot(pic_name , pic1, base_height = 10, base_aspect_ratio = 1)
        },
        error = function(e) {
            message("Error!!!!! plot2_1.png ")
            print(e)
            if (file.exists(pic_name)) {
                file.remove(pic_name)
            }
        }
    )
    # plot2_2
    pic_name <- paste("./file/", id, "/image/plot2_2_1.png", sep="")
    tryCatch(
        {
            pic11_1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
            save_plot(pic_name, pic11_1, base_height = 10, base_aspect_ratio = 1)
        },
        error = function(e) {
            message("Error!!!!! plot2_2_1.png ")
            print(e)
            if (file.exists(pic_name)) {
                file.remove(pic_name)
            }
        }
    )
    tryCatch(
        {
            pic11_2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
            save_plot(paste("./file/", id, "/image/plot2_2_2.png", sep=""), pic11_2, base_height = 10, base_aspect_ratio = 1)
        },
        error = function(e) {
            message("Error!!!!! plot2_2_2.png ")
            print(e)
        }
    )
    # plot2_3
    tryCatch(
        {
            rescale.AsIs <- function(x, ...){
                dropAsis <- function(x){
                    cls <- class(x)
                    structure(x, class = setdiff(cls, "AsIs"))
                }
                scales:::rescale(dropAsis(x), ...)
            }

            pic12_1 <- cnetplot(edox,categorySize="geneNum",foldChange=geneList)
            pic12_2 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
            pic12_3 <- cnetplot(edox, node_label="category")
            pic12_4 <- cnetplot(edox, node_label="all")
            pic12_5 <- cnetplot(edox, node_label="none")

            save_plot(paste("./file/", id, "/image/plot2_3_1.png", sep=""), pic12_1, base_height = 10, base_aspect_ratio = 1)
            save_plot(paste("./file/", id, "/image/plot2_3_2.png", sep=""), pic12_2, base_height = 10, base_aspect_ratio = 1)
            save_plot(paste("./file/", id, "/image/plot2_3_3.png", sep=""), pic12_3, base_height = 10, base_aspect_ratio = 1)
            save_plot(paste("./file/", id, "/image/plot2_3_4.png", sep=""), pic12_4, base_height = 10, base_aspect_ratio = 1)
            save_plot(paste("./file/", id, "/image/plot2_3_5.png", sep=""), pic12_5, base_height = 10, base_aspect_ratio = 1)
        },
        error = function(e) {
            message("Error!!!!! plot2_3.png ")
            print(e)
        }
    )
    # plot2_4
    tryCatch(
        {
            pic13_1 <- heatplot(edox)
            pic13_2 <- heatplot(edox, foldChange=geneList)
            save_plot(paste("./file/", id, "/image/plot2_4_1.png", sep=""), pic13_1)
            save_plot(paste("./file/", id, "/image/plot2_4_2.png", sep=""), pic13_2)

            save_plot(paste("./file/", id, "/image/heatplot_1.png", sep=""), pic13_1, base_height = 10, base_aspect_ratio = 3,limitsize = FALSE)
            save_plot(paste("./file/", id, "/image/heatplot_2.png", sep=""), pic13_2, base_height = 10, base_aspect_ratio = 3,limitsize = FALSE)

        },
        error = function(e) {
            message("Error!!!!! plot2_4.png ")
            print(e)
        }
    )
    # plot2_5_1
    tryCatch(
        {
            edo <- pairwise_termsim(edo)
            pic14_1 <- emapplot(edo)
            pic14_2 <- emapplot(edo, cex_category=1.5)
            pic14_3 <- emapplot(edo, layout="kk")
            pic14_4 <- emapplot(edo, cex_category=1.5,layout="kk")

            save_plot(paste("./file/", id, "/image/plot2_5_1.png", sep=""), pic14_1, base_height = 10, base_aspect_ratio = 1)
            save_plot(paste("./file/", id, "/image/plot2_5_2.png", sep=""), pic14_2, base_height = 10, base_aspect_ratio = 1)
            save_plot(paste("./file/", id, "/image/plot2_5_3.png", sep=""), pic14_3, base_height = 10, base_aspect_ratio = 1)
            save_plot(paste("./file/", id, "/image/plot2_5_4.png", sep=""), pic14_4, base_height = 10, base_aspect_ratio = 1)
        },
        error = function(e) {
            message("Error!!!!! plot2_5_1.png ")
            print(e)
        }
    )
    # plot2_6
    tryCatch(
        {
            xx <- compareCluster(data, fun="enrichKEGG",
                                    organism="hsa", pvalueCutoff=0.05) #pvalueCutoff正常為0.05
                xx <- pairwise_termsim(xx)
                pic15_1 <- emapplot(xx)
                pic15_2 <- emapplot(xx, legend_n=2)
                pic15_3 <- emapplot(xx, pie="count")
                pic15_4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")

                save_plot(paste("./file/", id, "/image/plot2_6_1.png", sep=""), pic15_1, base_height = 10, base_aspect_ratio = 1)
                save_plot(paste("./file/", id, "/image/plot2_6_2.png", sep=""), pic15_2, base_height = 10, base_aspect_ratio = 1)
                save_plot(paste("./file/", id, "/image/plot2_6_3.png", sep=""), pic15_3, base_height = 10, base_aspect_ratio = 1)
                save_plot(paste("./file/", id, "/image/plot2_6_4.png", sep=""), pic15_4, base_height = 10, base_aspect_ratio = 1)
        },
        error = function(e) {
            message("Error!!!!! plot2_6.png ")
            print(e)
        }
    )
    # plot2_7
    tryCatch(
        {
            save_plot(paste("./file/", id, "/image/plot2_7.png", sep=""), upsetplot(edo), base_height = 10, base_aspect_ratio = 1.5)
        },
        error = function(e) {
            message("Error!!!!! plot2_7.png ")
            print(e)
        }
    )
    # plot2_9
    tryCatch(
        {
            pic17 <- ridgeplot(edo2) + xlab("expression distributions of enriched genes (log2FC)")
            save_plot(paste("./file/", id, "/image/plot2_9.png", sep=""), pic17, base_height = 12, base_aspect_ratio = 0.65)
        },
        error = function(e) {
            message("Error!!!!! plot2_9.png ")
            print(e)
        }
    )
    # plot2_10
}
edo <- r_geneList_de_up_down(de_up_down, range, enrichment_analysis_methods, universal_enrichment_category, universal_enrichment_subcategory)
edo2 <- gseDO(geneList, pvalueCutoff=1)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

draw_plot2_1_to_2_7_2_9(edo, edo2, edox, geneList)

pic18_1 <- gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
pic18_2 <- gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])
pic18_3 <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])
pic18_4 <- gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])
pic18_5 <- gseaplot2(edo2, geneSetID = 1:3)
pic18_6 <- gseaplot2(edo2, geneSetID = 1:3, pvalue_table = TRUE, color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
pic18_7 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1) + xlab("Rank")
pic18_8 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1:2)# + xlab("Rank")
pic18_9 <- gsearank(edo2, 1, title = edo2[1, "Description"])
pic18_10 <- lapply(1:3, function(i){
        anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
        lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
        gsearank(edo2, i, edo2[i, 2]) + xlab("Rank")+
        annotate("text", 5000, edo2[i, "enrichmentScore"]* .75, label =lab, hjust=0, vjust=0)})


save_plot(paste("./file/", id, "/image/plot2_10_1.png", sep=""), pic18_1, base_height = 10, base_aspect_ratio = 1)
save_plot(paste("./file/", id, "/image/plot2_10_2.png", sep=""), pic18_2, base_height = 10, base_aspect_ratio = 1)

png(paste("./file/", id, "/image/plot2_10_3.png", sep=""))
pic18_3
dev.off()

png(paste("./file/", id, "/image/plot2_10_4.png", sep=""))
pic18_4
dev.off()

png(paste("./file/", id, "/image/plot2_10_5.png", sep=""))
pic18_5
dev.off()

png(paste("./file/", id, "/image/plot2_10_6.png", sep=""))
pic18_6
dev.off()

png(paste("./file/", id, "/image/plot2_10_7.png", sep=""))
pic18_7
dev.off()

png(paste("./file/", id, "/image/plot2_10_8.png", sep=""))
pic18_8
dev.off()

png(paste("./file/", id, "/image/plot2_10_9.png", sep=""))
pic18_9
dev.off()

png(paste("./file/", id, "/image/plot2_10_10.png", sep=""))
pic18_10
dev.off()