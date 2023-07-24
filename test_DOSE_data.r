library(DEP)
library(cowplot) # save_plot
library(dplyr)
library(ggplot2)
library(httr)
library(clusterProfiler)
library(DEP)
library(DOSE)
library(enrichplot)
library(NormalyzerDE)
library(SummarizedExperiment)
library(biomaRt)
library(ReactomePA)
library(AnnotationHub)
library(MeSHDbi)
library(msigdbr)
library(meshes)

# data <- read.csv('/app/file/dep_output_result.csv',header=TRUE ,fileEncoding ="UTF-8")
data <- read.csv('./_file_dep_output_result (1).csv',header=TRUE ,fileEncoding ="UTF-8")

## feature 1: numeric vector
geneList = data[, "B_vs_A_ratio"]

## feature 2: named vector
names(geneList) = as.character(data[, 'human_entrez'])

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)
print(head(geneList))

de_up_down <- "de"
range <- "0.5"
enrichment_analysis_methods <- "Universal"

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

    C3_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "MIR:MIRDB") %>%
                dplyr::select(gs_name, entrez_gene)
    edo <- enricher(gene=de, TERM2GENE=C3_t2g)
}
edo <- enrichDGN(de)
edo2 <- gseDO(geneList, pvalueCutoff=1)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')