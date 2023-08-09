library(biomaRt)
library(httr)
library(dplyr)
library(msigdbr)
source("uniprotAPI.r")

print("--------------dose1_config.r start--------------")
# ----------參數----------
Args <- commandArgs(TRUE)
id       <- Args[1]
species  <- Args[2]
ratioName <- Args[3]
# ----------參數----------

collections <- msigdbr_collections()
write.table(collections, file=paste("./file/", id, "/msigdbr_collections.csv", sep=""), sep=",",  row.names=FALSE)

dep_output_data = read.csv(file = paste("./file/", id, "/dep_output.csv", sep=""), header=TRUE, fileEncoding ="UTF-8")
dep_output_data[,'ID'] <- sub("-.*", "", dep_output_data$ID)

ratioCol <- grep( "ratio$" , colnames(dep_output_data) )
ratioColname <- colnames(dep_output_data[ratioCol])
numRatio <- length(ratioCol)

ID_string <- paste(dep_output_data$ID, collapse=",")
resultsTable <- generate_uniprot_resultTable(ID_string)

uniprot_entrez <- resultsTable[-c(1),]
#去除重複 (因為uniprot可能會有多個entrez => 排序後選最小的)
uniprot_entrez[, 2] <- as.numeric(uniprot_entrez[, 2])
uniprot_entrez <- uniprot_entrez [ order(uniprot_entrez$ID_uniprot, uniprot_entrez$GeneID_entrez),]
uniprot_entrez <- uniprot_entrez[!duplicated(uniprot_entrez$ID_uniprot),]
write.csv(uniprot_entrez, file="./file/uniprot_entrez.csv")

if(species != "human"){
    human <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
    print("******connect to ensembl (human)******")
    if(species == "mouse"){
        mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        print("******connect to ensembl => mouse******")
    }else{
        rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        print("******connect to ensembl => rat******")
    }
}

if(species != "human"){
    print("******not human******")
    if(species == "mouse"){
        human_mouse_rat <- getLDS(mart=mouse, attributes=c("entrezgene_id"), filters="entrezgene_id" , values=uniprot_entrez[,"GeneID_entrez"], attributesL= c("entrezgene_id"), martL = human, uniqueRows=T ) #values=uniprots_ID
        print("******species => mouse******")
    }else{
        human_mouse_rat <- getLDS(mart=rat, attributes=c("entrezgene_id"), filters="entrezgene_id" , values=uniprot_entrez[,"GeneID_entrez"], attributesL= c("entrezgene_id"), martL = human, uniqueRows=T ) #values=uniprots_ID
        print("******species => rat******")
    }
    colnames(human_mouse_rat) <- c('GeneID_entrez','human')
    human_mouse_rat <-human_mouse_rat[ order(human_mouse_rat$GeneID_entrez, human_mouse_rat$human),]
    human_mouse_rat <- human_mouse_rat[!duplicated(human_mouse_rat$GeneID_entrez),]

    temp <- inner_join( uniprot_entrez, human_mouse_rat, by="GeneID_entrez")
    colnames(temp) <- c('ID', 'GeneID', 'human_entrez')
    temp[, "ID"] <- as.character(temp[, "ID"])

    data <- inner_join(temp, dep_output_data, by="ID")
    print(dim(data))
    print(head(data))
}else{
    print("species => human******")
    colnames(uniprot_entrez) <- c("ID", "human_entrez")
    data <- inner_join(uniprot_entrez, dep_output_data, by="ID")
}
write.csv(data, file= paste("./file/", id,"/dep_output_result.csv", sep=""))