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

data <- UbiLength
data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Make SummarizedExperiment
columns <- grep("LFQ.", colnames(data_unique))
exp_design <- UbiLength_ExpDesign
data <- read.csv(file="2.proteinGroups_injury.txt", header=TRUE, fileEncoding ="UTF-8", sep = '\t')
cat("file's row * column =", dim(data), "\n")
colname_proteinIDs <- "Protein.IDs"
colname_geneNames <- "Gene.names"
cat('* Are there any duplicated gene names? ', data[, colname_geneNames] %>% duplicated() %>% any(), "\n")

                if ( data[, colname_geneNames] %>% duplicated() %>% any() ){
                    # Make a table of duplicated gene names
                    write.table(data %>% group_by_(.dots = colname_geneNames) %>% summarize(frequency = n()) %>%
                        arrange(desc(frequency)) %>% filter(frequency > 1), file = "./file/my_data1.txt", row.names =FALSE)

                    cat("a table of duplicated gene names: (table ", dim(table), "\n")
                    table <- read.csv(file="./file/my_data1.txt", header=TRUE, fileEncoding ="UTF-8", sep = ' ')
                    print(head(table, 7))
                }

                # Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
                data_unique <- make_unique(data, colname_geneNames, colname_proteinIDs, delim = ";")

                # Generate a SummarizedExperiment object using an experimental design
                LFQ_columns <- grep("Reporter.intensity.corrected.", colnames(data_unique)) # get LFQ column numbers

                experimental_design <- read.csv('2.exp_design_injury.csv',header=TRUE ,fileEncoding ="UTF-8")
                experimental_design$label = gsub(" ", ".", experimental_design$label)

                maxReplicate <- max(experimental_design$replicate) # max replicate
                experimental_design_condition <- unique(experimental_design$condition)

                data_se <- make_se(data_unique, LFQ_columns, experimental_design)

                # Let's have a look at the SummarizedExperiment object
                cat("* the SummarizedExperiment object:  \n")
                print(data_se)

proteins_unique <- data_unique
columns <- LFQ_columns
expdesign <- experimental_design
