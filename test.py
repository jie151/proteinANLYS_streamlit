from rpy2 import robjects
import streamlit as st
import pandas as pd
import os
from PIL import Image
from contextlib import contextmanager, redirect_stdout
from io import StringIO

os.system("lsb_release -a")
if not os.path.exists("./library/"):
    os.system("mkdir -m 777 library")

    robjects.r('''
        path = "./library"
        .libPaths(path)

        install.packages("ggupset", lib=path)
        install.packages("ggridges", lib=path)
        install.packages("rlang",lib=path)
        install.packages("BiocManager", repos = "http://cran.us.r-project.org", lib = path)
        library(BiocManager, lib = path)
        BiocManager::install("biomaRt", lib = path)

        install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz",repos = NULL,type = "source", lib=path)
        BiocManager::install("clusterProfiler", lib = path, force=TRUE, update=FALSE, ask=FALSE)

        BiocManager::install("mzR", lib=path, force=TRUE, update=FALSE, ask=FALSE)
        BiocManager::install("DEP", lib = path, force=TRUE, update=FALSE, ask=FALSE)
        BiocManager::install("SummarizedExperiment", lib = path, force=TRUE, update=FALSE, ask=FALSE)
        BiocManager::install("DOSE", lib = path, force=TRUE, update=FALSE, ask=FALSE)
        BiocManager::install("enrichplot", lib = path, force=TRUE, update=FALSE, ask=FALSE)
        BiocManager::install("NormalyzerDE", lib = path, force=TRUE, update=TRUE, ask=FALSE)
        BiocManager::install("org.Hs.eg.db", lib = path, force=TRUE, update=TRUE, ask=FALSE)
    ''')

# Load package in r
robjects.r('''
    path = "./library"
    library(cowplot) # save_plot
    library(dplyr)
    library(ggplot2)
    library(httr)

    .libPaths(path)
    library(biomaRt, lib = path)
    library(clusterProfiler, lib = path)
    library(DEP, lib = path)
    library(DOSE, lib = path)
    library(enrichplot, lib = path)
    library(NormalyzerDE, lib = path)
    library(SummarizedExperiment, lib = path)

''')

robjects.r('''

    data(geneList)
    de <- names(geneList)[abs(geneList) > 2]
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
        pic18_8 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1:2) + xlab("Rank")
        #pic9_7_8 <- cowplot::plot_grid(pic9_7, pic9_8, ncol=1, labels=LETTERS[1:2])

        pic18_9 <- gsearank(edo2, 1, title = edo2[1, "Description"])

        pic18_10 <- lapply(1:3, function(i){
            anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
            lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
            gsearank(edo2, i, edo2[i, 2]) + xlab("Rank")+
            annotate("text", 5000, edo2[i, "enrichmentScore"]* .75, label =lab, hjust=0, vjust=0)})

    save_plot("./file/image/plot2_10_1.png", pic18_1, base_height = 10, base_aspect_ratio = 1)
    save_plot("./file/image/plot2_10_2.png", pic18_2, base_height = 10, base_aspect_ratio = 1)
    save_plot("./file/image/plot2_10_7.png", pic18_7, base_height = 10, base_aspect_ratio = 1)
    save_plot("./file/image/plot2_10_9.png", pic18_9, base_height = 10, base_aspect_ratio = 1)

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

    png("./file/image/plot2_10_8.png")
    pic18_8
    dev.off()
''')
for i in range(1, 10):
    st.image(Image.open(f"./file/image/plot2_10_{i}.png"))

st.write("DONE")