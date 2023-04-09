import os

os.environ['R_HOME']= 'C:\\Program Files\\R\\R-4.2.3'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\bin\\X64\\'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\'


from rpy2 import robjects
import streamlit as st
import pandas as pd
import os
from PIL import Image
from contextlib import contextmanager, redirect_stdout
from io import StringIO
import streamlit.components.v1 as components


def r_test_file():
    robjects.r('''
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
        data(geneList)
        de <- names(geneList)[abs(geneList) > 2]
        edo <- enrichDGN(de)
        edo2 <- gseDO(geneList, pvalueCutoff=1)
        edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

    ''')

r_test_file()

def r_plot_gseaplot_2_10():
    robjects.r('''

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

        png("./test.png")
        pic18_10 <- lapply(1:3, function(i){
            anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
            lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
            gsearank(edo2, i, edo2[i, 2]) + xlab("Rank")+
            annotate("text", 5000, edo2[i, "enrichmentScore"]* .75, label =lab, hjust=0, vjust=0)})
        dev.off()

        save_plot("./file/image/plot2_10_1.png", pic18_1, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_2.png", pic18_2, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_3.png", pic18_3, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_4.png", pic18_4, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_5.png", pic18_5, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_6.png", pic18_6, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_7.png", pic18_7, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_8.png", pic18_8, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_9.png", pic18_9, base_height = 10, base_aspect_ratio = 1)
        save_plot("./plot2_10_10.png", pic18_10, base_height = 10, base_aspect_ratio = 1)
    ''')
    for i in range(1, 11):
        st.subheader(i)
        st.image(Image.open(f"./plot2_10_{i}.png"))

r_plot_gseaplot_2_10()

st.write("HI")