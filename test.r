packageVersion("biomaRt")
"C:/Users/jie/Desktop/pyEnv/web_proteinAnls/code/win-library/4.1/"


install.packages("BiocManager", lib = path)
BiocManager::install("DEP", lib = path)

.libPaths(path)

#library(biomaRt, lib=lib_path)
#library(clusterProfiler, lib=lib_path)
library(DEP, lib=lib_path)
#library(DOSE, lib=lib_path)
#library(enrichplot, lib=lib_path)
#library(NormalyzerDE, lib=lib_path)
#library(SummarizedExperiment, lib=lib_path)