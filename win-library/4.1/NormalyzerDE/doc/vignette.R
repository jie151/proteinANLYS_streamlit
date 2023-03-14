## ----setup, echo=FALSE, results="hide"----------------------------------------
# knitr::opts_chunk$set(tidy=FALSE, cache=TRUE, dev="png", message=FALSE,
# error=FALSE, warning=TRUE)

## -----------------------------------------------------------------------------
library(NormalyzerDE)
outDir <- tempdir()
designFp <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
dataFp <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
normalyzer(jobName="vignette_run", designPath=designFp, dataPath=dataFp, 
           outputDir=outDir)

## -----------------------------------------------------------------------------
normMatrixPath <- paste(outDir, "vignette_run/CycLoess-normalized.txt", sep="/")
normalyzerDE("vignette_run", 
             comparisons=c("4-5"),
             designPath=designFp, 
             dataPath=normMatrixPath,
             outputDir=outDir, 
             condCol="group")

## -----------------------------------------------------------------------------
dataMatrix <- read.table(dataFp, sep="\t", header = TRUE)
designMatrix <- read.table(designFp, sep="\t", header = TRUE)
designMatrix$sample <- as.character(designMatrix$sample)
dataOnly <- dataMatrix[, designMatrix$sample]
annotOnly <- dataMatrix[, !(colnames(dataMatrix) %in% designMatrix$sample)]

sumExpObj <- SummarizedExperiment::SummarizedExperiment(
    as.matrix(dataOnly),
    colData=designMatrix,
    rowData=annotOnly
)

normalyzer(jobName="sumExpRun", experimentObj = sumExpObj, outputDir=outDir)

## -----------------------------------------------------------------------------
fullDf <- read.csv(dataFp, sep="\t")
designDf <- read.csv(designFp, sep="\t")
head(fullDf, 1)
head(designDf, 1)

## -----------------------------------------------------------------------------
sampleNames <- as.character(designDf$sample)
typeof(sampleNames)

## -----------------------------------------------------------------------------
dataMat <- as.matrix(fullDf[, sampleNames])
retentionTimes <- fullDf$Average.RT

head(dataMat, 1)

## -----------------------------------------------------------------------------
typeof(dataMat)

print("Rows and columns of data")
dim(dataMat)

print("Number of retention times")
length(retentionTimes)

## -----------------------------------------------------------------------------
performCyclicLoessNormalization <- function(rawMatrix) {
    log2Matrix <- log2(rawMatrix)
    normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}

## -----------------------------------------------------------------------------
rtNormMat <- getRTNormalizedMatrix(dataMat, 
                                   retentionTimes, 
                                   performCyclicLoessNormalization, 
                                   stepSizeMinutes=1, 
                                   windowMinCount=100)

## -----------------------------------------------------------------------------
globalNormMat <- performCyclicLoessNormalization(dataMat)
dim(rtNormMat)
dim(globalNormMat)
head(rtNormMat, 1)
head(globalNormMat, 1)

## -----------------------------------------------------------------------------
layeredRtNormMat <- getSmoothedRTNormalizedMatrix(
    dataMat, 
    retentionTimes, 
    performCyclicLoessNormalization, 
    stepSizeMinutes=1, 
    windowMinCount=100, 
    windowShifts=3, 
    mergeMethod="mean")

dim(layeredRtNormMat)
head(layeredRtNormMat, 1)

## -----------------------------------------------------------------------------
jobName <- "vignette_run"
experimentObj <- setupRawDataObject(dataFp, designFp, "default", TRUE, "sample", "group")
normObj <- getVerifiedNormalyzerObject(jobName, experimentObj)

## -----------------------------------------------------------------------------
normResults <- normMethods(normObj)

## -----------------------------------------------------------------------------
normResultsWithEval <- analyzeNormalizations(normResults)

## -----------------------------------------------------------------------------
jobDir <- setupJobDir("vignette_run", tempdir())
writeNormalizedDatasets(normResultsWithEval, jobDir)

## -----------------------------------------------------------------------------
generatePlots(normResultsWithEval, jobDir)

## -----------------------------------------------------------------------------
bestNormMatPath <- paste(jobDir, "RT-Loess-normalized.txt", sep="/")
experimentObj <- setupRawContrastObject(bestNormMatPath, designFp, "sample")
nst <- NormalyzerStatistics(experimentObj, logTrans=FALSE)

## -----------------------------------------------------------------------------
comparisons <- c("4-5")
nst <- calculateContrasts(nst, comparisons, condCol="group", leastRepCount=2)

## -----------------------------------------------------------------------------
annotDf <- generateAnnotatedMatrix(nst)
utils::write.table(annotDf, file=paste(jobDir, "stat_table.tsv", sep="/"))
generateStatsReport(nst, "Vignette stats", jobDir)

## -----------------------------------------------------------------------------
sessionInfo()

