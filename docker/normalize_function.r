
logNorm <- function(dat) {
    logInt <- log2(dat)
    #logInt <- replace(is.infinite(logInt), NA)
    logInt[is.infinite(as.matrix(logInt))] <- NA
    return(as.matrix(logInt))
}
medianNorm <- function(logDat) {
    # Find medians of each sample
    # Divide by median
    # Multiply by mean of medians
    sampleMed <- apply(logDat, 2, median, na.rm=TRUE)
    meanMed <- mean(sampleMed, na.rm=TRUE)
    out <- t(t(logDat) / sampleMed)
    out <- out * meanMed
    return(as.matrix(out))
}
meanNorm <- function(logDat) {
    # Find means of each sample
    # Divide by mean
    # Multiply by mean of means
    sampleMean <- apply(logDat, 2, mean, na.rm=TRUE)
    meanMean <- mean(sampleMean, na.rm=TRUE)
    out <- t(t(logDat) / sampleMean)
    out <- out * meanMean
    return(as.matrix(out))
}
vsnNorm <- function(dat) {
    vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
    colnames(vsnNormed) <- colnames(dat)
    row.names(vsnNormed) <- rownames(dat)
    return(as.matrix(vsnNormed))
}
quantNorm <- function(logDat) {
    quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy=FALSE)
    colnames(quantNormed) <- colnames(logDat)
    row.names(quantNormed) <- rownames(logDat)
    return(as.matrix(quantNormed))
}
cycLoessNorm <- function(logDat) {
    cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method="fast")
    colnames(cycLoessNormed) <- colnames(logDat)
    row.names(cycLoessNormed) <- rownames(logDat)
    return(as.matrix(cycLoessNormed))
}
rlrNorm <- function(logDat) {
    rlrNormed <- NormalyzerDE::performGlobalRLRNormalization(as.matrix(logDat), noLogTransform=TRUE)
    colnames(rlrNormed) <- colnames(logDat)
    row.names(rlrNormed) <- rownames(logDat)
    return(as.matrix(rlrNormed))
}
giNorm <- function(logDat) {
    giNormed <- NormalyzerDE::globalIntensityNormalization(as.matrix(logDat), noLogTransform=TRUE)
    colnames(giNormed) <- colnames(logDat)
    row.names(giNormed) <- rownames(logDat)
    return(as.matrix(giNormed))
}

normalize_proteiNorm <- function(se, normalizeOption) {
    # Show error if inputs are not the required classes
    assertthat::assert_that(inherits(se, "SummarizedExperiment"))

    # Variance stabilization transformation on assay data
    se_vsn <- se

    if(normalizeOption == "VSN")

        assay(se_vsn)  <- vsnNorm(2 ^ assay(se_vsn))

    else if(normalizeOption == "Log2")

        assay(se_vsn)  <- logNorm(2 ^ assay(se_vsn))

    else if(normalizeOption == "Median")

        assay(se_vsn)  <- medianNorm(logNorm(2 ^ assay(se_vsn)))

    else if(normalizeOption == "Mean")

        assay(se_vsn)  <- meanNorm(logNorm(2 ^ assay(se_vsn)))

    else if(normalizeOption == "Quantile")

        assay(se_vsn)  <- quantNorm(logNorm(2 ^ assay(se_vsn)))

    else if(normalizeOption == "Cyclic Loess")

        assay(se_vsn)  <- cycLoessNorm(logNorm(2 ^ assay(se_vsn)))

    else if(normalizeOption == "RLR")

        assay(se_vsn) <- rlrNorm(logNorm(2 ^ assay(se_vsn)))

    else
        assay(se_vsn)  <- giNorm(logNorm(2 ^ assay(se_vsn)))

    return(se_vsn)
    }