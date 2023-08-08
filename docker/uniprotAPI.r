isJobReady <- function(jobId) {
    pollingInterval = 15
    nTries = 60
    for (i in 1:nTries) {
        url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
        r <- GET(url = url, accept_json())
        status <- httr::content(r, as = "parsed")
        if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
            return(TRUE)
        }
        if (!is.null(status[["messages"]])) {
            print(status[["messages"]])
            return (FALSE)
        }
            Sys.sleep(pollingInterval)
    }
        return(FALSE)
}
getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
        url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
        url <- gsub("/results/", "/results/stream/", redirectURL)
    }
}

generate_uniprot_resultTable <- function(ID_string){
    result_1 <- POST("https://rest.uniprot.org/idmapping/run",
                add_headers("user-agent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36"),
                body = list(ids =ID_string , 'from'="UniProtKB_AC-ID",'to'="GeneID")) #c("C0HK80,C0HKD9")
    submission <- httr::content(result_1, as = "parsed")
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- httr::content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])

    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable = read.table(text = httr::content(r), sep = "\t", header=TRUE, col.names=c('ID_uniprot', 'GeneID_entrez'), encoding = "UTF-8")
    return (resultsTable)
}