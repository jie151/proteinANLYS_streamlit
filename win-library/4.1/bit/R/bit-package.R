# Package documentation
# (c) 2008-2017 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk

#' bit: Classes and methods for fast memory-efficient boolean selections
#' 
#' Provided are classes for boolean and skewed boolean vectors, fast boolean 
#' methods, fast unique and non-unique integer sorting, fast set operations on 
#' sorted and unsorted sets of integers, and foundations for ff (range indices, 
#' compression, chunked processing).
#' 
#' For details view the vignettes \url{../doc/bit-usage.pdf} and
#' \url{../doc/bit-performance.pdf}
#'
#'@name bit-package
NULL

# devtools::use_vignette("bit-usage")
# devtools::use_vignette("bit-performance")

# require(rhub)
# rhub_bit_4.0.5 <- check_for_cran(
#   path = "../bit_4.0.5.tar.gz"
# , email = "Jens.Oehlschlaegel@truecluster.com"
# , check_args = "--as-cran"
# , env_vars = c('_R_CHECK_FORCE_SUGGESTS_'= "false",'_R_CHECK_CRAN_INCOMING_USE_ASPELL_'= "true", '_R_CHECK_XREFS_MIND_SUSPECT_ANCHORS_'="true")
# , platforms = NULL
# , show_status = FALSE
# )
# 
# ─  Uploading package
# ─  Preparing build, see status at
# https://builder.r-hub.io/status/bit_4.0.5.tar.gz-b6a1c13ef92c42d9b68257139666fb46
# https://builder.r-hub.io/status/bit_4.0.5.tar.gz-abc80afc9741404c9e86511dee4585c5
# https://builder.r-hub.io/status/bit_4.0.5.tar.gz-93409d19ee544acab83be7de1b380740
# https://builder.r-hub.io/status/bit_4.0.5.tar.gz-d4f141780d004fa388c82d027a6e40e2

# olddir <- "../revdepold"
# newdir <- "../revdepnew"
# tools::check_packages_in_dir(olddir,
#                              check_args = c("--as-cran", ""),
#                              reverse = list(repos = getOption("repos")["CRAN"]))
# tools::check_packages_in_dir(newdir, old=olddir
#                              check_args = c("--as-cran", ""),
#                              reverse = list(repos = getOption("repos")["CRAN"]))
# tools::summarize_check_packages_in_dir_results(newdir, all = FALSE, full = TRUE)
# tools::check_packages_in_dir_changes(newdir, olddir, outputs = TRUE, sources = FALSE)
