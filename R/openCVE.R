#' Open Cancer Variant Explorer (CVE) Shiny app
#' @description The openCVE function opens the CVE Shiny app.
#' The function to supplement the R package with the Shiny app was suggested by Dean Attali (http://deanattali.com). Currently, the only extension available is a melanoma co-expression network (WGCNAmelanoma).
#' @param x A dataframe (for single file) or list (for multiple oncotator output files)
#' @param sample_names A character vector with sample name(s)
#' @param extension A character vector of extention name
#' @examples
#' openCVE(oncotator_example, sample_names="case study")
#' openCVE(oncotator_example, sample_names="case study",extension="WGCNAmelanoma")
openCVE <- function(x, sample_names, extension=FALSE) {
  #set app directory
  if(extension=="WGCNAmelanoma"){
    appDir <-system.file("Shiny","CVE_WGCNA_melanoma",package="CVE")
  } else {
    appDir <- system.file("Shiny", "CVE", package = "CVE")
  }
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `CVE`.", call. = FALSE)
  }

  #if data frame -> turn into list
  if(class(x)=="data.frame"){
    z = vector("list",1)
    z[[1]] = x
    names(z) = sample_names
    v <<-z
  } else if (class(x)=="list"){
    z = x
    names(z) = sample_names
    v <<-z
  } else {
    stop("x is neither a data frame nor list. Please change format accordingly.")
  }

  shiny::runApp(appDir, display.mode = "normal")
}

