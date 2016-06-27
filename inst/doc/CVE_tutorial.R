## ----eval=FALSE----------------------------------------------------------
#  source('http://www.bioconductor.org/biocLite.R')
#  biocLite('CVE')

## ----eval=FALSE----------------------------------------------------------
#  library(CVE)

## ---- echo=F,results='hide', eval=F--------------------------------------
#  devtools::load_all()

## ---- eval=FALSE---------------------------------------------------------
#  openCVE(oncotator_example, sample_names='case study')

## ---- eval=FALSE---------------------------------------------------------
#  openCVE(oncotator_example, sample_names='case study', extension='WGCNAmelanoma')

## ------------------------------------------------------------------------
sessionInfo()

