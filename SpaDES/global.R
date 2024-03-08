if (!requireNamespace("Require")) {
  remotes::install_github("PredictiveEcology/Require")
}
## set lib path to project lib
invisible(Require::setLibPaths("packages/")) 

Require::Require("PredictiveEcology/SpaDES.project@transition (HEAD)")

## install/load packages
Require::Require(c(
  "bcgov/climr@devl (HEAD)", 
  "data.table",
  "foreach",
  "ggplot2", 
  "reproducible",
  "smotefamily", 
  "mlr3learners",
  "mlr3spatial", 
  "mlr3spatiotempcv",
  "mlr3viz",
  "terra", 
  "themis",
  "tidymodels"
))

## install, but don't load these.
Require::Require(c(
  "aws.s3",
  "future",
  "sf"), 
  require = FALSE)

## in 03_RunPredictHex -- not sure what we'll need
# require(foreach)
# require(dplyr)
# require(reshape2)
# library(doParallel)
# library(tidyr)
# require(sf)
# require(RPostgreSQL)
# library(disk.frame)
# require(RPostgres)

options(reproducible.cachePath = normalizePath("reproducible.cache/", winslash = "/"),
        reproducible.destinationPath = normalizePath("data/", winslash = "/"),
        climr.cache.path = "climr.cache/")

## source functions
source("SpaDES/R/utils.R")

## source scripts
source("R/02_Build_WNA_BGC_trainingset.R")
source("R/03_RunPredictHex.R")
source("R/04_CreateBGCfutsMap.R")