if (!requireNamespace("Require")) {
  remotes::install_github("PredictiveEcology/Require")
}

## set lib path to project lib
invisible(Require::setLibPaths("packages/")) 

## install/load packages
Require::Require(c(
  "bcgov/climr@devl (HEAD)", 
  "caret", 
  "clhs", 
  "data.table", 
  "foreach", 
  "ggplot2", 
  "randtoolbox", 
  "ranger", 
  "reproducible",
  "sf", 
  "smotefamily", 
  "mlr3learners",
  "mlr3spatial", 
  "mlr3spatiotempcv",
  "mlr3viz",
  "rmapshaper",
  "terra", 
  "themis",
  "tidymodels"
))

## install, but don't load these.
Require::Require(c(
  "aws.s3",
  require = FALSE)

options(reproducible.cachePath = normalizePath("reproducible.cache/", winslash = "/"),
        reproducible.destinationPath = normalizePath("data/", winslash = "/"),
        climr.cache.path = "climr.cache/")

## source functions
source("R/utils.R")

## source scripts
source("R/02_Build_WNA_BGC_trainingset.R")
source("R/03_RunPredictHex.R")
source("R/04_CreateBGCfutsMap.R")