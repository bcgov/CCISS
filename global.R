if (!requireNamespace("Require")) {
  remotes::install_github("PredictiveEcology/Require")
}

## set lib path to project lib
invisible(Require::setLibPaths("packages/")) 

## install/load packages
Require::Require(c(
  "bcgov/climr@devl", 
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
  "terra", 
  "themis",
  "tidymodels"
))


options(reproducible.cachePath = normalizePath("reproducible.cache/", winslash = "/"),
        reproducible.destinationPath = normalizePath("data/", winslash = "/"),
        climr.cache.path = "climr.cache/")

## source functions
source("R/utils.R")
