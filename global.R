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

## for 03 but to simplify.
require (smotefamily)
library(foreign)
library(RStoolbox)
library(maptools)
library(spatstat)
require(scales)
require(caret)
require(gstat)
require(purrr)
require(forcats)
require(StatMatch)
require(lwgeom)

## install, but don't load these.
Require::Require(c(
  "future"), 
  require = FALSE)

options(reproducible.cachePath = normalizePath("reproducible.cache/", winslash = "/"),
        reproducible.destinationPath = normalizePath("data/", winslash = "/"),
        climr.cache.path = "climr.cache/")

## source functions
source("R/utils.R")
