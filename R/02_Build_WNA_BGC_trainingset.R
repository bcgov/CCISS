## This script must be sourced

## FIT AND EVALUATE A BGC MODEL

## TODO: need proper documentation and references to sources/scripts that produced these objects
## for now, download files manually from object storage and put them in:
## filepath(options(reproducible.destinationPath), "WNA_BGC_v12_5Apr2022.gpkg")
## filepath(options(reproducible.destinationPath), "northamerica/northamerica_elevation_cec_2023.tif")
## then select 'y' when asked if objects are in the right place.

## TODO: MakeGrd. is getting bgcs and dem, but only at bc scale. get them again.
bgcs <- reproducible::prepInputs(url = "//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022.gpkg",
                                 targetFile = "WNA_BGC_v12_5Apr2022.gpkg",
                                 fun = "sf::st_read")

# elev <- reproducible::prepInputs(url = "//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif",
#                                  targetFile = "northamerica/northamerica_elevation_cec_2023.tif")
elev <- reproducible::prepInputs(url = "//objectstore2.nrs.bcgov/ffec/DEM/DEM_30m_WNA/WNA_DEM_SRT_30m_cropped.tif",
                                 targetFile = "DEM_30m_WNA/WNA_DEM_SRT_30m_cropped.tif")
coords <- makePointCoords(bgcs, elev) |>
  Cache()
coords <- coords[!is.na(elev),]

## full training area
trainingarea <- ext(c(-125, -112, 43, 55))
coords_train <- subsetByExtent(coords, trainingarea)
setnames(coords_train, c("lon", "lat", "elev", "id"))

## this is no longer necessary with automated spatial CV
## make subsets of the study area for hold-outs (gaps)
# studyarea <- ext(c(-123, -117, 49, 52.5))
# gapextents <- makeGapExtents(studyarea, 5L)
# ## convert to poly
# gap_poly <- lapply(gapextents, vect, crs = "EPSG:4326")
# gap_poly <- do.call(rbind, gap_poly)
# 
# coords_gaps <- subsetByExtent(coords_train, gap_poly)
# coords_trainWgaps <- coords_train[!coords_gaps, on = "id"]

vars_needed <- c("DD5", "DD_0_at", "DD_0_wt", "PPT05", "PPT06", "PPT07", "PPT08",
                 "PPT09", "CMD", "PPT_at", "PPT_wt", "CMD07", "SHM", "AHM", "NFFD", "PAS", "CMI")

clim_vars <- getClimate(coords_train, bgcs,
                        which_normal = "normal_composite", return_normal = TRUE, 
                        vars = vars_needed, cache = TRUE) |>
  Cache()

## subset coords objects to ids with data
coords_train <- coords_train[clim_vars[, .(id)], on = "id", nomatch = 0L]
# coords_trainWgaps <- coords_trainWgaps[clim_vars[, .(id)], on = "id", nomatch = 0L]

setDT(clim_vars)   ## this shouldn't be necessary, submit issue/reprex to reproducible.
addVars(clim_vars)

vs_final <- c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
              "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS")

trainData <- clim_vars[, .SD, .SDcols = c("BGC", "id", vs_final)]
trainData <- trainData[complete.cases(trainData)]

# BGC_counts <- trainData[, .(Num = .N), by = .(BGC)]   ## for inspection


## clean trainind dataset - remove bad BGCs, outliers and BGCs with low sample sizes.
badbgcs <- c("BWBSvk", "ICHmc1a", "MHun", "SBSun", "ESSFun", "SWBvk", "MSdm3",
             "ESSFdc3", "IDFdxx_WY", "MSabS", "FGff", "JPWmk_WY" ) #, "ESSFab""CWHws2", "CWHwm", "CWHms1" , 
trainData <- trainData[!BGC %in% badbgcs,]

trainData <- removeOutlier(as.data.frame(trainData), alpha = .025, vars = vs_final) |>
  Cache()

trainData <- rmLowSampleBGCs(trainData) |>
  Cache()


## subsampling
dataBalance_recipe <- recipe(BGC ~ ., data =  trainData) |>
  step_downsample(BGC, under_ratio = 90) |>  ## subsamples "oversampled" BGCs
  prep()

## extract data.table
trainData_balanced <- dataBalance_recipe |>
  juice() |>
  as.data.table()

# trainData_balanced[,.(Num = .N), by = BGC]   ## for inspection

trainData_balanced[, BGC := as.factor(BGC)]

## FIT MODEL WITH FULL DATA SET OR CROSS-VALIDATION
cols <- c("BGC", vs_final)

# trainData_balanced2 <- trainData_balanced[coords_trainWgaps, on = "id", nomatch = 0L]   ## better to evaluate with CV (below)

coords_train_balanced <- coords_train[trainData_balanced, on = "id"]
coords_train_balanced <- st_as_sf(coords_train_balanced, coords = c("lon", "lat"))
st_crs(coords_train_balanced) <- crs(elev, proj = TRUE)

## make a modelling task
tsk_bgc <- as_task_classif_st(coords_train_balanced, target = "BGC")

## choose model type and eval metric
lrn_rf <- lrn("classif.ranger", predict_type = "response",
              num.trees = 501,
              splitrule =  "extratrees",
              mtry = 4,
              min.node.size = 2,
              importance = "permutation",
              write.forest = TRUE)
measure_acc <- msrs(c("classif.acc", "classif.ce", "oob_error"))

## fit full model --------------------------------------------
lrn_rf$train(tsk_bgc)
lrn_rf$model
## this plot is useless, but code is kept here for fut reference
# autoplot(lrn_rf$predict(tsk_bgc)) +
#   theme(legend.position = "none")

## evaluate model with CV ------------------------------------
## define and apply cv strategy
folds <- 10
cv_strategy <- rsmp("repeated_spcv_coords", folds = folds, repeats = 1)

future::plan("multisession", 
             workers = ifelse(folds <= future::availableCores(), folds, future::availableCores()))
RF_cv <- mlr3::resample(tsk_bgc, lrn_rf, cv_strategy, store_models = TRUE) |>
  Cache()
future:::ClusterRegistry("stop")

RF_cv$score(measure_acc)  ## eval of each fold
RF_cv$aggregate(measure_acc)  ## aggregated scores
RF_cv$prediction()$score(msrs(c("classif.acc", "classif.ce"),
                              average = "micro"))  ## pool predictions across resampling iterations into one Prediction object and then compute the measure on this directly
## confusion matrix
RF_cv$prediction()$confusion

