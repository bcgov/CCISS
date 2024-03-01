## This script must be sourced

## FIT AND EVALUATE A BGC MODEL

## TODO: need proper documentation and references to sources/scripts that produced input objects

dPath <- unlist(options("reproducible.destinationPath"))

dwnldFromObjSto(prefix = "~/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022",
                bucket = "gmrtde",
                path = dPath)
bgcs <- vect(file.path(dPath, "WNA_BGC_v12_5Apr2022.gpkg"))  ## for poly validation if need be
.gc()

dwnldFromObjSto(prefix = "~/DEM/DEM_NorAm/NA_Elevation/data/northamerica",
                bucket = "gmrtde",
                path = dPath)
elev <- rast(file.path(dPath, "northamerica_elevation_cec_2023.tif"))

cacheExtra <- list(summary(elev), table(bgcs$BGC))
elev <- Cache(postProcessTerra,
              from = elev,
              cropTo = bgcs,
              projectTo = bgcs,
              maskTo = NA,
              .cacheExtra = cacheExtra,
              userTags = "elev", 
              omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
.gc()
              
coords_train <- Cache(makePointCoords,
                bgc_poly = bgcs,
                elev = elev,
                      .cacheExtra = cacheExtra,
                userTags = "coords",
                omitArgs = c("userTags", "bgc_poly", "elev"))
coords_train <- coords_train[!is.na(elev),]
coords_train[, id := seq_along(id)]  ## re-do to have contiguous IDs

vars_needed <- c("DD5", "DD_0_at", "DD_0_wt", "PPT05", "PPT06", "PPT07", "PPT08",
                 "PPT09", "CMD", "PPT_at", "PPT_wt", "CMD07", "SHM", "AHM", "NFFD", "PAS", "CMI")

clim_vars <- Cache(
  getClimate,
  coords = coords_train, 
  bgcs = bgcs,
  which_normal = "normal_composite", 
  return_normal = TRUE, 
  vars = vars_needed, 
  cache = TRUE,
  .cacheExtra = list(summary(bgcs), summary(coords_train)),
  userTags = "clim_vars",
  omitArgs = c("userTags", "coords", "bgcs"))

## subset coords objects to ids with data
coords_train <- coords_train[clim_vars[, .(id)], on = "id", nomatch = 0L]

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
lrn_rf_resp <- lrn("classif.ranger", 
                   predict_type = "response",
                   num.trees = 501,
                   splitrule =  "extratrees",
                   mtry = 4,
                   min.node.size = 2,
                   importance = "permutation",
                   write.forest = TRUE)

## we'll also fit a probability type model for current/normal period predictions
lrn_rf_prob <- lrn("classif.ranger", 
                   predict_type = "probability",
                   num.trees = 501,
                   splitrule =  "extratrees",
                   mtry = 4,
                   min.node.size = 2,
                   importance = "permutation",
                   write.forest = TRUE)

measure_acc <- msrs(c("classif.acc", "classif.ce", "oob_error"))

## fit full models --------------------------------------------
lrn_rf_resp$train(tsk_bgc)
lrn_rf_resp$model

lrn_rf_prob$train(tsk_bgc)
lrn_rf_prob$model
## this plot is useless, but code is kept here for fut reference
# autoplot(lrn_rf_resp$predict(tsk_bgc)) +
#   theme(legend.position = "none")

## evaluate model with CV ------------------------------------
## define and apply cv strategy
folds <- 10
cv_strategy <- rsmp("repeated_spcv_coords", folds = folds, repeats = 1)

future::plan("multisession", 
             workers = ifelse(folds <= future::availableCores(), folds, future::availableCores()))
RF_cv_resp <- mlr3::resample(tsk_bgc, lrn_rf_resp, cv_strategy, store_models = TRUE) |>
  Cache()
future:::ClusterRegistry("stop")

future::plan("multisession", 
             workers = ifelse(folds <= future::availableCores(), folds, future::availableCores()))
RF_cv_prob <- mlr3::resample(tsk_bgc, lrn_rf_prob, cv_strategy, store_models = TRUE) |>
  Cache()
future:::ClusterRegistry("stop")

## eval metrics
## TODO: find a pretty way to export these
RF_cv_prob$score(measure_acc)  ## eval of each fold
RF_cv_prob$aggregate(measure_acc)  ## aggregated scores
RF_cv_prob$prediction()$score(msrs(c("classif.acc", "classif.ce"),
                              average = "micro"))  ## pool predictions across resampling iterations into one Prediction object and then compute the measure on this directly
## confusion matrix
RF_cv_prob$prediction()$confusion

## TODO: save/export models and evaluation metrics