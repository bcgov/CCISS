#' Create a set of training points with associated
#'   elevation and BGC values.
#'
#' @param bgc_poly an `sf` object (or one cohercible to an `sf` object) of BGC polygons.
#' @param elev a `SpatRaster` or `RasterLayer` of elevation covering the extent of `bgc_poly`.
#' @param gridSize numeric. Distance in m between points.
#' @param crs passed to [sf::st_as_sf()]
#' 
#' @details Points are sampled regularly from a grid with cell size
#'  defined by `gridSize` that covers `bgc_poly`.
#'
#' @return a `data.table` of point coordinates with associated IDs,
#'   elevation and BGCs.
#'   
#' @export
#'
#' @importFrom terra extract vect geom
#' @importFrom sf st_make_grid st_transform
#' @importFrom data.table setDT
makePointCoords <- function(bgc_poly, elev, gridSize = 2000, crs = "EPSG:4326") {
  if (!inherits(bgc_poly, "sf")) {
    bgc_poly <- tryCatch(st_as_sf(bgc_poly), error = function(e) e)
    if (is(bgc_poly, "error")) {
      stop("Can't coherce bgc_poly to an sf object. Please pass an sf object or another cohercible object class.")
    }
  }
  
  if (!is(elev, "SpatRaster")) {
    elev <- tryCatch(rast(elev), error = function(e) e)
    if (is(elev, "error")) {
      stop("Can't coherce elev to a SpatRaster. Please pass SpatRaster or RasterLayer.")
    }
  }
  
  bgc_grid <- st_make_grid(bgc_poly, cellsize = gridSize, what = "centers") |>
    st_transform(crs = st_crs(crs))
  bgc_grid2 <- vect(bgc_grid)
  tmp_elev <- terra::extract(elev, bgc_grid2)
  
  coords <- geom(bgc_grid2, df = TRUE)
  
  setDT(coords)
  coords[, c("part","hole","geom") := NULL]
  coords[, elev := tmp_elev[,2]]
  coords[, id := seq_along(elev)]
  
  return(coords)
}


#' Subset a group of points by spatial extent
#'
#' @param coords `data.table` of point coordinates with columns "x" (longitude)
#'   and "y" (latitude) (and any additional columns), a `SpatVector` or object 
#'   cohersible to `SpatVector`.
#' @param cropExt `SpatExtent` or `SpatVector` (whose extent will be used)
#'   to subset the data to. Defaults to the `SpatExtent` of an area in Southern BC.
#' @param crs passed to [terra::vect()] to coerce coords to a `SpatVector` 
#'   if it is not one already.
#'
#' @return a subset `coords` object.
#' @export
#'
#' @importFrom terra vect crop
#' @importFrom data.table as.data.table
subsetByExtent <- function(coords, cropExt = ext(c(-123, -118, 49, 52)), crs = "EPSG:4326") {
  
  isSpatial <- FALSE
  
  if (!inherits(cropExt, c("SpatVector", "SpatExtent"))) {
    stop("cropExt must be a SpatVector or a SpatExtent")
  } 
  
  if (is(coords, "data.table")) {
    coords <- tryCatch(as.data.table(coords), error = function(e) e)
    if (is(coords, "error")) {
      stop("Can't coherce 'coords' to a data.table.",
           "  \nPlease pass a data.table, a SpatVector, or object to data.table or SpatVector")
    }
    coords_poly <- suppressWarnings(vect(coords, geom = c("x", "y"), crs = crs))
  } else {
    isSpatial <- TRUE
    if (is(coords, "SpatVector")) {
      coords_poly <- coords
    } else {
      coords_poly <- tryCatch(vect(coords), error = function(e) e)
      if (is(coords_poly, "error")) {
        stop("Can't coherce 'coords' to a SpatVector.",
             "  \nPlease pass a data.table, a SpatVector, or object to data.table or SpatVector")
      }
    }
  }

  coords_out <- crop(coords_poly, cropExt)
  
  if (!nrow(coords_out)) {
    warning("No points left after cropping to cropExt. Do projections match?")
  }
  
  if (!isSpatial) {
    coords_out <- as.data.table(coords_out, geom = "XY")
    coords_out <- coords_out[, .SD, .SDcols = names(coords)]  ## re-order
  }
  
  return(coords_out)
}



#' Make generate extents for gaps used in hold-out samples
#'   for model cross-validation.
#'
#' @param studyarea `SpatExtent` of study area where gaps should be generated.
#'   Defaults to an area in Southern BC. 
#' @param ngaps integer. Number of gaps wanted..
#'
#' @return a `list` of extents of length `ngaps`.
#' 
#' @export
#'
#' @importFrom terra ext
makeGapExtents <- function(studyarea = ext(c(-123, -118, 49, 52)), ngaps = 5L) {
  centre <- c(mean(studyarea[1:2]), mean(studyarea[3:4]))
  range <- c(diff(studyarea[1:2]), diff(studyarea[3:4]))
  gapcentre <- matrix(c(-1,1,1,1,1,-1,-1,-1), 4, byrow=T)
  gapextents <- ext(c(centre[1]+c(-1,1)/8*range[1], centre[2]+c(-1,1)/8*range[2]))
  
  ngaps <- ngaps - 1  ## we have one already
  extragaps <- lapply(seq_len(ngaps), function(i) {
    gap <- ext(c(centre[1]+sort(gapcentre[i,1]*c(1,3)/8)*range[1], centre[2]+sort(gapcentre[i,2]*c(1,3)/8)*range[2]))
    gap
  })

  return(append(gapextents, extragaps))
}

#' Prepares coordinates and obtains climate normals
#'  using `climr_downscale`
#'
#' @param coords a `data.table` with point coordinates ("x" = longitude,
#'   "y" = latitude), elevation ("elev") and point IDs ("id").
#' @param bgcs a polygon `sf` object of biogeoclimatic zones to train the model. 
#' @param ... further arguments passed to [climr_downscale()].
#' 
#' @seealso [climr_downscale()]
#'
#' @return climate normals as a `data.table`
#' @export
#'  
#' @importFrom sf st_as_sf  st_join
#' @importFrom data.table data.table
getClimate <- function(coords, bgcs, crs = "EPSG:4326", ...) {
  dots <- list(...)
  
  if (any(!c("lon", "lat", "elev", "id") %in% names(coords))) {
    stop("coords must contain columns 'lon', 'lat', 'elev' and 'id'")
  }
  
  coords_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = crs)
  coords_sf$elev <- NULL
  coords_sf <- st_transform(coords_sf, 3005)
  
  coords_bgc <- st_join(coords_sf, bgcs)
  coords_bgc <- data.table(coords_bgc[,c("id","BGC")])
  coords_bgc[, geometry := NULL]
  coords_bgc <- coords_bgc[!is.na(BGC),]
  
  coords <- as.data.frame(coords)
  
  args <- append(list(coords = coords, coords_bgc = coords_bgc), dots)
  out <- do.call(.getClimVars, args)

  return(out)
}


#' Wrapper for `climr_downscale`
#'
#' @inheritParams getClimate 
#'
#' @return climate normals as a `data.table`.
#' 
#' @importFrom climr climr_downscale
#' @importFrom data.table setDT
.getClimVars <- function(coords, coords_bgc, ...) {
  clim_vars <- climr_downscale(coords, ...) |>
    Cache()
  setDT(clim_vars)
  clim_vars <- clim_vars[!is.nan(PPT05),] ##lots of points in the ocean
  clim_vars[coords_bgc, BGC := i.BGC, on = "id"]
  clim_vars <- clim_vars[!is.na(BGC), ]
  clim_vars[, PERIOD := NULL]
  
  return(clim_vars)
}

#' Add extra climate variables
#'
#' @param dat a `data.table` with columns:
#'    * PPT05, PPT06, PPT07, PPT08, PPT09 (May, ..., September precip), 
#'    PPT_at (autumn precip), PPT_wt (winter precip)
#'    * CMD07 (July climate moisture deficit), CMD (annual CMD)
#'    * DD_0_at (autumn degree-days below 0 deg), DD_0_wt (winter degree-days below 0 deg)
#'    
#' @details This function calculates more climate variables derived from those
#'   output by `climr_downscale`. Presently it adds the following:
#'   *   May to June precip.: \eqn{PPT_MJ = PPT05 + PPT06}
#'   *   July to September precip \eqn{PPT_JAS = PPT07 + PPT08 + PPT09}
#'   *   Precipitation during vegetation dormant period: \eqn{PPT.dormant = PPT_at + PPT_wt}
#'   *   CMD deficit \eqn{CMD.def = 500 - PPT.dormant} (bounded to 0)
#'   *   \eqn{CMDMax = CMD07}
#'   *   \eqn{CMD.total = CMD.def + CMD}
#'   *   \eqn{DD_delayed = ((DD_0_at + DD_0_wt)*0.0238) - 1.8386} bounded to 0)
#'  
#' @return `dat` with all of the above added climate variables.
#' @export
#'
#' @examples
addVars <- function(dat) {
  dat[, PPT_MJ := PPT05 + PPT06]
  dat[, PPT_JAS := PPT07 + PPT08 + PPT09]
  dat[, PPT.dormant := PPT_at + PPT_wt]
  dat[, CMD.def := pmax(0, 500 - PPT.dormant)]
  dat[, CMDMax := CMD07]   ## TODO: THIS IS NOT NECESSARILY CMD MAX
  dat[, CMD.total := CMD.def + CMD]
  dat[, DD_delayed := pmax(0, ((DD_0_at + DD_0_wt)*0.0238) - 1.8386)]
}

#' Log-transform climate variables
#'
#' @param dat a `data.table` with columns of climate variables corresponding to 
#'   the selected climate `elements`.
#' @param elements character. Climate elements to search for in `dat`.
#' @param base numeric. Log base.
#' @param add.vars logical. If `TRUE`, the new logged variables are added to `dat` 
#'   (TRUE). Otherwise, original column values are replaced with the logs (FALSE).
#' @param zero_adjust logical. If `TRUE` adjusts zeroes in raw data as:
#'   \eqn{base^{\log_base{x_min} - 1}.
#'   where \eqn{x_min} is the minimum non-zero, non-NA value.
#'
#' @details
#'   All column names that partially match strings in `elements` will be 
#'   log-transformed.
#' 
#' @return
#' @export
#'
#' @examples
logVars <- function(dat,
                    elements = c("AHM", "DD", "Eref", "FFP", "NFFD", "PAS", "PPT", "SHM", "CMI"),
                    base = exp(1),
                    add.fields = FALSE,
                    zero_adjust = FALSE) {
  
  dat <- copy(dat)
  
  # Fields to operate on (generally these should be ratio (zero-limited) variable)
  logFields <- grep(paste(elements, collapse = "|"), names(dat), value = TRUE)
  dat.log <- dat[, .SD, .SDcols = logFields]
  
  # If specified by the user, give zero values a positive value that is one order of magnitude less than the minimum positive value
  if (zero_adjust) {
    dat.log[, lapply(.SD, function(x) {
      x[x <= 0] <- base^(log(min(x[x > 0], na.rm = TRUE), base = base) - 1)
      return(x)
    })]
  }
  
  # Perform log transformation
  dat.log <- dat.log[, lapply(.SD, function(x) log(x, base = base))]
  
  # Add 
  if(add.fields){
    setnames(dat.log, logFields, paste0(logFields, "_log"))
    dat <- cbind(dat, dat.log)
  } else {
    dat[, (logFields) := Map(x =.SD, xname = logFields, f = function(x, xname) {
      x <- dat.log[[xname]]
      return(x)
    }), .SDcols = logFields]
  }
  return(dat)
}


#' Remove outliers from data
#'
#' @param dat a `data.table` target data to "clean"
#' @param alpha numeric. The alpha value used to determine the cutoff for outliers.
#' @param IDvars character. Names of columns from which outliers should be excluded.
#' 
#' @details TODO. Parallelizes computations internally with `foreach`.
#'
#' @return
#' @seealso [foreach::foreach()]
#' 
#' @importFrom foreach foreach %do%
#' @importFrom stats mahalanobis qchisq cov
#' 
#' @export
removeOutlier <- function(dat, alpha, vars){
  out <- foreach(curr = unique(as.character(dat$BGC)), .combine = rbind) %do% {
    temp <- dat[dat$BGC == curr,]
    md <- tryCatch(mahalanobis(temp[, vars],
                               center = colMeans(temp[, vars]),
                               cov = cov(temp[, vars])), error = function(e) e)
    if (!inherits(md,"error")){
      ctf <- qchisq(1-alpha, df = ncol(temp)-1)
      outl <- which(md > ctf)
      message(paste("Removing", length(outl), "outliers from", curr, "; "), sep = " ")
      if (length(outl) > 0){
        temp <- temp[-outl,]
      }
    }
    temp
  }
  return(out)
}

#' Remove BGCs with low sample sizes
#'
#' @param dat a `data.table` with a "BGC" column and number of rows
#'   being the sampled points.
#' @param cutoff minimum number of points necessary to retain a 
#'   given BGC.
#'
#' @return `dat` without the excluded BGCs and their samples.
#' @importFrom data.table as.data.table
#' @export
rmLowSampleBGCs <- function(dat, cutoff = 30) {
  dat <- as.data.table(dat)
  BGC_Nums <- dat[,.(Num = .N), by = BGC]
  BGC_good <- dat[!BGC %in% BGC_Nums[Num < cutoff, BGC],] 
  return(BGC_good)
}