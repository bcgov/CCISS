## This script must be sourced

##read in provincial outline, tile, and output for climatebc
##Kiri Daust

## create siteno, bgcs, dist_code table for preselected by BEC
hexgrd <- reproducible::prepInputs(url = "//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/HexGrid400m_Sept2021.gpkg",
                                   targetFile = "HexGrid400m_Sept2021.gpkg",
                                   fun = "sf::st_read")

bgcs <- reproducible::prepInputs(url = "//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022.gpkg",
                                 targetFile = "WNA_BGC_v12_5Apr2022.gpkg",
                                 fun = "sf::st_read",
                                 projectTo = st_crs(4326))

bgcs <- st_zm(bgcs)
bgcs <- st_cast(bgcs, "MULTIPOLYGON")

bgcs_hexgrd <- st_join(hexgrd, bgcs)

regions <- st_read("~/CommonTables/ForestRegions.gpkg","ForestDistricts2")
regions <- regions["ORG_UNIT"]
bgcs_hexgrd <- st_join(bgcs_hexgrd, regions)
grd2 <- as.data.table(st_drop_geometry(bgcs_hexgrd))
grd2 <- na.omit(grd2)
setnames(grd2, c("siteno","zone","bgcs","dist_code"))
dbWriteTable(con,"siteidx",grd2, row.names = F)
dbExecute(con,"create index on siteidx(bgcs)")
dbExecute(con,"create index on siteidx(dist_code,bgcs)")
#bgcs <- bgcs[is.na(bgcs$State),c("bgcs")]

##This is what I used to create the new hex grid
dem <- raster("./BigDat/BC_25m_DEM_WGS84.tif")
BC <- st_read(dsn = "./BigDat/BC_Province_Outline_Clean4.gpkg")
BC <- ms_simplify(BC, keep = 0.1,sys = T)
library(mapview)
mapview(BC)
st_write(BC,dsn = "BC_Simplified.gpkg")
st_is_valid(BC)

BC <- st_read("./BigDat/TempBC2.gpkg")##very simple outline
st_crs(BC) <- 3005
bgcs_hexgrd <- st_make_grid(BC,cellsize = 4000, square = F, flat_topped = F) ##make grid
st_write(bgcs_hexgrd, dsn = "Grid4km.gpkg")
bgcs_hexgrd <- st_as_sf(data.frame(ID = 1:length(bgcs_hexgrd),geom = bgcs_hexgrd))
bgcs_hexgrd <- st_join(bgcs_hexgrd,BC)
bgcs_hexgrd <- bgcs_hexgrd[!is.na(bgcs_hexgrd$State),]
st_write(bgcs_hexgrd, dsn = "Grid4km.gpkg", delete_dsn = T)
pts_all <- st_centroid(bgcs_hexgrd)
pts_all <- pts_all["ID"]
pts_all$ID <- seq(1:nrow(pts_all))
colnames(pts_all)[1] <- "siteno"
st_write(pts_all,con,"pts_4km")

BC <- st_read("./BC_Simplified.gpkg") ##actual outline
bcRast <- raster(BC,resolution = c(2000,2000))
centres <- coordinates(bcRast)
centres <- as.data.table(centres)
centres[,RastID := 1:nrow(centres)]
rast_points <- st_as_sf(centres,coords = c("x","y"), crs = 3005)
rast_pointsBC <- st_join(rast_points,BC,left = F)
rast_pointsBC$State <- NULL
colnames(rast_pointsBC)[1] <- "rast_id"
bcRast[rast_pointsBC$RastID] <- 5
values(bcRast) <- 1
writeRaster(bcRast,"BC_Raster.tif", format = "GTiff",overwrite = T)
st_write(rast_pointsBC,"BC_Raster_Centroids.gpkg")

library(rpostgis)
library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user = "postgres", 
                 host = "138.197.168.220",
                 password = "PowerOfBEC", port = 5432, 
                 dbname = "cciss") ### for local use

pgWriteRast(con,name = "bc_raster", raster = bcRast, overwrite = T)
st_write(rast_pointsBC,con,"rast_centroids")

dbExecute(con,"create table pts2km_ids as
          select rast_centroids.rast_id,
          hex_grid.siteno
          from rast_centroids
          inner join hex_grid
          on ST_Intersects(rast_centroids.geometry,hex_grid.geom)")
dbExecute(con,
          "create table pts2km_future as
          select pts2km_ids.rast_id,
          cciss_future12.gcm, cciss_future12.scenario, 
          cciss_future12.futureperiod, cciss_future12.bgcs, 
          cciss_future12.bgc_pred
          from pts2km_ids
          inner join cciss_future12 
          on pts2km_ids.siteno = cciss_future12.siteno")

dbExecute(con,"create index on pts2km_future (futureperiod)")
dbExecute(con, "create index on pts2km_future (gcm,scenario,futureperiod)")
square_grd <- st_make_grid(BC,cellsize = 4000, square = T)
square_grd <- st_as_sf(data.frame(geom = square_grd))
square_grd <- st_join(square_grd,BC,left = F)
square_grd$siteno <- 1:nrow(square_grd)
square_grd$State <- NULL
bgcs_hexgrd <- st_read("~/HexGrid_Tile/BC_BiggerHex.gpkg")

# bgcs_hexgrd$geom <- bgcs_hexgrd$geom + c(-167,-95) ##adjust so matches old centroids
# bgcs_hexgrd$siteno <- 1:nrow(bgcs_hexgrd)
# bgcs_hexgrd <- bgcs_hexgrd["siteno"]
# st_write(bgcs_hexgrd,dsn = "./BC_HexPoly400m.gpkg", driver = "GPKG")

###intersect with old points
Rcpp::sourceCpp("./C_Helper.cpp")
bgcs_hexgrd <- st_read(dsn = "BC_HexPoly400m.gpkg")

oldPoints <- st_read("./RCB_Hex400_Points.gpkg")##old centre points
colnames(oldPoints)[1] <- "OldID"
st_crs(bgcs_hexgrd) <- 3005
temp <- st_intersects(bgcs_hexgrd,oldPoints,sparse = T) ##faster than join
test <- unlist_sgbp(temp) ##c++ function
fwrite(test,"RCB_CrosswalkTable.csv")

crosswalk <- data.table(NewID = bgcs_hexgrd$ID, OldID = test)
crosswalk[OldID == 0, OldID := NA]
pntsNeeded <- crosswalk[is.na(OldID),NewID]

###now this is Will's part
#dem <- raster("./BigDat/BC_25m_DEM_WGS84.tif")
dem <- raster("D:/CommonTables/DEMs/BC_25m_DEM_WGS84.tif")
#BC <- st_read(dsn = "./BigDat/BC_Province_Outline_Clean.gpkg")
BC <- st_read(dsn = "D:/CommonTables/BC_AB_US_Shp/BC_Province_Outline_Clean.gpkg")
BC <- st_buffer(BC, dist = 0)
BC <- ms_simplify(BC, keep = 0.2)
st_layers("~/CommonTables/ForestRegions.gpkg")
regions <- st_read("~/CommonTables/ForestRegions.gpkg","ForestRegions_clipped")
rcb <- regions[regions$ORG_UNIT == "RCB","ORG_UNIT"]
rcb <- ms_simplify(rcb,keep = 0.05,sys = T)

bb <- st_as_sfc(st_bbox(rcb))
bgcs_hexgrd <- st_make_grid(bb,cellsize = 400, square = F, flat_topped = F)
bgcs_hexgrd <- st_as_sf(data.frame(ID = 1:length(bgcs_hexgrd),geom = bgcs_hexgrd))
ptsAll <- st_centroid(bgcs_hexgrd)
grdPts <- st_sf(ID = seq(length(ptsAll)), geometry = ptsAll)
st_write(ptsAll, dsn = "./RCB_Hex400_Points.gpkg", layer = "HexPts400", driver = "GPKG", overwrite = T, append = F)

ptsAll <- st_centroid(bgcs_hexgrd)
st_write(bgcs_hexgrd, dsn = "./RCB_Hex400_Poly.gpkg", driver = "GPKG", overwrite = T, append = F)

ptsNew <- ptsAll[ptsAll$siteno %in% pntsNeeded,]

grdPts <- st_read(dsn = "BC_HexPoints400m.gpkg")


tiles <- st_make_grid(BC, cellsize = c(250000,vertDist(307000)))
plot(tiles)
plot(BC, add = T)

tilesID <- st_as_sf(data.frame(tID = 1:length(tiles)),geom = tiles)
tilesUse <- st_join(tilesID,BC)
tilesUse <- tilesUse[!is.na(tilesUse$State),]
tilesUse <- tilesUse[,"tID"]
tilesUse <- unique(tilesUse)
library(mapview)
mapview(BC)+
  tilesUse

tilesUse$tID <- 1:nrow(tilesUse)
st_write(tilesUse,dsn = "TileOutlines.gpkg",overwrite = T, append = F)
bgcs <- st_zm(bgcs)
bgcs <- st_cast(bgcs,"MULTIPOLYGON")
testGrd <- st_zm(testGrd)

datOut <- foreach(tile = tilesUse$tID, .combine = rbind) %do% {
  cat("Processing tile",tile,"... \n")
  testTile <- tilesUse[tile,]
  testGrd <- st_intersection(grdPts, testTile)
  if(nrow(testGrd) > 1){
    grdBGC <- st_join(testGrd,bgcs)
    grdBGC <- st_transform(grdBGC, 4326)
    grdBGC$el <- raster::extract(dem, grdBGC)
    out <- cbind(st_drop_geometry(grdBGC),st_coordinates(grdBGC)) %>% as.data.table()
    out <- out[,.(ID1 = ID, ID2 = bgcs, lat = Y, long = X, el)]
    out[,TileNum := tile]
    out
  }else{
    NULL
  }
  
}

datOut <- out
dat <- unique(datOut, by = "ID1")
dat <- dat[!is.na(ID2),]
fwrite(dat, "RCB_ClimBC.csv",eol = "\r\n")
tileID <- dat[,.(ID1,TileNum)]
fwrite(tileID,"TileIDs.csv")

for(i in unique(tileID$TileNum)){
  dat2 <- dat[TileNum == i,]
  dat2[,TileNum := NULL]
  dat2 <- dat2[complete.cases(dat2),]
  fwrite(dat2, paste0("./Output/Tile",i,"_In.csv"), eol = "\r\n")
}

### which IDs have been successfully download by climBC
# ids <- fread("/media/data/ClimateBC_Data/Tile12_In_280 GCMsMSY.csv", select = "ID1")
tileNum <- 13
idDownload <- fread("/media/data/ClimateBC_Data/Tile12_CutTest.csv")
ids <- idDownload$ID1
length(ids[ids == 889155])

##
maxSize <- 200000
tiles <- c(13:22)
for(tile in tiles){
  dat <- fread(paste0("./Output/Tile",tile,"_In.csv"))
  n = nrow(dat)
  brks <- seq(1,n,by = maxSize)
  brks <- c(brks,n)
  i = 0
  for(j in 1:(length(brks)-1)){
    i = i+1
    temp <- dat[brks[j]:brks[j+1],]
    fwrite(temp, file = paste0("./Output/Tile",tile,"_",i,".csv"))
  }
  
}

