
## Points from a 8km hex grid of western north america are generated in R
## and submitted to ClimateNA to extract annual and seasonal variables for 
## the historic normal period (1961-90) and an ensemble future climate scenario 
## (rcp45 2040-2070). These data sets are combined. Several additional climate 
## variables are generated including several monthly climate sums for precipitation
## and growing degree days. All other monthly variables are removed. A winter rechange 
## moisture deficit was calculated and summed with the climatic moisture deficit to account 
## for regions that begin the growing season in soil moisture deficit.

region <- "RCB"
Y <- fread(paste0("./inputs/", region, "_400m_HexPts_Normal_1991_2020MSY.csv"))
Y <- addVars(Y)

### load ranger  model
load("./BGC_models/WNAv12_Subzone_11_Var_ranger_5Apr22.Rdata")
#load("./BGC_models/WNAv12_Subzone_11_Var_tidy_5Apr22.Rdata")
BGCmodel <- BGCmodel2
model_vars <- as.data.frame(BGCmodel$variable.importance) %>% tibble::rownames_to_column()

covcount <- nrow(model_vars)
idLocs <- Y[,.(ID1,Latitude,Longitude)] 
Y <- Y %>% dplyr::select(ID1, model_vars [,1])
Y1 <- data.table(Y)
Y2 <- Y1[,lapply(.SD,mean),
         by = .(ID1),]
#--------unremark for 1991-2019 period
X2 <- Y2 ## use this line if wanting to run the modern climate change period
timeperiod = "1991-2020"

###Predict for ranger model

grid.pred <- predict(BGCmodel, data = X2[,-c(1)])
BGC <- as.data.frame(grid.pred$predictions)%>% tibble::rownames_to_column() %>% dplyr::rename("BGC" = "grid.pred$predictions")
X1.pred <- cbind(X2, BGC) %>% select(ID1, BGC)
X1.pred$BGC <-  fct_explicit_na(X1.pred$BGC , na_level = "(None)")
X1.pred$ID1 <- as.character(X1.pred$ID1)

# Attribute hex grid with subzone/variant call


##############link predicted Zones to Polygons and write shape file

hexpoly <- st_read(dsn = paste0("./hexpolys/", region, "_bgc_hex400.gpkg"))#, layer = "USA_bgc_hex_800m")
hexpoly$hex_id <- as.character(hexpoly$hex_id)
hexZone <- left_join(hexpoly, X1.pred, by = c("hex_id" = "ID1"))%>% st_transform(3005) %>% st_cast()
temp <- hexZone %>% select(BGC, geom)
temp2 <- st_zm(temp, drop=T, what='ZM') 
##unremark for undissolved layer
#st_write(temp2, dsn = paste0("./outputs/", region, "_", "SubZoneMap_hex400_undissolved.gpkg"), driver = "GPKG", delete_dsn = TRUE)

## unremark for Dissolved
##hexZone <- st_read(dsn = "./outputs/WA_bgc_hex8000_ungrouped.gpkg")#, layer = "USA_bgc_hex_800m") ## need to read it back in to ensure all type Polygon is consistent
temp3 <- hexZone
temp3$BGC <- droplevels(temp3$BGC)
temp3 <-  st_as_sf(temp3)#
st_precision(temp3) <- .5
temp3$BGC <- forcats::fct_explicit_na(temp3$BGC,na_level = "(None)")
temp3 <- temp3[,c("BGC","geom")]
t2 <- aggregate(temp3[,-1], by = list(temp3$BGC), do_union = T, FUN = mean) %>% dplyr::rename(BGC = Group.1)
t2 <- st_zm(t2, drop=T, what='ZM') %>% st_transform(3005) %>% st_cast()

t2 <- t2 %>% st_buffer(0) ## randomly fixes geometry problems
#mapView(t2)
BGC_area <- t2 %>%
  mutate(Area = st_area(.)) %>% mutate (Area = Area/1000000) %>%
  mutate(ID = seq_along(BGC)) %>% dplyr::select(BGC, Area) %>% st_set_geometry(NULL)
#write.csv(BGC_area, paste0("./outputs/", region, "_", timeperiod, "_", covcount, "_BGC_area_predicted_mahbgcreduced.csv"))
st_write(t2, dsn = paste0("./outputs/", covcount, "_", "Subzone_Map_hex400_dissolved_reduced_ranger_6Apr2022.gpkg"), layer = paste0(region, "_", timeperiod, "_", region ,"_", covcount, "_vars_ranger"), driver = "GPKG", delete_layer = TRUE)



### Kiri Variable investigation


#bgc <- st_read("~/CommonTables/BC_BGCv12.gpkg")
bgc <- st_read(paste0(cloud_dir, "BC_BGCv12.gpkg"))
idLocs2 <- st_as_sf(idLocs, coords = c("Longitude","Latitude"), crs = 4326, agr = "constant")
idLocs2 <- st_transform(idLocs2, 3005)
bgc <- st_zm(bgc)
bgc <- st_cast(bgc,"MULTIPOLYGON")
bgcID <- st_join(idLocs2,bgc,left = T)
bgcID <- as.data.table(st_drop_geometry(bgcID))
BGC2 <- as.data.table(BGC)
setnames(BGC2, c("ID1","BGC.pred"))
BGC2[,ID1 := as.integer(ID1)]
BGC2[bgcID, BGC := i.BGC, on = "ID1"]
BGC2[,Change := paste0(BGC,"-",BGC.pred)]
allChanges <- table(BGC2$Change)
sort(allChanges)
allChanges2 <- as.data.frame(allChanges)
##lets try SBSdk -> IDFdk1



origUnit <- "IDFdm1"
changeUnit <-"IDFdm_MT"

changeIDs <- BGC2[BGC == origUnit & BGC.pred == changeUnit,ID1]
datChange <- Y[ID1 %in% changeIDs,]
datChange <- melt(datChange,id.vars = "ID1", variable.name = "ModVar")
datChange[,Unit := "AOIChange"]
datChange[,ID1 := NULL]
modVars <- c("ID2", model_vars$rowname)

datFull <- fread("~/BGC-Climate-Summaries/ClimateBC_Out/WNA_v12_Points_1901-2020MSY.csv")
datFull <- datFull[Year > 1960 & Year < 1991,]
datFull[,`:=`(Latitude = NULL,Longitude = NULL,Elevation = NULL)]
datFull[,Period := fifelse(Year > 1990,"Current","Historic")]
datAvg <- datFull[,lapply(.SD,mean), by = .(ID2,ID1), .SDcols = -"Year"]
datAvg <- addVars(datAvg)
wnaAll <- rbind(datAvg,datAvg2)
datAvg2 <- fread("./BC_BGCv12_Normal6190.csv")

# datFull <- fread("~/BGC-Climate-Summaries/ClimateBC_Out/BC_v12_Points_1901-2020MSY.csv")
# datFull <- datFull[Year > 1960 & Year < 1991,]
#datFull <- datFull[Year > 1959 & ,]
#datFull[,`:=`(Latitude = NULL,Longitude = NULL,Elevation = NULL)]
#datFull[,Period := fifelse(Year > 1990,"Current","Historic")]
#datAvg <- datFull[,lapply(.SD,mean), by = .(ID2,ID1), .SDcols = -"Year"]
#datAvg <- addVars(datAvg)


datAvg <- fread(paste0(cloud_dir, "WNA_BGCv12_Normal6190.csv"))

#datAvg <- fread(paste0(cloud_dir, "BC_BGCv12_Normal6190.csv"))
datAvg <- datAvg[ID2 %chin% c(origUnit,changeUnit),..modVars]

#datAvg <- fread("./BC_BGCv12_Normal6190.csv")
#datAvg <- datAvg[ID2 %chin% c(origUnit,changeUnit),..modVars]

datAvg <- melt(datAvg, id.vars = c("ID2"),variable.name = "ModVar")
setnames(datAvg, old = "ID2", new = "Unit")

datAll <- rbind(datAvg,datChange)
datAll <- datAll[value < 2000 & value > -500,]
datAll[,Unit := factor(Unit, levels = c(origUnit,"AOIChange",changeUnit))]


# temp <- datAll[ModVar == "DD_delayed",]
# temp[,Unit := factor(Unit, levels = c(origUnit,"AOIChange",changeUnit))]

ggplot(datAll,aes(x = Unit, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ ModVar, ncol = 3,scales = "free_y")


# End Kiri



###now cleanup and remove crumbs
library(units)
t3 <- st_cast(t2, "MULTIPOLYGON") %>% st_cast("POLYGON")
t3 <- t3 %>%
  mutate(Area = st_area(.)) %>%
  mutate(ID = seq_along(BGC))
#unique(t3$Area)



size <- 300000
size <- set_units(size, "m^2")
t3$Area <- set_units(size, "m^2")
tSmall <- t3[t3$Area <= size,]
t3$BGC <- as.character(t3$BGC)

require(doParallel)
coreNum <- as.numeric(detectCores()-1)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)

###loop through each polygon < size, determine intersects, and assign to zone with most edge touching
###all the built in functions I found only dealt with holes in the middle of polygons
i = 1
new <- foreach(i = 1:length(tSmall$ID), .combine = rbind, .packages = c("foreach","sf")) %dopar% {
  ID <- tSmall$ID[i]
  nbrs <- st_intersects(tSmall[i,],t3)[[1]]
  nbrs <- nbrs[!nbrs %in% ID]
  if(length(nbrs) == 0){return(NULL)}
  lines <- st_intersection(t3[ID,],t3[nbrs,])
  lines <- st_cast(lines)
  l.len <- st_length(lines)
  names(l.len) <- lines$BGC.1
  zn <- names(l.len)[l.len == max(l.len)][1]
  newDat <- t3[ID,]
  newDat$BGC <- zn
  newDat
}

stopCluster(coreNo)
gc()
temp <- t3[!t3$ID %in% new$ID,]
t3 <- rbind(temp, new) %>%
  mutate(Zone = as.factor(BGC))

###now have to combine crumbs with existing large polygons
temp2 <- t3
st_precision(temp2) <- 0.5
t3 <- temp2 %>%
  group_by(BGC) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

#mapview(t2, zcol = "BGC")
t3 <- st_zm(t3, drop=T, what='ZM')
t3 <- t3 %>% st_buffer (0)
st_write(t3, dsn = paste0("./outputs/", region, "_SubZoneMap_1991_2019_eliminated.gpkg"), driver = "GPKG")




