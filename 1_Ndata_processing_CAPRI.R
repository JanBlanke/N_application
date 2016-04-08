########################################################
### Name: 1_Ndata_processing_script.R                ###
### Author: Jan Blanke                               ###
### Description: Process nitrogen data for LPJ-GUESS ###
########################################################

library(raster)
library(rgdal)
library(gdata)
library(rgdal)
library(simecol)
library(rgeos)
eu.countries.mask <- readOGR("/home/jan/GIS_data", "EU_match_clue")

## Read NUTS2 shapefile
eu.nuts <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/", "NUTS_ETRS_1989_LAEA")
eu.nuts <- spTransform(eu.nuts, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
eu.nuts <- eu.nuts[!is.na(eu.nuts$CAPRI_NUTS), ]

## Gridlist for extracting
grid <- read.table("/home/jan/Dropbox/30arcmin_europe_masked_gridlist.txt", header=F, sep="")
for (i in 1:2) grid[, i] <- grid[, i] + 0.25 
grid <- grid[order(grid[, 1], grid[, 2]), ]
coordinates(grid) <- ~ V1 + V2
proj4string(grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# delete two cells
grid <- grid[-c(1183, 2171)]

## Calculate nearest point matrix
d <- gDistance(grid, spgeom2=NULL, byid=T)
min.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])
dist.df <- cbind(as.data.frame(grid), as.data.frame(grid)[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2]))
colnames(dist.df) <- c(colnames(as.data.frame(grid)), 'n.lat', 'n.long', 'distance')

## Read nitrogen data from Julia
setwd("/home/jan/Dropbox/Paper_2_nitro/Data/CAPRI")
list.files()

## Read Elliot Nitrogen data
sp.nappl <- readRDS("/home/jan/Dropbox/Paper_2_nitro/Data/Elliot/Elliot_spdf.rds")


### 1 ### Loop through years for HISTORIC period ### -----------------------------------------------------------------------------------

# read CAPRI data
#nappl <- read.xls("ToClue_no_header.xlsx", sheet=3) # historic data with from 1990 to 2007
#save(nappl, file = "ToClue_no_header_sheet3.Rdata")
load("ToClue_no_header_sheet3.Rdata")

colnames(nappl)[1] <- c("CAPRI_NUTS")

pb   <- txtProgressBar(min(1990), max(2007), style=3)

wheat.list <- list()
maiz.list <- list()

for (i in 1990:2008){ # loop through historic years, fill wheat and maize lists
  
  nappl.sub <- subset(nappl, nappl$YEAR == i)  
  # join NUTS2 polygons and nitrogen data
  merged <- merge(eu.nuts, nappl.sub, by="CAPRI_NUTS")
  
  swheat.poly <- merged[, "SWHE"]
  maiz.poly <- merged[, "MAIZ"] 
     
  ###################################################
  ## Extract gridlist from CAPRI shapefile for WHEAT!
  over.wheat.capri <- over(grid, as(swheat.poly, "SpatialPolygons"))
  over.wheat.elliot <- over(grid, sp.nappl)[, 2]
    
  ## NAs from gridlist not matching polygon
  na.idx <- as.vector(which(is.na(over.wheat.capri))) # rows from gridlist that are NA
  over.wheat.capri[na.idx] <- 1 # temporarily set to 1 since NA is not allowed
  wheat.df <- cbind(as.data.frame(grid), swheat.poly[over.wheat.capri, ]) # return N values that correspond to each point!!!
  
  for (k in na.idx) {
    n.coords <- dist.df[k, 3:4]
    wheat.df[k, 3] <- subset(wheat.df, V1 == n.coords[[1]] & V2 == n.coords[[2]])[3] #take value from nearest point
  } 
  
  ## NAs N dataset
  na.idx.2 <- which(is.na(wheat.df[, 3])) # rows from df/gridlist which are now NA after returning CAPRI value
  elliot.vals <- over.wheat.elliot[na.idx.2] * 10000
  elliot.vals[elliot.vals > 160] <- 100 # Set Sweden to 100 since the Elliot data for Sweden is unrealistic
  wheat.df[na.idx.2, 3] <- elliot.vals
    
  wheat.df <- na.omit(wheat.df)
  
  rownames(wheat.df) <- 1:nrow(wheat.df)
  colnames(wheat.df)[1:2] <- c("x", "y")
  wheat.df$Year <- paste(i)
  if (i == 1990) j <- 1
  wheat.list[[j]] <- wheat.df # fill list
   
  ###################################################
  ## Extract gridlist from CAPRI shapefile for MAIZE!
  
  over.maiz.capri <- over(grid, as(maiz.poly,"SpatialPolygons"))
  over.maiz.elliot <- over(grid, sp.nappl)[, 1]    
  
  na.idx <- as.vector(which(is.na(over.maiz.capri))) # rows from gridlist that are NA
  over.maiz.capri[na.idx] <- 1 # temporarily set to 1 since NA is not allowed
  maiz.df <- cbind(as.data.frame(grid), maiz.poly[over.maiz.capri, ])
  
  for (k in na.idx) {
    n.coords <- dist.df[k, 3:4]
    maiz.df[k, 3] <- subset(maiz.df, V1 == n.coords[[1]] & V2 == n.coords[[2]])[3] #take value from nearest point
  } 
  
  na.idx.2 <- which(is.na(maiz.df[, 3])) # rows from df/gridlist which are now NA after returning CAPRI value
  
  maiz.df[na.idx.2, 3] <- over.maiz.elliot[na.idx.2] * 10000
  
  maiz.df <- na.omit(maiz.df)
  rownames(maiz.df) <- 1:nrow(maiz.df) 
  maiz.df$Year <- paste(i)
  maiz.list[[j]] <- maiz.df
  j <- j + 1  
  
  ## Temporary: make and write raster 
  #sp.wheat <- SpatialPixelsDataFrame(wheat.df[, 1:2], as.data.frame(wheat.df[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  #sp.maize <- SpatialPixelsDataFrame(maiz.df[, 1:2], as.data.frame(maiz.df[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
  #r.wheat <- raster(sp.wheat)  
  #r.maize <- raster(sp.maize)  
  
  #writeRaster(r.wheat, paste("/home/jan/Dropbox/Nitrogen_maps/CAPRI/hist_data_processed/", "r.wheat.", i, ".tif", sep=""))
  #writeRaster(r.maize, paste("/home/jan/Dropbox/Nitrogen_maps/CAPRI/hist_data_processed/", "r.maize.", i, ".tif", sep=""))
  
  setTxtProgressBar(pb, i)  

}

##################################################################################
## match dimensions
ref.coords <- readRDS("/home/jan/Dropbox/Nitrogen_maps/ref_coords_wheat_maiz.rds")

# WHEAT
for (i in 1:19) print(dim(wheat.list[[i]]))
for (i in 1:19){
  temp <- duplicated(rbind(ref.coords, wheat.list[[i]][, 1:2]))
  temp <- temp[2200:length(temp)]
  idx <- which(temp == FALSE)
  if (length(idx) >= 1) wheat.list[[i]] <- wheat.list[[i]][- idx, ]  
}

# MAIZ
for (i in 1:19) print(dim(maiz.list[[i]]))
for (i in 1:19){
  temp <- duplicated(rbind(ref.coords, maiz.list[[i]][, 1:2]))
  temp <- temp[2200:length(temp)]
  idx <- which(temp == FALSE)
  if (length(idx) >= 1) maiz.list[[i]] <- maiz.list[[i]][- idx, ]  
}

# combine data
library(plyr)
wheat.df.hist <- ldply(wheat.list)
maiz.df.hist <- ldply(maiz.list)

crops.combined.hist <- cbind(wheat.df.hist[, c(1:2, 4, 3)], maiz.df.hist[, 3])
colnames(crops.combined.hist) <- c("lon", "lat", "Year", "SWHE", "MAIZE")
head(crops.combined.hist)
tail(crops.combined.hist)

saveRDS(crops.combined.hist, "/home/jan/Dropbox/Nitrogen_maps/CAPRI/hist_data_processed/CAPRI_combined_hist.rds")
write.table(crops.combined.hist, "/home/jan/Dropbox/Nitrogen_maps/CAPRI/hist_data_processed/CAPRI_combined_hist.txt", quote=F, row.names=F,  sep = "\t")


### 2 ### Loop through years for SCENARIO period  ### -----------------------------------------------------------------------------------
crops.combined.hist <- readRDS("/home/jan/Dropbox/Nitrogen_maps/CAPRI/hist_data_processed/CAPRI_combined_hist.rds")

### Currently: B1

nappl <- read.xls("ToClueB1_no_header.xls", sheet=2) # here: A2
#save(nappl, file = "ToClueB1_no_header_sheet2.Rdata")
#load("ToClueA2_no_header_sheet2.Rdata")

colnames(nappl)[1:2] <- c("CAPRI_NUTS", "Year")

wheat.list <- list()
maiz.list <- list()

for (i in unique(nappl$Year)){
  
  nappl.sub <- subset(nappl, nappl$Year == i)
  
  # join NUTS2 polygons and nitrogen data
  merged <- merge(eu.nuts, nappl.sub, by="CAPRI_NUTS")
  swheat.poly <- merged[, "SWHE"]
  maiz.poly <- merged[, "MAIZ"] 
   
  ## Extract gridlist from CAPRI shapefile for WHEAT!
  over.wheat.capri <- over(grid, as(swheat.poly,"SpatialPolygons"))
  over.wheat.elliot <- over(grid, sp.nappl)[, 2]    
  
  na.idx <- as.vector(which(is.na(over.wheat.capri))) # rows from gridlist that are NA
  over.wheat.capri[na.idx] <- 1 # temporarily set to 1 since NA is not allowed
  wheat.df <- cbind(as.data.frame(grid), swheat.poly[over.wheat.capri, ])
  
  #for (k in na.idx) wheat.df[k, 3] <- modal(wheat.df[(k - 2):(k + 2), 3], na.rm=T)
  for (k in na.idx) {
    n.coords <- dist.df[k, 3:4]
    wheat.df[k, 3] <- subset(wheat.df, V1 == n.coords[[1]] & V2 == n.coords[[2]])[3] #take value from nearest point
  }   
    
  na.idx.2 <- which(is.na(wheat.df[, 3])) # rows from df/gridlist which are now NA after returning CAPRI value
  elliot.vals <- over.wheat.elliot[na.idx.2] * 10000
  elliot.vals[elliot.vals > 160] <- 100 # Set Sweden to 100 since the Elliot data for Sweden is unrealistic
  wheat.df[na.idx.2, 3] <- elliot.vals
  
  wheat.df <- na.omit(wheat.df)
  rownames(wheat.df) <- 1:nrow(wheat.df)
  colnames(wheat.df)[1:2] <- c("x", "y")
  year.substr <- substr(paste(i), 4,7)
  wheat.df$Year <- year.substr
  if (year.substr == 2020) j <- 1
  wheat.list[[j]] <- wheat.df
  
  ## Extract gridlist from CAPRI shapefile for MAIZ!
  over.maiz.capri <- over(grid, as(maiz.poly,"SpatialPolygons"))
  over.maiz.elliot <- over(grid, sp.nappl)[, 1]    
  na.idx <- as.vector(which(is.na(over.maiz.capri))) # rows from gridlist that are NA
  over.maiz.capri[na.idx] <- 1 # temporarily set to 1 since NA is not allowed
  maiz.df <- cbind(as.data.frame(grid), maiz.poly[over.maiz.capri, ])
  
  #for (k in na.idx) maiz.df[k, 3] <- modal(maiz.df[(k - 2):(k + 2), 3], na.rm=T)
  for (k in na.idx) {
    n.coords <- dist.df[k, 3:4]
    maiz.df[k, 3] <- subset(maiz.df, V1 == n.coords[[1]] & V2 == n.coords[[2]])[3] #take value from nearest point
  }
  
  na.idx.2 <- which(is.na(maiz.df[, 3])) # rows from df/gridlist which are now NA after returning CAPRI value
  maiz.df[na.idx.2, 3] <- over.maiz.elliot[na.idx.2] * 10000
  maiz.df <- na.omit(maiz.df)
  rownames(maiz.df) <- 1:nrow(maiz.df) 
  maiz.df$Year <- year.substr
  maiz.list[[j]] <- maiz.df
  j <- j + 1 
  
  ## Temporary: make and write raster 
  #sp.wheat <- SpatialPixelsDataFrame(wheat.df[, 1:2], as.data.frame(wheat.df[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  #sp.maize <- SpatialPixelsDataFrame(maiz.df[, 1:2], as.data.frame(maiz.df[, 3]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  #r.wheat <- raster(sp.wheat)  
  #r.maize <- raster(sp.maize)  
  
  #writeRaster(r.wheat, paste("/home/jan/Dropbox/Nitrogen_maps/CAPRI/b1_data_processed/", "r.wheat.", i, ".tif", sep=""), overwrite=T)
  #writeRaster(r.maize, paste("/home/jan/Dropbox/Nitrogen_maps/CAPRI/b1_data_processed/", "r.maize.", i, ".tif", sep=""), overwrite=T)
    
}

## match dimensions
ref.coords <- readRDS("/home/jan/Dropbox/Nitrogen_maps/ref_coords_wheat_maiz.rds")

# WHEAT
for (i in 1:3) print(dim(wheat.list[[i]]))
for (i in 1:3){
  temp <- duplicated(rbind(ref.coords, wheat.list[[i]][, 1:2]))
  temp <- temp[2200:length(temp)]
  idx <- which(temp == FALSE)
  if (length(idx) >= 1) wheat.list[[i]] <- wheat.list[[i]][- idx, ]  
}

# MAIZ
for (i in 1:3) print(dim(maiz.list[[i]]))
for (i in 1:3){
  temp <- duplicated(rbind(ref.coords, maiz.list[[i]][, 1:2]))
  temp <- temp[2200:length(temp)]
  idx <- which(temp == FALSE)
  if (length(idx) >= 1) maiz.list[[i]] <- maiz.list[[i]][- idx, ]  
}

# combine data
library(plyr)
wheat.df.scen <- rbind.fill(wheat.list)
maiz.df.scen <- rbind.fill(maiz.list)

crops.combined.scen <- cbind(wheat.df.scen[, c(1:2, 4, 3)], maiz.df.scen[, 3])
colnames(crops.combined.scen) <- c("lon", "lat", "Year", "SWHE", "MAIZ")
head(crops.combined.scen)
tail(crops.combined.scen)

saveRDS(crops.combined.scen, "/home/jan/Dropbox/Nitrogen_maps/CAPRI/b1_data_processed/CAPRI_combined_scen.rds")
write.table(crops.combined.scen, "/home/jan/Dropbox/Nitrogen_maps/CAPRI/b1_data_processed/CAPRI_combined_scen.txt", quote=F, row.names=F,  sep = "\t")


### 3 ### Combine, reformat and interpolate HIST and SCEN  ### ----------------------------------------------------------------------------------------------------
# use crops.combined.hist, crops.combined.scen from above

# HIST
crops.hist.rf <- crops.combined.hist[, 3:5]
temp.list <- list(); j <- 1
for (i in 1990:2008) {
  temp.list [[j]] <- subset(crops.hist.rf, crops.hist.rf$Year == i)
  j <- j + 1
}
crops.hist.rf <- do.call(cbind, temp.list)

crops.hist.rf <- crops.hist.rf[, - which(colnames(crops.hist.rf) == "Year")]
crops.hist.rf <- data.frame(t(crops.hist.rf))
crops.hist.rf <- cbind(rep(1990:2008, each=2), crops.hist.rf)
names(crops.hist.rf) <- paste("X", 1:ncol(crops.hist.rf), sep="")

# SCEN
crops.scen.rf <- crops.combined.scen[, 3:5]
temp.list <- list(); j <- 1
for (i in c(2020, 2030, 2040)) {
  temp.list [[j]] <- subset(crops.scen.rf, crops.scen.rf$Year == i)
  j <- j + 1
}
crops.scen.rf <- do.call(cbind, temp.list)

crops.scen.rf <- crops.scen.rf[, - which(colnames(crops.scen.rf) == "Year")]
crops.scen.rf <- data.frame(t(crops.scen.rf))
crops.scen.rf <- cbind(rep(c(2020, 2030, 2040), each=2), crops.scen.rf)
names(crops.scen.rf) <- paste("X", 1:ncol(crops.scen.rf), sep="")

# combine and interpolate
crops.all.rf <- rbind(crops.hist.rf, crops.scen.rf)
crops.all.rf[, 1:6]

wheat.rf <- crops.all.rf[seq(1, nrow(crops.all.rf), by=2), ]
wheat.rf[, 1:6]
maiz.rf <- crops.all.rf[seq(2, nrow(crops.all.rf), by=2), ]
maiz.rf[, 1:6]

wheat.interp <- approxTime(wheat.rf, 1990:2040, rule = 2)
wheat.interp[, 1:6]

maiz.interp <- approxTime(maiz.rf, 1990:2040, rule = 2)
maiz.interp[, 1:6]

#plot(maiz.interp[, 100], type="l")
wheat.t <- data.frame(t(wheat.interp[, - 1]))
maiz.t <- data.frame(t(maiz.interp[, - 1]))

# WHEAT
temp.list <- list()
for (i in 1:ncol(wheat.t)) temp.list [[i]] <- as.data.frame(wheat.t[, i])
wheat.1col <- rbind.fill(temp.list)

# MAIZ
temp.list <- list()
for (i in 1:ncol(maiz.t)) temp.list [[i]] <- as.data.frame(maiz.t[, i])
maiz.1col <- rbind.fill(temp.list)

# Combine again
coords.rep <- crops.combined.hist[1:2199, 1:2]
coords.rep <- coords.rep[rep(seq_len(nrow(coords.rep)), 51), ]

all.interp <- cbind(coords.rep, rep(1990:2040, each=2199), wheat.1col, maiz.1col)

colnames(all.interp) <- c("lon", "lat", "Year", "SWHE", "MAIZ")
head(all.interp)
tail(all.interp)


######################################
### Combine hist and scen and finalize
# TeCo  TeWWs	TeWWw	TeCoi	TeWWsi	TeWWwi
n.all <- cbind(all.interp, all.interp[, c(4, 5, 4, 4)])
colnames(n.all) <- c("lon", "lat", "year", "TeSW", "TeCo", "TeSWi", "TeCoi", "TeWW", "TeWWi")
n.all[, 4:9] <- n.all[, 4:9] / 10000

#reorder for read-in
n.all[, 1:2] <- n.all[, 1:2] - 0.25
n.all <- n.all[order(n.all[, 1], n.all[, 2]), ]

n.all.m <- within(n.all, {
  lon <- sprintf("%.2f", lon)
  lat <- sprintf("%.2f", lat)
  
  TeSW <- sprintf("%f", TeSW)
  TeCo  <- sprintf("%f", TeCo)
  TeSWi  <- sprintf("%f", TeSWi)
  TeCoi  <- sprintf("%f", TeCoi)
  TeWW  <- sprintf("%f", TeWW)
  TeWWi  <- sprintf("%f", TeWWi)
}) 

#write final file
write.table(n.all.m, "/home/jan/Dropbox/Nitrogen_maps/CAPRI/hist_b1_processed/CAPRI_hist_b1_interp.dat", quote=F, row.names=F,  sep = "\t")




