##############################################################
### Name: Fig5_relative_comb_change_linegraph_facetgrid.R  ###
### Author: Jan Blanke                                     ###
### Description: Plot relative change (yield, leaching)    ###
### Source: standalone, allDynamic simulations             ### 
##############################################################

library(rgdal)
library(rasterVis)
library(plyr)

eu.countries.mask <- readOGR("/home/jan/GIS_data", "EU_match_clue")
eu.countries.underlay <- readOGR("/home/jan/GIS_data", "EU_27_underlay")
source("/home/jan/Dropbox/R/functions/colorramps.R")

### Process simulation data for YIELD -------------------------------------------------------------------


## Settings 
#setwd("/media/jan/AddData/Simulations/Paper_2_nitrogen/1_sims_predictions/") 
setwd("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/allDynamic/output")
modality <- "allDynamic" # allDynamic, constN, constCO2, constCLIM, Elliot, Stocker
data.out <- "yield" # yield, nflux, rcp45, rcp85
crops <- c("TeWW", "TeCo") # TeCo, TeWW
years.ts <- c(2008, 2020, 2030, 2040)
years.avg <- 1 # year +- 1 
firstyear <- 1850
##


list.files() # 5 GCMs will be used
file.out <- grep(paste(data.out, ".out", sep=""), list.files(), value=T)

gcm.list <- list()

for (j in file.out){ # read whole table
  
  out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
  
  cft.stack.list <- list()
  
  for (i in years.ts){ # subset years, crops, ...
    
    years <- (i - years.avg):(i + years.avg)
    #modelyears <- years - firstyear + 500
    df.out.sub <- subset(out.df, out.df[3] >= min(years) & out.df[3] <= max(years))
    
    if (length(years) > 1 && grepl("yield", data.out)){
      
      df.out.sub <- ddply(df.out.sub, c(.(Lon), .(Lat)), summarise, mean.teww = mean(TeWW, na.rm = TRUE), mean.tewwi = mean(TeWWi, na.rm = TRUE), mean.tesw = mean(TeSW, na.rm = TRUE), mean.teswi = mean(TeSWi, na.rm = TRUE), mean.teco = mean(TeCo, na.rm = TRUE), mean.tecoi = mean(TeCoi, na.rm = TRUE))
      df.out.sub <- cbind(df.out.sub[, 1:2], rep(mean(years), 2199), df.out.sub[, 3:ncol(df.out.sub)]); colnames(df.out.sub)[3] <- "Year"
      
    } else { # N leaching!
      
      df.out.sub <- ddply(df.out.sub, c(.(Lon), .(Lat)), summarise, mean.leach = mean(leach, na.rm = TRUE))
      df.out.sub <- cbind(df.out.sub[, 1:2], rep(mean(years), 2199), df.out.sub$mean.leach) 
      colnames(df.out.sub)[3:4] <- c("Year", "nleach")
      
    }
    
    sp.df.out <- SpatialPixelsDataFrame(df.out.sub[, 1:2], as.data.frame(df.out.sub[, 4:ncol(df.out.sub)]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    out.stack <- stack(sp.df.out)
    names(out.stack) <- names(sp.df.out)
    
    if (j == file.out[1] && grepl("yield", data.out)){ # load crop fracs in first iteration
      # read crop frac file for weighing
      crop.fracs <- read.table("/media/jan/AddData/Data/Stefan_data/cropland_hurtt_mirca_2000_corrsum.out", head=T)
      crop.fracs[, 1:2] <- crop.fracs[, 1:2] + 0.25
      head(crop.fracs)
      cfracs.sp <- SpatialPixelsDataFrame(crop.fracs[, 1:2], as.data.frame(crop.fracs[, 3:8]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      cfracs.stack <- stack(cfracs.sp)
      frac.mask <- crop(cfracs.stack, out.stack[[1]])
      frac.mask <- mask(frac.mask, out.stack[[1]])
      frac.mask <- frac.mask[[c(2, 5, 1, 4, 3, 6)]]
    }
    
    if (grepl("yield", data.out)) { # for yield
      out.stack <- out.stack * 10 / (0.85 * 2.0 * 0.446) # convert to tonnes per hectar and wet matter
      yield.wt.tew <- ((out.stack[[1]] * frac.mask[[1]]) + (out.stack[[2]] * frac.mask[[2]]) + (out.stack[[3]] * frac.mask[[3]]) + (out.stack[[4]] * frac.mask[[4]])) / (sum(frac.mask[[1:4]]))
      yield.wt.teco <- (out.stack[[5]] * frac.mask[[5]] + out.stack[[6]] * frac.mask[[6]]) / (sum(frac.mask[[5:6]]))
    } else { # for N leaching
      r.nleach <- out.stack
    }  
    
    idx <- which(i == years.ts)
    
    if (grepl("yield", data.out)) { # for yield
      cft.stack.list[[idx]] <- yield.wt.tew # fill head of list with wheat
      cft.stack.list[[idx + length(years.ts)]] <- yield.wt.teco # fill tail with maize
    } else { # for N leaching
      cft.stack.list[[idx]] <- r.nleach 
    } 
  }
  
  all.stack <- stack(cft.stack.list)
  
  if (grepl("yield", data.out)) {
    names(all.stack) <- paste(rep(crops, each=4), ".", years.ts, sep="")
  } else {
    names(all.stack) <- paste(rep(data.out, each=4), ".", years.ts, sep="")
  }
  
  idxy <- which(j == file.out)
  gcm.list[[idxy]] <- all.stack
  
} # run 5 models * 2 rcps


### Calculate YIELD ensemble mean and sd for each time slice. Note: sd not based on area, only on GCMs --------------------------------
if (grepl("yield", data.out)) {
  
  for (i in c("rcp45", "rcp85")){ 
    
    idx.list <- grep(i, file.out, value=F)
    # 5 gcms, 12 layers(crop * year)
    
    # 2008
    s.2008 <- stack(gcm.list[[idx.list[1]]][[1]], gcm.list[[idx.list[2]]][[1]], gcm.list[[idx.list[3]]][[1]], gcm.list[[idx.list[4]]][[1]], gcm.list[[idx.list[5]]][[1]], gcm.list[[idx.list[1]]][[5]], gcm.list[[idx.list[2]]][[5]], gcm.list[[idx.list[3]]][[5]], gcm.list[[idx.list[4]]][[5]], gcm.list[[idx.list[5]]][[5]])
    
    # 2020
    s.2020 <- stack(gcm.list[[idx.list[1]]][[2]], gcm.list[[idx.list[2]]][[2]], gcm.list[[idx.list[3]]][[2]], gcm.list[[idx.list[2]]][[2]], gcm.list[[idx.list[5]]][[2]], gcm.list[[idx.list[1]]][[6]], gcm.list[[idx.list[2]]][[6]], gcm.list[[idx.list[3]]][[6]], gcm.list[[idx.list[4]]][[6]], gcm.list[[idx.list[5]]][[6]])
    
    # 2030
    s.2030 <- stack(gcm.list[[idx.list[1]]][[3]], gcm.list[[idx.list[2]]][[3]], gcm.list[[idx.list[3]]][[3]], gcm.list[[idx.list[4]]][[3]], gcm.list[[idx.list[5]]][[3]], gcm.list[[idx.list[1]]][[7]], gcm.list[[idx.list[2]]][[7]], gcm.list[[idx.list[3]]][[7]], gcm.list[[idx.list[4]]][[7]], gcm.list[[idx.list[5]]][[7]])
    
    # 2040
    s.2040 <- stack(gcm.list[[idx.list[1]]][[4]], gcm.list[[idx.list[2]]][[4]], gcm.list[[idx.list[3]]][[4]], gcm.list[[idx.list[4]]][[4]], gcm.list[[idx.list[5]]][[4]], gcm.list[[idx.list[1]]][[8]], gcm.list[[idx.list[2]]][[8]], gcm.list[[idx.list[3]]][[8]], gcm.list[[idx.list[4]]][[8]], gcm.list[[idx.list[5]]][[8]])
    
    # TODO: calculate sd for each area from all 5 gcms, the same for mean !!!!!!
    setwd("/home/jan/GIS_data/Agroclimatic_zones")
    ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")
    
    ## load zone rasters
    files.r <- list.files(pattern="*.tif$")[1:9]
    
    # use extract
    mean.list.2008 <- mean.list.2020 <- mean.list.2030 <- mean.list.2040 <- list()
    sd.list.2008 <- sd.list.2020 <- sd.list.2030 <- sd.list.2040 <- list()
    
    for (z in files.r){
      ptmp <- rasterToPoints(raster(z), spatial=T)
      zone.extr.2008 <- as.data.frame(extract(s.2008, ptmp))
      zone.extr.2020 <- as.data.frame(extract(s.2020, ptmp))
      zone.extr.2030 <- as.data.frame(extract(s.2030, ptmp))
      zone.extr.2040 <- as.data.frame(extract(s.2040, ptmp))
      
      #mean and sd over area and gcms
      #wheat.2008.mean <- median(zone.extr.2008[, 1:5], na.rm=T) # median is better
      #wheat.2008.sd <- sd(zone.extr.2008[, 1:5], na.rm=T) 
      
      # ONLY mean and sd over gcms
      idx <- which(z == files.r)
      mean.list.2008[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2008[, 1:5]))), median(as.numeric(colwise(median, na.rm=T)(zone.extr.2008[, 6:10]))))
      mean.list.2020[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2020[, 1:5]))), median(as.numeric(colwise(median, na.rm=T)(zone.extr.2020[, 6:10]))))
      mean.list.2030[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2030[, 1:5]))), median(as.numeric(colwise(median, na.rm=T)(zone.extr.2030[, 6:10]))))
      mean.list.2040[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2040[, 1:5]))), median(as.numeric(colwise(median, na.rm=T)(zone.extr.2040[, 6:10]))))
      
      sd.list.2008[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2008[, 1:5]))), sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2008[, 6:10])))) 
      sd.list.2020[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2020[, 1:5]))), sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2020[, 6:10]))))
      sd.list.2030[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2030[, 1:5]))), sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2030[, 6:10]))))
      sd.list.2040[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2040[, 1:5]))), sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2040[, 6:10]))))
    }
    
    all.df.sd <- as.data.frame(cbind(unlist(sd.list.2008)[seq(1, 18, 2)], unlist(sd.list.2008)[seq(2, 18, 2)],
                                     unlist(sd.list.2020)[seq(1, 18, 2)], unlist(sd.list.2020)[seq(2, 18, 2)],
                                     unlist(sd.list.2030)[seq(1, 18, 2)], unlist(sd.list.2030)[seq(2, 18, 2)],
                                     unlist(sd.list.2040)[seq(1, 18, 2)], unlist(sd.list.2040)[seq(2, 18, 2)]))
    
    colnames(all.df.sd) <- paste(rep(crops, times=4), ".", rep(years.ts, each=2), sep="") # each row = region
    rownames(all.df.sd) <- substr(files.r,1,nchar(files.r) - 4)
    
    all.df.mean <- as.data.frame(cbind(unlist(mean.list.2008)[seq(1, 18, 2)], unlist(mean.list.2008)[seq(2, 18, 2)],
                                       unlist(mean.list.2020)[seq(1, 18, 2)], unlist(mean.list.2020)[seq(2, 18, 2)],
                                       unlist(mean.list.2030)[seq(1, 18, 2)], unlist(mean.list.2030)[seq(2, 18, 2)],
                                       unlist(mean.list.2040)[seq(1, 18, 2)], unlist(mean.list.2040)[seq(2, 18, 2)]))
    
    colnames(all.df.mean) <- paste(rep(crops, times=4), ".", rep(years.ts, each=2), sep="") # each row = region
    rownames(all.df.mean) <- substr(files.r,1,nchar(files.r) - 4)
    
    setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/2_predictions/allDynamic/bugfix_branch/")
    
    file.name <- paste("Yield", "_ensemble_mean_", modality, "_" , i, ".rds", sep="") 
    saveRDS(all.df.mean, file=file.name)
    
    file.name <- paste("Yield", "_ensemble_sd_", modality, "_" , i, ".rds", sep="")  
    saveRDS(all.df.sd, file=file.name)
    
  } # loop through c("rcp45", "rcp85")
  
} 




### Process NLEACH simulation data -------------------------------------------------------------------------------------------------------

## Settings 
#setwd("/media/jan/AddData/Simulations/Paper_2_nitrogen/1_sims_predictions/") # N switched off
setwd("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/allDynamic/output")
modality <- "allDynamic" # allDynamic, constN, constCO2, constCLIM, Elliot, Stocker
data.out <- "nflux" # yield, nflux, rcp45, rcp85
#crops <- c("TeWW", "TeCo") # TeCo, TeWW
years.ts <- c(2008, 2020, 2030, 2040)
years.avg <- 1 # year +- 1 
firstyear <- 1850
##


list.files() # 5 GCMs will be used
file.out <- grep(paste(data.out, ".out", sep=""), list.files(), value=T)

gcm.list <- list()

for (j in file.out){ # read whole table
  
  out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
  
  cft.stack.list <- list()
  
  for (i in years.ts){ # subset years, crops, ...
    
    years <- (i - years.avg):(i + years.avg)
    #modelyears <- years - firstyear + 500
    df.out.sub <- subset(out.df, out.df[3] >= min(years) & out.df[3] <= max(years))
    
    df.out.sub <- ddply(df.out.sub, c(.(Lon), .(Lat)), summarise, mean.leach = mean(leach, na.rm = TRUE))
    df.out.sub <- cbind(df.out.sub[, 1:2], rep(mean(years), 2199), df.out.sub$mean.leach) 
    colnames(df.out.sub)[3:4] <- c("Year", "nleach")
    
    sp.df.out <- SpatialPixelsDataFrame(df.out.sub[, 1:2], as.data.frame(df.out.sub[, 4:ncol(df.out.sub)]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    out.stack <- stack(sp.df.out)
    names(out.stack) <- names(sp.df.out)
    
    r.nleach <- out.stack
    
    
    idx <- which(i == years.ts)
    cft.stack.list[[idx]] <- r.nleach 
    
  }
  
  all.stack <- stack(cft.stack.list)
  names(all.stack) <- paste(rep(data.out, each=4), ".", years.ts, sep="")
  
  idxy <- which(j == file.out)
  gcm.list[[idxy]] <- all.stack
  
} # run 5 models and 2 rcps


### Calculate NLEACH ensemble mean and sd for each time slice. Note: sd not based on area, only on GCMs ----------------------------

for (i in c("rcp45", "rcp85")){ 
  
  idx.list <- grep(i, file.out, value=F)
  
  # 2008
  s.2008 <- stack(gcm.list[[idx.list[1]]][[1]], gcm.list[[idx.list[2]]][[1]], gcm.list[[idx.list[3]]][[1]], gcm.list[[idx.list[4]]][[1]], gcm.list[[idx.list[5]]][[1]])
  
  # 2020
  s.2020 <- stack(gcm.list[[idx.list[1]]][[2]], gcm.list[[idx.list[2]]][[2]], gcm.list[[idx.list[3]]][[2]], gcm.list[[idx.list[2]]][[2]], gcm.list[[idx.list[5]]][[2]])
  
  # 2030
  s.2030 <- stack(gcm.list[[idx.list[1]]][[3]], gcm.list[[idx.list[2]]][[3]], gcm.list[[idx.list[3]]][[3]], gcm.list[[idx.list[4]]][[3]], gcm.list[[idx.list[5]]][[3]])
  
  # 2040
  s.2040 <- stack(gcm.list[[idx.list[1]]][[4]], gcm.list[[idx.list[2]]][[4]], gcm.list[[idx.list[3]]][[4]], gcm.list[[idx.list[4]]][[4]], gcm.list[[idx.list[5]]][[4]])
  
  # TODO: calculate sd for each area from all 5 gcms, the same for mean !!!!!!
  setwd("/home/jan/GIS_data/Agroclimatic_zones")
  ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")
  
  ## load zone rasters
  files.r <- list.files(pattern="*.tif$")[1:9]
  
  # use extract
  mean.list.2008 <- mean.list.2020 <- mean.list.2030 <- mean.list.2040 <- list()
  sd.list.2008 <- sd.list.2020 <- sd.list.2030 <- sd.list.2040 <- list()
  
  for (z in files.r){
    ptmp <- rasterToPoints(raster(z), spatial=T)
    zone.extr.2008 <- as.data.frame(extract(s.2008, ptmp))
    zone.extr.2020 <- as.data.frame(extract(s.2020, ptmp))
    zone.extr.2030 <- as.data.frame(extract(s.2030, ptmp))
    zone.extr.2040 <- as.data.frame(extract(s.2040, ptmp))
    
    #mean and sd over area and gcms
    #wheat.2008.mean <- median(zone.extr.2008[, 1:5], na.rm=T) # median is better
    #wheat.2008.sd <- sd(zone.extr.2008[, 1:5], na.rm=T) 
    
    # ONLY mean and sd over gcms
    idx <- which(z == files.r)
    
    mean.list.2008[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2008[, 1:5]))))
    mean.list.2020[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2020[, 1:5]))))
    mean.list.2030[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2030[, 1:5]))))
    mean.list.2040[[idx]] <- c(median(as.numeric(colwise(median, na.rm=T)(zone.extr.2040[, 1:5]))))
    
    sd.list.2008[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2008[, 1:5])))) 
    sd.list.2020[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2020[, 1:5]))))
    sd.list.2030[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2030[, 1:5]))))
    sd.list.2040[[idx]] <- c(sd(as.numeric(colwise(median, na.rm=T)(zone.extr.2040[, 1:5]))))
  }
  
  all.df.sd <- as.data.frame(cbind(unlist(sd.list.2008),
                                   unlist(sd.list.2020),
                                   unlist(sd.list.2030),
                                   unlist(sd.list.2040)))
  
  colnames(all.df.sd) <- paste(years.ts) # each row = region
  rownames(all.df.sd) <- substr(files.r,1,nchar(files.r) - 4)
  
  all.df.mean <- as.data.frame(cbind(unlist(mean.list.2008),
                                     unlist(mean.list.2020),
                                     unlist(mean.list.2030),
                                     unlist(mean.list.2040)))
  
  colnames(all.df.mean) <- paste(years.ts) # each row = region
  rownames(all.df.mean) <- substr(files.r,1,nchar(files.r) - 4)
  
  setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/2_predictions/allDynamic/bugfix_branch/")
  
  file.name <- paste("Nleach", "_ensemble_mean_zonal_", modality, "_" , i, ".rds", sep="") 
  saveRDS(all.df.mean, file=file.name)
  
  file.name <- paste("Nleach", "_ensemble_sd_zonal_", modality, "_" , i, ".rds", sep="")  
  saveRDS(all.df.sd, file=file.name)
  
} # loop through c("rcp45", "rcp85")





###### Make line graph for BOTH (load data) -----------------------------------------------------------------------------------

### Process yield ###
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/2_predictions/allDynamic/bugfix_branch/")
list.files(pattern="mean_allDynamic")

mean.45 <- readRDS("Yield_ensemble_mean_allDynamic_rcp45.rds")[, c(1,3,5,7,2,4,6,8)]
mean.85 <- readRDS("Yield_ensemble_mean_allDynamic_rcp85.rds")[, c(1,3,5,7,2,4,6,8)]
zs.mean <- as.data.frame(cbind(mean.45, mean.85))

sd.45 <- readRDS("Yield_ensemble_sd_allDynamic_rcp45.rds")[, c(1,3,5,7,2,4,6,8)]
sd.85 <- readRDS("Yield_ensemble_sd_allDynamic_rcp85.rds")[, c(1,3,5,7,2,4,6,8)]
zs.sd.gcms <- as.data.frame(cbind(sd.45, sd.85))

## convert to relative change to 2008 
idx.subtr <- grep("2008", names(zs.mean))
zs.mean.rel.change <- as.data.frame(matrix(nrow = 9, ncol = 16))
upper.rel.change <- as.data.frame(matrix(nrow = 9, ncol = 16))
lower.rel.change <- as.data.frame(matrix(nrow = 9, ncol = 16))

idx <- idx.subtr[1]

for (i in 1:16){
  zs.mean.rel.change[, i] <- ((zs.mean[, i] - zs.mean[, idx]) / zs.mean[, idx] * 100)
  upper.rel.change[, i] <- (((zs.mean[, i] + zs.sd.gcms[, i]) - (zs.mean[, idx])) / (zs.mean[, idx]) * 100) 
  lower.rel.change[, i] <- (((zs.mean[, i] - zs.sd.gcms[, i]) - (zs.mean[, idx])) / (zs.mean[, idx]) * 100)
  
  if (i == 4) idx <- idx.subtr[2]
  if (i == 8) idx <- idx.subtr[3]
  if (i == 12) idx <- idx.subtr[4]
}

colnames(zs.mean.rel.change) <- colnames(zs.mean)
rownames(zs.mean.rel.change) <- rownames(zs.mean)
colnames(upper.rel.change) <- colnames(zs.mean)
rownames(upper.rel.change) <- rownames(zs.mean)
colnames(lower.rel.change) <- colnames(zs.mean)
rownames(lower.rel.change) <- rownames(zs.mean)

# set sd for 2008 to 0
upper.rel.change[, c(1,5,9,13)] <- 0
lower.rel.change[, c(1,5,9,13)] <- 0


### Process Nleaching ###
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/2_predictions/allDynamic/bugfix_branch/")
list.files()
list.files(pattern="mean_allDynamic")

mean.45 <- readRDS("Nleach_ensemble_mean_zonal_allDynamic_rcp45.rds")
mean.85 <- readRDS("Nleach_ensemble_mean_zonal_allDynamic_rcp85.rds")
zs.mean <- as.data.frame(cbind(mean.45, mean.85))

sd.45 <- readRDS("Nleach_ensemble_sd_zonal_allDynamic_rcp45.rds")
sd.85 <- readRDS("Nleach_ensemble_sd_zonal_allDynamic_rcp85.rds")
zs.sd.gcms <- as.data.frame(cbind(sd.45, sd.85))
# delete sd for 2008
zs.sd.gcms[, c(1, 5)] <- 0

## convert to relative change to 2008 
idx.subtr <- grep("2008", names(zs.mean))
zs.mean.rel.change.nleach <- as.data.frame(matrix(nrow = 9, ncol = 8))
upper.rel.change.nleach <- as.data.frame(matrix(nrow = 9, ncol = 8))
lower.rel.change.nleach <- as.data.frame(matrix(nrow = 9, ncol = 8))

idx <- idx.subtr[1]

for (i in 1:8){
  zs.mean.rel.change.nleach[, i] <- ((zs.mean[, i] - zs.mean[, idx]) / zs.mean[, idx] * 100)
  upper.rel.change.nleach[, i] <- (((zs.mean[, i] + zs.sd.gcms[, i]) - (zs.mean[, idx])) / (zs.mean[, idx]) * 100) 
  lower.rel.change.nleach[, i] <- (((zs.mean[, i] - zs.sd.gcms[, i]) - (zs.mean[, idx])) / (zs.mean[, idx]) * 100)
  
  if (i == 4) idx <- idx.subtr[2]
  if (i == 8) idx <- idx.subtr[3]
}

zs.mean.rel.change.nleach[is.na(zs.mean.rel.change.nleach)] <- 0
upper.rel.change.nleach[is.na(upper.rel.change.nleach)] <- 0
lower.rel.change.nleach[is.na(lower.rel.change.nleach)] <- 0

colnames(zs.mean.rel.change.nleach) <- colnames(zs.mean)
rownames(zs.mean.rel.change.nleach) <- rownames(zs.mean)
colnames(upper.rel.change.nleach) <- colnames(zs.mean)
colnames(lower.rel.change.nleach) <- colnames(zs.mean)


### Combine yield and leaching ###
zs.mean.comb <- cbind(zs.mean.rel.change, zs.mean.rel.change.nleach)
zs.upper.comb <- cbind(upper.rel.change, upper.rel.change.nleach)
zs.lower.comb <- cbind(lower.rel.change, lower.rel.change.nleach)
# order: mean wheat, mean maize, upper wheat, upper maize, lower wheat, lower maize, mean leaching, upper leach, lower leach
names(zs.mean.comb) <- paste("lala", 1:24, sep="_")

### Melt and make facet grid plot ###
library(scales)
library(ggthemes)
show_col(tableau_color_pal("tableau20")(20))
cols <- c("#AEC7E8", "#1F77B4","#7F7F7F", "black","#FF9896", "#D62728")
cols <- cols[c(5,6,1,2,3,4)]

### melt for relative change for years 2008:2040 !!!
zones.shrt <- c("Alpine", "Atlantic_C", "Atlantic_N", "Atlantic_S", "Boreal", "Continental_N", "Continental_S", "Mediterranean_N", "Mediterranean_S")

library(reshape)
mean.melt <- reshape2::melt(zs.mean.comb)
mean.melt$year <- rep(rep(c("2008", "2020", "2030", "2040"), each=9), times=6)
mean.melt$zones <- rep(zones.shrt, times=24)
mean.melt$output <- c(rep(rep(c("teww", "teco"), each=36), 2),  rep("leaching", 72))
mean.melt$output <- factor(mean.melt$output, levels=c("teww", "teco", "leaching"))
mean.melt$rcps <- c(rep(c("rcp4.5", "rcp8.5"), each=72), rep(c("rcp4.5", "rcp8.5"), each=36))
mean.melt$groups <- paste(mean.melt$output, mean.melt$rcps, sep="_")
mean.melt$upper <- melt(zs.upper.comb[,])[, 2]
mean.melt$lower <- melt(zs.lower.comb[,])[, 2]

#delete atlantic north for maize
idx.del <- which(mean.melt$zones == "Atlantic_N" & mean.melt$output == "teco")  
mean.melt <- mean.melt[- idx.del, ]
#delete leaching for boreal
#idx.del <- which(mean.melt$zones == "Boreal" & mean.melt$output == "leaching")  
#mean.melt <- mean.melt[- idx.del, ]


# Lineplot, facet wrap for zones
lnplt <- ggplot(data=mean.melt) + 
  geom_ribbon(aes(x=year, ymin=lower, ymax=upper, fill=groups, group=groups), alpha=.15) + 
  geom_line(aes(x=year, y=value, color=groups, group=groups), lwd=1) + 
  facet_grid(output ~ zones, scales="free_y") + 
  geom_point(aes(x=year, y=value, color=groups), size=2.3, shape=21, fill="white", position=position_dodge(0.07)) + 
  theme_bw() + 
  ggtitle("Relative change in yield") + 
  ylab("Relative response (%)") + 
  scale_fill_manual(values=cols, name="RCPs") + 
  geom_hline(aes(yintercept=0), linetype = "dotted") + 
  guides(col=guide_legend(ncol=4)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(strip.background = element_rect(fill = "grey40"), strip.text = element_text(colour = 'white', face="bold")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position = "top") + 
  scale_color_manual(values=cols, name="RCPs")

lnplt

#+ geom_errorbar(aes(x=year, ymin=lower, ymax=upper, color=groups), width=0, position=position_dodge(0.07))


## Save to pdf
# + geom_ribbon(aes(x=year, ymin=lower, ymax=upper, fill=groups, group=groups), alpha=.25)
ggsave(lnplt, file="/home/jan/Dropbox/Paper_2_nitro/Figures/lineplot_allDynamic.pdf", width=12, height=8, units="in")



