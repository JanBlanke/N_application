########################################################
### Name: 2_prepare_ensemble.R                       ###
### Author: Jan Blanke                               ###
### Description: Process LPJ-GUESS sims for analysis ###
########################################################

library(rgdal)
library(rasterVis)
library(plyr)

eu.countries.mask <- readOGR("/home/jan/GIS_data", "EU_match_clue")
eu.countries.underlay <- readOGR("/home/jan/GIS_data", "EU_27_underlay")
asseng.points <- readOGR("/home/jan/Dropbox/GIS_layers", "asseng_locations")

source("/home/jan/Dropbox/R/functions/colorramps.R")

######## Modify ########
#setwd("/media/jan/AddData/Simulations/Paper_2_nitrogen/ensemble/") # ensemble run with N and CO2 switched on
#setwd("/media/jan/AddData/Simulations/Paper_2_nitrogen/ensemble_constantN") # N switched off
setwd("/media/jan/AddData/Simulations/Paper_2_nitrogen/ensemble_constantCO2N") # CO2 switched off
modality <- "constantCO2" # allDynamic, constantN, constantCO2
data.out <- "rcp85_nflux" # yield, nflux, rcp45, rcp85
crop <- "TeWW" # TeCo, TeWW
years.ts <- c(2008, 2020, 2030, 2040)
years.avg <- 2 # year +- 2 
firstyear <- 1850
########################

list.files() # 5 GCMs will be used
file.out <- grep(paste(data.out, ".out", sep=""), list.files(), value=T)

gcm.list <- list()

for (j in file.out){

  out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
  cft.stack.list <- list()
  
  
  for (i in years.ts){
    
    years <- (i - years.avg):(i + years.avg)
    modelyears <- years - firstyear + 500
    df.out.sub <- subset(out.df, out.df[3] >= min(modelyears) & out.df[3] <= max(modelyears))
    
    if (length(years) > 1 && grepl("yield", data.out)) {
      
      df.out.sub <- ddply(df.out.sub, c(.(Lon), .(Lat)), summarise, mean.teww = mean(TeWW, na.rm = TRUE), mean.tewwi = mean(TeWWi, na.rm = TRUE), mean.tesw = mean(TeSW, na.rm = TRUE), mean.teswi = mean(TeSWi, na.rm = TRUE), mean.teco = mean(TeCo, na.rm = TRUE), mean.tecoi = mean(TeCoi, na.rm = TRUE))
      df.out.sub <- cbind(df.out.sub[, 1:2], rep(mean(years), 2199), df.out.sub[, 3:ncol(df.out.sub)]); colnames(df.out.sub)[3] <- "Year"
    
    } else {
      
      # currently: N leaching!
      df.out.sub <- ddply(df.out.sub, c(.(Lon), .(Lat)), summarise, mean.leach = mean(leach, na.rm = TRUE))
      df.out.sub <- cbind(df.out.sub[, 1:2], rep(mean(years), 2199), df.out.sub[, 3:ncol(df.out.sub)]); colnames(df.out.sub)[3] <- "Year" 
      
    }
    
    sp.df.out <- SpatialPixelsDataFrame(df.out.sub[, 1:2], as.data.frame(df.out.sub[, 4:ncol(df.out.sub)]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    out.stack <- stack(sp.df.out)
    out.stack <- reclassify(out.stack, matrix(c(0, NA), ncol=2, byrow=TRUE))
    names(out.stack) <- names(sp.df.out)
    
    if (crop == "TeWW" && grepl("yield", data.out)) {
      cft <- sum(out.stack[[c(1, 2, 3, 4)]])
      cft <- (cft / 0.45) * (1 / 0.85)
    } else if (crop == "TeCo" && grepl("yield", data.out)) {
      cft <- sum(out.stack[[c(5, 6)]])
      cft <- (cft / 0.45) * (1 / 0.85)      
    }  else {
      cft <- out.stack
    }  
    
    idx <- which(i == years.ts)
    cft.stack.list[[idx]] <- cft
  }
  
  all.stack <- stack(cft.stack.list)
  names(all.stack) <- paste(crop, ".", years.ts, sep="")
  idxy <- which(j == file.out)
  gcm.list[[idxy]] <- all.stack
  
}

### calculate ensemble mean for each time slice
s.2008 <- stack(gcm.list[[1]][[1]], gcm.list[[2]][[1]], gcm.list[[3]][[1]], gcm.list[[4]][[1]], gcm.list[[5]][[1]])
s.2020 <- stack(gcm.list[[1]][[2]], gcm.list[[2]][[1]], gcm.list[[3]][[2]], gcm.list[[4]][[2]], gcm.list[[5]][[2]])
s.2030 <- stack(gcm.list[[1]][[3]], gcm.list[[2]][[1]], gcm.list[[3]][[3]], gcm.list[[4]][[3]], gcm.list[[5]][[3]])
s.2040 <- stack(gcm.list[[1]][[4]], gcm.list[[2]][[1]], gcm.list[[3]][[4]], gcm.list[[4]][[4]], gcm.list[[5]][[4]])

mean.2008 <- mean(s.2008, na.rm=T) 
mean.2020 <- mean(s.2020, na.rm=T) 
mean.2030 <- mean(s.2030, na.rm=T) 
mean.2040 <- mean(s.2040, na.rm=T) 

sd.fun = function(x) sqrt(var(x))
sd.2008 <- calc(s.2008, sd.fun)
sd.2020 <- calc(s.2020, sd.fun)
sd.2030 <- calc(s.2030, sd.fun)
sd.2040 <- calc(s.2040, sd.fun)

ensemble.mean.stack <- stack(mean.2008, mean.2020, mean.2030, mean.2040)
names(ensemble.mean.stack) <- paste(crop, ".", years.ts, sep="")
ensemble.sd.stack <- stack(sd.2008, sd.2020, sd.2030, sd.2040)
names(ensemble.sd.stack) <- paste(crop, ".", years.ts, sep="")

### Plots
setwd("/home/jan/Results/Paper_2_intensity/nleach_rasters/")

# ensemble mean
#maxval.mean <- max(maxValue(ensemble.mean.stack))
#minval.mean <- min(minValue(ensemble.mean.stack))
#levelplot(ensemble.mean.stack, par.settings=yieldTheme, main="Climate: Ensemble RCP4.5, Nitrogen: VOLANTE B1", at=seq(minval.mean, maxval.mean, length=100)) + layer(sp.polygons(eu.countries.underlay, lwd=0.5, alpha=1, fil="grey90", col="grey75"), under=T) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black")) 

file.name <- paste(crop, "_", substr(data.out, 7, nchar(data.out)), "_ensemble_mean_", modality, "_" , substr(data.out, 1, 5), ".nc", sep="") 
writeRaster(ensemble.mean.stack, file=file.name)


# ensemble sd
#maxval.sd <- max(maxValue(ensemble.sd.stack))
#levelplot(ensemble.sd.stack, par.settings=purpleTheme, main="Climate: SD RCP4.5, Nitrogen: VOLANTE B1", at=seq(0, maxval.sd, length=100)) + layer(sp.polygons(eu.countries.underlay, lwd=0.5, alpha=1, fil="grey75", col="white"), under=T) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black"), under=F)

file.name <- paste(crop, "_", substr(data.out, 7, nchar(data.out)), "_ensemble_sd_", modality, "_" , substr(data.out, 1, 5), ".nc", sep="") 
writeRaster(ensemble.mean.stack, file=file.name)

# difference raster
#r.diff <- ensemble.mean.stack[[4]] - ensemble.mean.stack[[1]]
#divThemeOpt <- rasterTheme(region=divRampFun(r.diff))
#levelplot(r.diff, par.settings=divThemeOpt, at=seq(min(r.diff[], na.rm = T), max(r.diff[], na.rm = T), length=100)) + layer(sp.polygons(eu.countries.underlay, #lwd=1, alpha=1, fil="grey85", col="white"), under=T) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black")) 

# Points from Asseng
#+ layer(sp.points(asseng.points, col="red", lwd=5))



