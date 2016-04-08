#****************************************************************************************
#* Name: Fig6_sen_slope_factors_NUTS2.R                                                 *
#* Author: Jan Blanke                                                                   *
#* Description: Plot spatial factors attribution (after Chang 2015)                     *
#* Dependencies: /home/jan/Dropbox/R/2nd_paper/N_application/prepare_5gcms_sens_slope.R * 
#*****************************************************************************************

# run following script first: "/home/jan/Dropbox/R/2nd_paper/N_application/prepare_5gcms_sens_slope.R"

library(raster)
library(rgdal)
library(EcoGenetics)

eu.nuts2.mask <- readOGR("/home/jan/GIS_data", "NUTS2_wgs84_masked")
eu.countries.underlay <- readOGR("/home/jan/GIS_data", "EU_27_underlay")
world.borders <- readOGR("/home/jan/GIS_data", "world_borders")
world.borders <- crop(world.borders, extent(-20, 35, 30, 70))#crop world borders

eu.nuts <- readOGR("/home/jan/GIS_data", "NUTS2_wgs84_masked")
eu.nuts <- spTransform(eu.nuts, CRS("+proj=longlat +ellps=WGS84 +no_defs"))
rcp <- "rcp85"

loopstr <- c("allDynamic", "Climconst", "CO2const", "Nconst")


###### Create sen slope rasters ------------------------------------------------------------------------------------------------------

### loop through combinations for WHEAT
for (i in loopstr){

  # load stack
  setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/5GCMs_netcdf/bugfix_33y")
  
  file.out <- grep(paste(i), list.files(), value=T)
  file.out <- grep("wheat", file.out, value=T)
  file.out <- grep(rcp, file.out, value=T)
  
  tempstack <- stack(file.out)
  rasters <- tempstack
  dates <- seq(2008, 2040, by=1)
  
  setwd("/home/jan/tmp/")
  idx.is <- na.idx <- which(is.na(rasters[[1]][]))
  rasters <- brick(rasters)
  
  for (j in 1:length(dates)) rasters[[j]][idx.is] <- 999
    
  eco.theilsen(rasters, dates)
  
  slope <- raster("slope.tif"); slope[idx.is] <- NA
  plot(slope * 10)
  pvalue <- raster("pvalue.tif"); pvalue[idx.is] <- NA
  plot(pvalue)
  
  sig.idx <- which(getValues(pvalue) <= 0.1)
  pvalue[sig.idx] <- 2
  sig.idx <- which(getValues(pvalue) < 2)
  pvalue[sig.idx] <- 0
  
  setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/bugfix_slope/")
  writeRaster(slope * 10, paste("slope_wheat_", i, "_", rcp, ".tif", sep=""), overwrite=T)
  writeRaster(pvalue, paste("pval_wheat_", i, "_", rcp,".tif", sep=""), overwrite=T)
}


### loop through combinations for MAIZE
for (i in loopstr){
  
  # load stack
  setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/5GCMs_netcdf/bugfix_33y")
  
  file.out <- grep(paste(i), list.files(), value=T)
  file.out <- grep("maize", file.out, value=T)
  file.out <- grep(rcp, file.out, value=T)
  
  tempstack <- stack(file.out)
  rasters <- tempstack
  dates <- seq(2008, 2040, by=1)
  
  setwd("/home/jan/tmp/")
  idx.is <- na.idx <- which(is.na(rasters[[1]][]))
  rasters <- brick(rasters)
  
  for (j in 1:length(dates)) rasters[[j]][idx.is] <- 999
  
  eco.theilsen(rasters, dates)
  
  slope <- raster("slope.tif"); slope[idx.is] <- NA
  plot(slope * 10)
  pvalue <- raster("pvalue.tif"); pvalue[idx.is] <- NA
  plot(pvalue)
  
  sig.idx <- which(getValues(pvalue) <= 0.1)
  pvalue[sig.idx] <- 2
  sig.idx <- which(getValues(pvalue) < 2)
  pvalue[sig.idx] <- 0
  
  setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/bugfix_slope/")
  writeRaster(slope * 10, paste("slope_maize_", i, "_", rcp, ".tif", sep=""), overwrite=T)
  writeRaster(pvalue, paste("pval_maize_", i, "_", rcp, ".tif", sep=""), overwrite=T)
  
}


### loop through combinations for NLEACH
for (i in loopstr){
  
  # load stack
  setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/5GCMs_netcdf/bugfix_33y")
  
  file.out <- grep(paste(i), list.files(), value=T)
  file.out <- grep("leach", file.out, value=T)
  file.out <- grep(rcp, file.out, value=T)
  
  tempstack <- stack(file.out)
  rasters <- tempstack 
  dates <- seq(2008, 2040, by=1)
  
  setwd("/home/jan/tmp/")
  idx.is <- na.idx <- which(is.na(rasters[[1]][]))
  rasters <- brick(rasters)
  
  for (j in 1:length(dates)) rasters[[j]][idx.is] <- 999
  
  eco.theilsen(rasters, dates)
  
  slope <- raster("slope.tif"); slope[idx.is] <- NA
  plot(slope * 10)
  pvalue <- raster("pvalue.tif"); pvalue[idx.is] <- NA
  plot(pvalue)
  
  sig.idx <- which(getValues(pvalue) <= 0.1)
  pvalue[sig.idx] <- 2
  sig.idx <- which(getValues(pvalue) < 2)
  pvalue[sig.idx] <- 0
  
  setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/bugfix_slope/")
  writeRaster(slope * 10, paste("slope_leach_", i, "_", rcp, ".tif", sep=""), overwrite=T)
  writeRaster(pvalue, paste("pval_leach_", i, "_", rcp, ".tif", sep=""), overwrite=T)
  
}


###### Subract and plot trend rasters to attribute effects (according to Chang et al. 2015) --------------------------------------

library(rasterVis)
library(rgdal)
source("/home/jan/Dropbox/R/functions/color_ramps.R")

### Process WHEAT and MAIZE
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/bugfix_slope/")
lsf.slope.wheat <- list.files(pattern="slope_wheat"); lsf.slope.maize <- list.files(pattern="slope_maize")
lsf.slope.wheat <- grep(rcp, lsf.slope.wheat, value=T); lsf.slope.maize <- grep(rcp, lsf.slope.maize, value=T)

lsf.pval.wheat <- list.files(pattern="pval_wheat"); lsf.pval.maize <- list.files(pattern="pval_maize")
lsf.pval.wheat <- grep(rcp, lsf.pval.wheat, value=T); lsf.pval.maize <- grep(rcp, lsf.pval.maize, value=T)

s.slope.wheat <- stack(lsf.slope.wheat); s.slope.maize <- stack(lsf.slope.maize)
s.p.wheat <- stack(lsf.pval.wheat); s.p.maize <- stack(lsf.pval.maize)

# prepare effect rasters for factors
print(loopstr)

all.dyn.w <- s.slope.wheat[[1]]; all.dyn.m <- s.slope.maize[[1]]
n.eff.w <- s.slope.wheat[[1]] - s.slope.wheat[[4]]; n.eff.m <- s.slope.maize[[1]] - s.slope.maize[[4]]
co2.eff.w <- s.slope.wheat[[1]] - s.slope.wheat[[3]]; co2.eff.m <- s.slope.maize[[1]] - s.slope.maize[[3]]
clim.eff.w <- s.slope.wheat[[1]] - s.slope.wheat[[2]]; clim.eff.m <- s.slope.maize[[1]] - s.slope.maize[[2]]

eff.stack.w <- stack(all.dyn.w, n.eff.w, co2.eff.w, clim.eff.w); eff.stack.m <- stack(all.dyn.m, n.eff.m, co2.eff.m, clim.eff.m)
names(eff.stack.w) <- c("dynamic", "n.eff", "co2.eff", "clim.eff"); names(eff.stack.m) <- c("dynamic", "n.eff", "co2.eff", "clim.eff")


### Process NLEACHING
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/bugfix_slope/")
list.files()

lsf.slope.nleach <- list.files(pattern="slope_leach")
lsf.slope.nleach <- grep(rcp, lsf.slope.leach, value=T)

s.slope.leach <- stack(lsf.slope.nleach)

# prepare effect rasters for factors
all.dyn.l <- s.slope.leach[[1]]
n.eff.l <- s.slope.leach[[1]] - s.slope.leach[[4]]
co2.eff.l <- s.slope.leach[[1]] - s.slope.leach[[3]]
clim.eff.l <- s.slope.leach[[1]] - s.slope.leach[[2]]

eff.stack.l <- stack(all.dyn.l, n.eff.l, co2.eff.l, clim.eff.l)
names(eff.stack.l) <- c("dynamic", "n.eff", "co2.eff", "clim.eff")


### Calculate WHEAT zonal statistics:
zonal.mean.mat <- extract(eff.stack.w, eu.nuts, weights=T)
zonal.df <- as.data.frame(matrix(nrow = 270, ncol = 4))
for (i in 1:270){ #loop through regions
  df <- as.data.frame(zonal.mean.mat[[i]])
  #product
  df.t <- (df[, 1:4] * df$weight) 
  #sum
  zonal.df[i, 1:4]  <- colSums(df.t[, 1:4], na.rm=T) / sum(df$weight, na.rm=T) #
}

eu.nuts$wheat_dynamic <- zonal.df[, 1]
eu.nuts$wheat_neff <- zonal.df[, 2]
eu.nuts$wheat_co2eff <- zonal.df[, 3]
eu.nuts$wheat_climeff <- zonal.df[, 4]

### Map for WHEAT
my.settings <- list(strip.background=list(col="grey40"), background=list(col="white"), layout.heights=list(top.padding=-2, bottom.padding=-1))

maxval <- max(as.data.frame(eu.nuts[20:23]))
layout.1 <- list(world.borders, fill="grey90", col="grey55", lwd=0.5)
layout.2 <- list("sp.text", c(-7, 68), "a)",  which = 1) 

wp <- spplot(eu.nuts[20:23], 
             col.regions=rev(rdbuPal), 
             layout=c(4,1), 
             col = "grey25", 
             main=NULL, 
             sp.layout=list(layout.1, layout.2), 
             colorkey = list(space = "right", width = 2, height = 1),  
             xlim = c(-11, 32), ylim = c(34, 70), 
             lwd=0.07, 
             par.settings = my.settings, 
             par.strip.text = list(col = "white", font=2),
             at=seq(-maxval, maxval, length=15), 
             names.attr=c("dynamic", "n.eff", "co2.eff", "clim.eff"))
wp


### Calculate MAIZE zonal statistics:
zonal.mean.mat <- extract(eff.stack.m, eu.nuts, weights=T)
zonal.df <- as.data.frame(matrix(nrow = 270, ncol = 4))
for (i in 1:270){ #loop through regions
  df <- as.data.frame(zonal.mean.mat[[i]])
  #product
  df.t <- (df[, 1:4] * df$weight) 
  #sum
  zonal.df[i, 1:4]  <- colSums(df.t[, 1:4], na.rm=T) / sum(df$weight, na.rm=T) #
}

eu.nuts$maize_dynamic <- zonal.df[, 1]
eu.nuts$maize_neff <- zonal.df[, 2]
eu.nuts$maize_co2eff <- zonal.df[, 3]
eu.nuts$maize_climeff <- zonal.df[, 4]

### MAP for MAIZE
my.settings <- list(strip.background=list(col="grey40"), background=list(col="white"), layout.heights=list(top.padding=-2, bottom.padding=-1))

maxval <- max(as.data.frame(eu.nuts[24:27]))
layout.1 <- list(world.borders, fill="grey90", col="grey55", lwd=0.5)
layout.2 <- list("sp.text", c(-7, 68), "b)",  which = 1) 

mp <- spplot(eu.nuts[24:27],
             col.regions=rev(rdbuPal), 
             widths=0.07, 
             layout=c(4,1), 
             col = "grey25", 
             main=NULL, 
             sp.layout=list(layout.1, layout.2), 
             colorkey = list(space = "right", width = 2, height = 1),  
             xlim = c(-11, 32), ylim = c(34, 70), 
             lwd=0.07, 
             par.settings = my.settings, 
             par.strip.text = list(col = "white", font=2),
             at=seq(-maxval, maxval, length=15), 
             names.attr=c("dynamic", "n.eff", "co2.eff", "clim.eff"))
mp


### Calculate NLEACHING zonal statistics:
zonal.mean.mat <- extract(eff.stack.l, eu.nuts, weights=T)
zonal.df <- as.data.frame(matrix(nrow = 270, ncol = 4))
for (i in 1:270){ #loop through regions
  df <- as.data.frame(zonal.mean.mat[[i]])
  #product
  df.t <- (df[, 1:4] * df$weight) 
  #sum
  zonal.df[i, 1:4]  <- colSums(df.t[, 1:4], na.rm=T) / sum(df$weight, na.rm=T) 
}

#delete leaching values above 20
zonal.df[zonal.df>=20] <- NA

eu.nuts$leaching_dynamic <- zonal.df[, 1]; 
eu.nuts$leaching_neff <- zonal.df[, 2]
eu.nuts$leaching_co2eff <- zonal.df[, 3]
eu.nuts$leaching_climeff <- zonal.df[, 4]

### MAP for NLEACHING
my.settings <- list(strip.background=list(col="grey40"), background=list(col="white"), layout.heights=list(top.padding=-2, bottom.padding=-1))

maxval <- max(as.data.frame(eu.nuts[28:31]), na.rm=T)
layout.1 <- list(world.borders, fill="grey90", col="grey55", lwd=0.5)
layout.2 <- list("sp.text", c(-7, 68), "c)",  which = 1) 

#pdf("/home/jan/Dropbox/Paper_2_nitro/effect_wheat_zonal.pdf", width=20/2.54 , height=29/2.54)
lp <- spplot(eu.nuts[28:31], 
             col.regions=rev(rdbuPal), 
             layout=c(4,1), 
             col = "grey25", main=NULL, sp.layout=list(layout.1, layout.2), 
             colorkey = list(space = "right", width = 2, height = 1),  
             xlim = c(-11, 32), ylim = c(34, 70), 
             lwd=0.07, 
             par.settings = my.settings, 
             par.strip.text = list(col = "white", font=2),
             at=seq(-maxval, maxval, length=15),
             names.attr=c("dynamic", "n.eff", "co2.eff", "clim.eff"))
lp
#dev.off()


### MAP WHEAT + MAIZE + NLEACH
#width = 5.9; pw = list(x=width, units="cm") 
png("/home/jan/Dropbox/Paper_2_nitro/Figures/attribution_trend_ensemble_nut2_rcp85_test.png", width = 27 , height = 27  , units = 'cm', res = 600)
print(wp, split=c(1,1,1,3), more=T) 
print(mp, split=c(1,2,1,3), more=T)
print(lp, split=c(1,3,1,3), more=F)
dev.off()


### get percentages positive trends for paper ---------------------------------------------------------------------------------------------------
tmp <- as.data.frame(eu.nuts[20])
(length(which(tmp >= 0)) / nrow(tmp)) * 100

tmp <- as.data.frame(eu.nuts[24])
(length(which(tmp >= 0)) / nrow(tmp)) * 100

tmp <- as.data.frame(eu.nuts[28])
(length(which(tmp < 0)) / nrow(tmp)) * 100








