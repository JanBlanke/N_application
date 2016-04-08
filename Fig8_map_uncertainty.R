#******************************************
#* Name: Fig8_plot_tradeoff_NUTS2.R       *
#* Author: Jan Blanke                     *
#* Description: Plot categorical tradeoff *
#******************************************

### NOTE: do RCPs and CAPRI/Stocker individually!

library(sp)
library(raster)
library(rasterVis)
library(plyr)
library(rgdal)
library(EcoGenetics)

### Process yield ------------------------------------------------------------------------------------------------

crops <- c("TeWW", "TeCo") # TeCo, TeWW
years.ts <- seq(2008, 2040, by=1) 
years.avg <- 0 # year +- 1 
firstyear <- 1850
output <- "yield"
years <- 2040 # 
rcp <- "rcp85"

## Initialize
df.list.yield <- list()

## read crop frac file for weighing
crop.fracs <- read.table("/media/jan/AddData/Data/Stefan_data/cropland_hurtt_mirca_2000_corrsum.out", head=T)
crop.fracs[, 1:2] <- crop.fracs[, 1:2] + 0.25
head(crop.fracs)
cfracs.sp <- SpatialPixelsDataFrame(crop.fracs[, 1:2], as.data.frame(crop.fracs[, 3:8]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
cfracs.stack <- stack(cfracs.sp)


#setwd("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/allDynamic/output")
setwd("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/Stocker/output")
list.files() # 5 GCMs will be used
file.out <- grep(paste("yield", ".out", sep=""), list.files(), value=T)
file.out <- grep(rcp, file.out, value=T)

gcm.list <- list()

for (j in file.out){ # loop through GCM
  
  #file.to.read <- list.files(pattern = paste(rcp, "_yield", sep=""))
  out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
  
  # Subset according to years
  time.list <- list()
  for (k in years.ts){
    idx.k <- which(years.ts == k)
    
    out.df.sub <- subset(out.df, out.df[3] == k)
    
    out.df.sub <- SpatialPixelsDataFrame(out.df.sub[, 1:2], as.data.frame(out.df.sub[, 4:9]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    out.stack <- stack(out.df.sub) 
    out.stack <- out.stack * 10 / (0.85 * 2.0 * 0.446) # convert to tonnes per hectar and wet matter
    
    if (idx.k == 1){
      frac.mask <- crop(cfracs.stack, out.stack[[1]])
      frac.mask <- mask(frac.mask, out.stack[[1]])
      frac.mask <- frac.mask[[c(2, 5, 1, 4, 3, 6)]]
    }
    
    yield <- ((out.stack[[1]] * frac.mask[[1]]) + (out.stack[[2]] * frac.mask[[2]]) + (out.stack[[3]] * frac.mask[[3]]) + (out.stack[[4]] * frac.mask[[4]]) + (out.stack[[5]] * frac.mask[[5]]) + (out.stack[[6]] * frac.mask[[6]])) / (sum(frac.mask[[1:6]]))
  
    time.list[[idx.k]] <- yield
  }
  
  all.years <- stack(time.list)
  idx.j <- which(file.out == j)
  gcm.list[[j]] <- all.years
  
  
} # end loop GCMs

# calculate mean of different GCms
all.gcms <- stack(gcm.list)

# calculate mean over gcms with equal weights
all.yield <- (all.gcms[[1:length(years.ts)]] + all.gcms[[(length(years.ts) + 1):(length(years.ts) * 2)]] + all.gcms[[(1 + (length(years.ts) * 2)):(length(years.ts) * 3)]] + all.gcms[[(1 + (length(years.ts) * 3)):(length(years.ts) * 4)]] + all.gcms[[(1 + (length(years.ts) * 4)):(length(years.ts) * 5)]]) / 5

# write out
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/6_Tradeoff/")
file.name <- paste("stocker_yield_5gcms_", rcp, "_", length(years.ts), "y", ".nc", sep="") 
writeRaster(all.yield, file=file.name)
  
  
### Process N leaching ----------------------------------------------------------------------------------------------------------

crops <- "leach"
years.ts <- seq(2008, 2040, by=1)  
years.avg <- 0 # year +- 1 
firstyear <- 1850
output <- "nflux" # lai, cpool, cflux, cmass
years <- 2040 # 
rcp <- "rcp85"

## Initialize
df.list.leach <- list()

#setwd("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/allDynamic/output")
setwd("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/Stocker/output")
list.files() # 5 GCMs will be used
file.out <- grep(paste(output, ".out", sep=""), list.files(), value=T)
file.out <- grep(rcp, file.out, value=T)

gcm.list <- list()

for (j in file.out){ # loop through GCM
  
  out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 

  # Subset according to years
  time.list <- list()
  for (k in years.ts){
    
    out.df.sub <- subset(out.df, out.df[3] ==k)
    out.df.sub <- SpatialPixelsDataFrame(out.df.sub[, 1:2], as.data.frame(out.df.sub[, 8]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    out.stack <- stack(out.df.sub) 
    
    idx.k <- which(years.ts == k)
    time.list[[idx.k]] <- out.stack
    
  }
  
  all.years <- stack(time.list)
  idx.j <- which(file.out == j)
  gcm.list[[j]] <- all.years
  
} # end loop GCMs

# calculate mean of different GCms
all.leach <- stack(gcm.list)
 
# average over gcms with equal weights
all.leach <- (all.leach[[1:length(years.ts)]] + all.leach[[(length(years.ts) + 1):(length(years.ts) * 2)]] + all.leach[[(1 + (length(years.ts) * 2)):(length(years.ts) * 3)]] + all.leach[[(1 + (length(years.ts) * 3)):(length(years.ts) * 4)]] + all.leach[[(1 + (length(years.ts) * 4)):(length(years.ts) * 5)]]) / 5

# write out
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/6_Tradeoff/")

file.name <- paste("stocker_leach_5gcms_", rcp, "_", length(years.ts), "y", ".nc", sep="") 
writeRaster(all.leach, file=file.name)


### Create sen slope and p-value rasters -------------------------------------------------------------------------------
rcp <- "rcp85"

### calc trend for Yield
# load stack
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/6_Tradeoff/")

file.out <- grep("stocker_yield", list.files(), value=T)
file.out <- grep(rcp, file.out, value=T)

rasters <- stack(file.out)
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

setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/6_Tradeoff/")
writeRaster(slope * 10, paste("stocker_yield_slope_", rcp, "_", length(years.ts), "y", ".nc", sep=""), overwrite=T)
writeRaster(pvalue, paste("stocker_yield_pval_", rcp, "_", length(years.ts), "y", ".nc", sep=""), overwrite=T)


### calc trend for Leaching
# load stack
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/6_Tradeoff/")

file.out <- grep("stocker_leach", list.files(), value=T)
file.out <- grep(rcp, file.out, value=T)

rasters <- stack(file.out)
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

setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/6_Tradeoff/")
writeRaster(slope * 10, paste("stocker_leach_slope_", rcp, "_", length(years.ts), "y", ".nc", sep=""), overwrite=T)
writeRaster(pvalue, paste("stocker_leach_pval_", rcp, "_", length(years.ts), "y", ".nc", sep=""), overwrite=T)


##############################################################################################################
### Calculate zonal statistics for CAPRI and Stocker##########################################################

# do the processing above for CAPRI and Stocker to just load the below

#eu.nuts2.mask <- readOGR("/home/jan/GIS_data", "NUTS2_wgs84_masked")
#eu.countries.underlay <- readOGR("/home/jan/GIS_data", "EU_27_underlay")
world.borders <- readOGR("/home/jan/GIS_data", "world_borders")
world.borders <- crop(world.borders, extent(-20, 35, 30, 70))
eu.nuts <- readOGR("/home/jan/GIS_data", "NUTS2_wgs84_masked")
eu.nuts <- spTransform(eu.nuts, CRS("+proj=longlat +ellps=WGS84 +no_defs"))

rcp <- "rcp85"

setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/6_Tradeoff/")
slope.yield <- list.files(pattern="yield_slope"); slope.leach <- list.files(pattern="leach_slope")
slope.yield <- grep(rcp, slope.yield, value=T); slope.leach <- grep(rcp, slope.leach, value=T)

s.slope.yield <- stack(slope.yield); s.slope.leach <- stack(slope.leach)
slope.comb <- stack(s.slope.yield, s.slope.leach) #order: stocker_yield, capri_yield, capri_leach, stocker_leach
slope.comb <- slope.comb[[c(2,3,1,4)]] #order: capri_yield, capri_leach, stocker_yield, stocker_leach


## zonal statistics:
zonal.mean.mat <- extract(slope.comb, eu.nuts, weights=T)
zonal.df <- as.data.frame(matrix(nrow = 270, ncol = 4))
for (i in 1:270){ #loop through regions
  df <- as.data.frame(zonal.mean.mat[[i]])
  #product
  df.t <- (df[, 1:4] * df$weight) 
  #sum
  zonal.df[i, 1:4]  <- colSums(df.t[, 1:4], na.rm=T) / sum(df$weight, na.rm=T) #
}

dim(eu.nuts)
eu.nuts$capri_yield <- zonal.df[, 1]
eu.nuts$capri_leach <- zonal.df[, 2]
eu.nuts$stocker_yield <- zonal.df[, 3]
eu.nuts$stocker_leach <- zonal.df[, 4]

spplot(eu.nuts[20:23])
zonal.mean.mat <- eu.nuts[20:23]

#reclassify yield_capri
zonal.mean.mat@data[, 1][zonal.mean.mat@data[, 1] > 0] <- 1
zonal.mean.mat@data[, 1][zonal.mean.mat@data[, 1] < 0] <- 0
#reclassify leaching_capri
zonal.mean.mat@data[, 2][zonal.mean.mat@data[, 2] > 0] <- 1
zonal.mean.mat@data[, 2][zonal.mean.mat@data[, 2] < 0] <- 0
#reclassify yield_stocker
zonal.mean.mat@data[, 3][zonal.mean.mat@data[, 3] > 0] <- 1
zonal.mean.mat@data[, 3][zonal.mean.mat@data[, 3] < 0] <- 0
#reclassify leaching_stocker
zonal.mean.mat@data[, 4][zonal.mean.mat@data[, 4] > 0] <- 1
zonal.mean.mat@data[, 4][zonal.mean.mat@data[, 4] < 0] <- 0

## conditionals for yield
df <- as.data.frame(zonal.mean.mat)
temp <- rep(NA, 270)
temp[df[, 1] == 1 & df[, 3] == 1] <- "increase"
temp[df[, 1] == 0 & df[, 3] == 0] <- "decrease"
temp[df[, 1] == 1 & df[, 3] == 0] <- "uncertain"
temp[df[, 1] == 0 & df[, 3] == 1] <- "uncertain"
zonal.mean.mat@data[, 1] <- as.factor(temp)

## conditionals for leaching
df <- as.data.frame(zonal.mean.mat)
temp <- rep(NA, 270)
temp[df[, 2] == 1 & df[, 4] == 1] <- "increase"
temp[df[, 2] == 0 & df[, 4] == 0] <- "decrease"
temp[df[, 2] == 1 & df[, 4] == 0] <- "uncertain"
temp[df[, 2] == 0 & df[, 4] == 1] <- "uncertain"
zonal.mean.mat@data[, 2] <- as.factor(temp)

zonal.mean.mat <- zonal.mean.mat[1:2]
names(zonal.mean.mat) <- c("yield", "leaching")
levels(zonal.mean.mat@data[, 2]) <- list(increase="increase", uncertain="uncertain", decrease="decrease")
levels(zonal.mean.mat@data[, 1]) <- list(increase="increase", uncertain="uncertain", decrease="decrease")

#eu.nuts <- crop(eu.nuts, eff.stack.m)
library(scales)
library(ggthemes)

desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

show_col(tableau_color_pal("tableau20")(20))
myColors <- c("#FF9896", "grey40", "#AEC7E8") 
myColors <- c("#D62728", "grey40", "#1F77B4") 
cols <- colorRampPalette(brewer.pal(9,"RdBu"))(15)
myColors <- c(cols[4], "grey40", cols[12])

my.settings <- list(strip.background=list(col="grey40"), background=list(col="white"))

## MAP: yield and leaching
png("/home/jan/Dropbox/Paper_2_nitro/Figures/sp_uncertainty.png", width = 20 , height = 13, units = 'cm', res = 600)
spplot(zonal.mean.mat, 
       col.regions=desat(myColors, 1), 
       widths=0.07, 
       col = "grey25", 
       main=NULL, 
       sp.layout=list(world.borders, fill="grey90", col="grey55", lwd=0.25), 
       colorkey = list(space = "bottom", width = 2, height = 0.33),  
       xlim = c(-11, 32), ylim = c(34, 70), 
       lwd=0.07, 
       par.settings = my.settings,
       par.strip.text = list(col = "white", font=2)
)
dev.off()





