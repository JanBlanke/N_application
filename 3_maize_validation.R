############################################################
### Name: 3_maize_validation.R                           ###
### Author: Jan Blanke                                   ###
### Description: Regional maize validation for 1999-2006 ###
############################################################


#--- Settings ------------------------------------------------------------
filename <- "val_results.pdf" #file name of pdf containing all plots
dir <- "/home/jan/LPJ_GUESS/Paper_2/crop_climate_run/" # dir of LPJG sims
# -----------------------------------------------------------------------

library(rgdal)
library(maptools)
library(rasterVis)

# color settings
source("/home/jan/Dropbox/R/functions/colorramps.R")
yield.cols <- colorRampPalette(c(brewer.pal(9,"RdYlGn")[c(3, 4, 5, 6, 7, 8, 9)], "darkgreen"))(100)
  
# shapefiles
eu.ext <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/NUTS_2_shapefiles/", "NUTS2_eurostat_2010")

setwd("/home/jan/Downloads/")
list.files(pattern="apro")


#--- EUROSTAT data --------------------------------------------------------------------------------------
# read shape file
med.nuts <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/NUTS_2_shapefiles", "NUTS2_maizeval_1999_2006")


#--- LPJ-GUESS data --------------------------------------------------------------------------------------

df <- read.table(paste(dir, "CMIP5_NORESM1-M_rcp85_yield.out", sep=""), head=T)
df.sub <- subset(df, Year >= 598 & Year <= 605) # 1999:2006
df.sub[, 1:2] <- df.sub[, 1:2] + 0.25

maize.irr.list <- list(); maize.list <- list(); j <- 1

for (i in 598:605){  
  df.temp <- subset(df.sub, Year == i)
    maize.sp <- SpatialPixelsDataFrame(df.temp[, 1:2], as.data.frame(df.temp[, 8]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  maize.irr.sp <- SpatialPixelsDataFrame(df.temp[, 1:2], as.data.frame(df.temp[, 9]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  maize.irr.list[[j]] <- raster(maize.irr.sp)
  maize.list[[j]] <- raster(maize.sp)
  j <- j + 1  
}
  
maize.stack <- stack(maize.list) * 10 / (0.85 * 2.0 * 0.446)
maize.irr.stack <- stack(maize.irr.list) * 10 / (0.85 * 2.0 * 0.446)
maize.stack <- reclassify(maize.stack, matrix(c(0, NA), ncol=2, byrow=TRUE))
maize.irr.stack <- reclassify(maize.irr.stack, matrix(c(0, NA), ncol=2, byrow=TRUE))
names(maize.stack) <- paste("TeCo", 1999:2006, sep="")
names(maize.irr.stack) <- paste("TeCoI", 1999:2006, sep="")

## Read and process fraction data
crop.fracs <- read.table("/media/jan/AddData/Data/Stefan_data/cropland_hurtt_mirca_2000_corrsum.out", head=T)
crop.fracs[, 1:2] <- crop.fracs[, 1:2] + 0.25
head(crop.fracs)
cfracs.sp <- SpatialPixelsDataFrame(crop.fracs[, 1:2], as.data.frame(crop.fracs[, c(5,8)]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
cfracs.stack <- stack(cfracs.sp)
frac.mask <- crop(cfracs.stack, maize.stack[[1]])
frac.mask <- mask(frac.mask, maize.stack[[1]])

## calculate weighted maize yield based on crop fractions
yield.wt <- (maize.stack * frac.mask[[1]] + maize.irr.stack * frac.mask[[2]]) / (sum(frac.mask)) # to early to cal weighted yield already here?


#--- Calculate zonal statistic and combine data -------------------------------------------------------------------------------------------------

# zonal stats for fraction weighted yield
weights <- extract(yield.wt, med.nuts, weights=T)
zonal.df <- as.data.frame(matrix(nrow = 105, ncol = 8))

for (i in 1:105){
  df <- as.data.frame(weights[[i]])
  df <- df[, 1:8] * df[, 9]
  zonal.df[i, 1:8] <- colSums(df, na.rm=T)  
}

# combine data
zonal.data <- med.nuts
row.names(zonal.df) <- as.character(zonal.data$NUTS2)
row.names(zonal.data) <- as.character(zonal.data$NUTS2)
zonal.data <- spCbind(zonal.data, zonal.df)
names(zonal.data)[17:24] <- paste("r", 1999:2006, sep="")

# calculate mean obs and mean mod yield
zonal.data$meanObs <- rowMeans(as.data.frame(zonal.data[, 7:14]))
zonal.data$meanMod <- rowMeans(as.data.frame(zonal.data[, 17:24]))
zonal.data$meanDiff <- zonal.data$meanObs - zonal.data$meanMod


#--- Make maps of mean difference ------------------------------------------------------------------------------------------------------------

PDFPath <- paste("/home/jan/Results/Paper_2_intensity/", filename, sep="")
pdf(file=PDFPath) 

spplot(zonal.data["meanObs"], col.regions=yield.cols, main=list(label="Maize EUROSTAT"), sp.layout=list(eu.ext, col.regions= "grey85"), at=seq(2,15, length=20))
spplot(zonal.data["meanMod"], col.regions=yield.cols, main=list(label="Maize CRU LPJG"), sp.layout=eu.ext, at=seq(2,15, length=20))
spplot(zonal.data["meanDiff"], col.regions=rdbuPal, main=list(label="Maize CRU LPJG"), at=seq(-10,10, length=20))

#pw = list(x= 12, units="cm") 
#png("//home/jan/Results/Paper_2_intensity/maize_val_mean_maps.png", width = 27 , height = 15  , units = 'cm', res = 600)
#print(p1, split = c(1, 1, 2, 1), panel.width=pw, more = TRUE)
#print(p2, split = c(2, 1, 2, 1), panel.width=pw, more = FALSE)
#dev.off()

#--- RMSE and Willmotts between all years -----------------------------------------------------------------------------------------------------
library(openair)

RMSE <- vector()
Willmotts <- vector()
j <- 1

for (i in 1999:2006){
  ModStats <- modStats(mydata=as.data.frame(zonal.data), mod=paste("X", i, sep=""), obs=paste("r", i, sep=""))
  RMSE[j] <- ModStats$RMSE
  Willmotts[j] <- ModStats$IOA 
  j <- j + 1
}

plot(RMSE, pch=19)
plot(Willmotts, pch=19)
  
#---Line plots --------------------------------------------------------------------------------------------------------------------------------

# 1-1 line plot
#png("//home/jan/Results/Paper_2_intensity/maize_val_lineplots.png", width = 27 , height = 15  , units = 'cm', res = 600)
par(mfrow=c(2,4))
plot(zonal.data$X1999, zonal.data$r1999); abline(1,1)
plot(zonal.data$X2000, zonal.data$r2000); abline(1,1)
plot(zonal.data$X2001, zonal.data$r2001); abline(1,1)
plot(zonal.data$X2002, zonal.data$r2002); abline(1,1)
plot(zonal.data$X2003, zonal.data$r2003); abline(1,1)
plot(zonal.data$X2004, zonal.data$r2004); abline(1,1)
plot(zonal.data$X2005, zonal.data$r2005); abline(1,1)
plot(zonal.data$X2006, zonal.data$r2006); abline(1,1)

dev.off() #pdf()











