############################################################
### Name: 4_model_validation.R                           ###
### Author: Jan Blanke                                   ###
### Description: Regional maize validation for 1999-2006 ###
############################################################

library(rgdal)
library(maptools)
library(rasterVis)

read.preprocessed <- TRUE # read preprocessed or process from scratch

# color settings
source("/home/jan/Dropbox/R/functions/color_ramps.R")
yield.cols <- colorRampPalette(c(brewer.pal(9,"RdYlGn")[c(3, 4, 5, 6, 7, 8, 9)], "darkgreen"))(20)
# shapefiles
eu.ext <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/NUTS_2_shapefiles/", "NUTS2_eurostat_2010")
map_nuts2 <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/", "Europe_country_borders")
bestLPJG <- "/media/jan/AddData/Simulations/Paper_2_nitrogen/TeCo_CRU_validation_simple_manure/" # dir of best LPJG sims

if (read.preprocessed == FALSE){ 
  
  ###--- EUROSTAT data --------------------------------------------------------------------------------------
  # read shape file
  med.nuts <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/NUTS_2_shapefiles", "NUTS2_maizeval_1999_2006_area")
  
  ###--- LPJ-GUESS data --------------------------------------------------------------------------------------
  df <- read.table(paste(bestLPJG, "CMIP5_NORESM1-M_rcp85_yield.out", sep=""), head=T)
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
  yield.wt <- (maize.stack * frac.mask[[1]] + maize.irr.stack * frac.mask[[2]]) / (sum(frac.mask)) 
  
  lu.fracs <- read.table("/media/jan/AddData/Data/crop/lu_for_maize_validation", head=T)
  lu.fracs[, 1:2] <- lu.fracs[, 1:2] + 0.25
  lufracs.sp <- SpatialPixelsDataFrame(lu.fracs[, 1:2], as.data.frame(lu.fracs[, -1:-2]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  lufracs.stack <- stack(lufracs.sp)[[3]]
  lufracs.mask <- crop(lufracs.stack, maize.stack[[1]])
  lufracs.mask <- mask(lufracs.mask, maize.stack[[1]])
  
  yield.wt.lufracs <- stack(yield.wt, lufracs.mask)
  
  ###--- Calculate zonal statistic and combine data ----------------------------------------------------------------------------------
  
  # zonal stats for fraction weighted yield
  weights <- extract(yield.wt.lufracs, med.nuts, weights=T)
  zonal.df <- as.data.frame(matrix(nrow = 105, ncol = 8))
  
  for (i in 1:105){ #loop through regions
    df <- as.data.frame(weights[[i]])
    #product
    df.t <- (df[, 1:8] * df$weight * df$CROPLAND) #
    #sum
    zonal.df[i, 1:8]  <- colSums(df.t[, 1:8], na.rm=T) / sum(df$weight * df$CROPLAND, na.rm=T) #
  }
  
  # combine data
  zonal.data <- med.nuts
  row.names(zonal.df) <- as.character(zonal.data$NUTS2)
  row.names(zonal.data) <- as.character(zonal.data$NUTS2)
  zonal.data <- spCbind(zonal.data, zonal.df)
  names(zonal.data)[18:25] <- paste("r", 1999:2006, sep="")
  
  # calculate mean obs and mean mod yield
  zonal.data$meanObs <- rowMeans(as.data.frame(zonal.data[, 7:14]))
  zonal.data$meanMod <- rowMeans(as.data.frame(zonal.data[, 18:25]))
  zonal.data$meanDiff <- zonal.data$meanMod - zonal.data$meanObs
  
  idx <- which(zonal.data$meanDiff >= 6)
  zonal.data[idx,] <- NA
  
  ## Write data to hard drive
  save(zonal.data, file="/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/validation/data_to_plot.RData")

}

#--- Make maps for paper ------------------------------------------------------------------------------------------------------------
load("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/1_validation/data_to_plot.RData")

mO <- zonal.data["meanObs"] 
pObs <- spplot(mO, lwd=0.3, col.regions=yield.cols, main=NULL, sp.layout=list(eu.ext, fill="grey90", col="grey55", lwd=0.15), at=seq(2, 12.2, length=20), colorkey = list(space = "bottom", width = 2, height = 1), xlim = bbox(mO)[1, ] + c(0, 4), ylim = bbox(mO)[2, ] + c(0, 5))

maxval <- max(abs(range(zonal.data["meanDiff"]@data, na.rm=T)))
maxval <- 4.3
mD <- zonal.data["meanDiff"]
pDiff <- spplot(mD, lwd=0.3, col.regions=rdbuPal, main=NULL, at=seq(-maxval, maxval, length=20), colorkey = list(space = "bottom", width = 2, height = 1), xlim = bbox(mD)[1, ] + c(0, 4), ylim = bbox(mD)[2, ] + c(0, 5))

library(grid)
zonal.data.df <- as.data.frame(zonal.data)
pScatter <- ggplot(zonal.data.df, aes(x=meanObs, y=meanMod, weight=meanArea, size=meanArea)) + geom_abline(intercept=0, slope=1, linetype = "dotted") + geom_smooth(method="lm", color="darkorange", lwd=0.4) + geom_point(alpha = 0.5, color="white", shape=21, fill="black") + scale_size_area() + theme_bw() + scale_y_continuous(limits=c(2.3, 14.7)) + scale_x_continuous(limits=c(2.3, 12.1)) + theme(plot.margin=unit(c(-0.5,-0.5,-0.5,-0.5),"cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "transparent",colour = NA)) + theme(legend.position = "none") + scale_size(range = c(1,8)) + xlab("") + ylab("")
pScatter


# viewport approach
#pdf("/home/jan/Results/Paper_2_intensity/maize_val_current_bubble.pdf", width=27/2.54 , height=20/2.54)
png("/home/jan/Dropbox/Paper_2_nitro/Figures/maize_val.png", width = 27 , height = 18, units = 'cm', res = 600)
grid.newpage()
vp1 <- viewport(x = 0, y = 1, height = 1, width = 0.5, just = c("left", "top"))
pushViewport(vp1)
print(pObs, newpage = FALSE)
upViewport(1)
vp2 <- viewport(x = 0.5, y = 1, height = 1, width = 0.5, just = c("left", "top"))
pushViewport(vp2)
print(pDiff, newpage = FALSE) 
upViewport(1)
vp3 <- viewport(x = 0.95, y = 0.85, height = 0.28, width = 0.19, just = c("right", "top"))
pushViewport(vp3)
print(pScatter, newpage = FALSE) 
upViewport(1)
dev.off() 










