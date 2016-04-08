### Name: 6_bargraph_rel_importance.R
### Author: Jan Blanke
### Description: Make plots of anova analysis
########################################################################
library(reshape)
library(ggplot2)
library(sp)
library(raster)
library(rasterVis)
library(plyr)

### calculate ANOVA sensitivity for yield (save as .rds). do each RCP individually! ------------------------------------------------------------

## Settings
output <- "yield" # lai, cpool, cflux, cmass
years <- 2040 # 
rcp <- "rcp85"
calc.mean <- "zonal" #zonal

## Initialize
df.list.wheat <- list()
df.list.maize <- list()

loopstr <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/factorial_design", pattern="clim")

for (i in loopstr){ # loop through 8 combinations 
  
  lsf <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/factorial_design/nOFF_co2OFF_climOFF/output")
  file.out <- grep(paste("yield", ".out", sep=""), lsf, value=T)
  file.out <- grep(rcp, file.out, value=T)
  
  gcm.list <- list()
  
  for (j in file.out){ # loop through GCM
    
    setwd(paste("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/factorial_design/", i, "/output/", sep=""))
    
    file.to.read <- list.files(pattern = paste(rcp, "_yield", sep=""))
    
    out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
    
    ### Convert to spatial objects
    # Subset according to years
    #years.tmp <- years - 1850 + 500
    
    out.df <- subset(out.df, out.df[3] >= min(years) & out.df[3] <= max(years))
    
    out.df <- SpatialPixelsDataFrame(out.df[, 1:2], as.data.frame(out.df[, 4:9]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    out.stack <- stack(out.df) 
    
    if (i == loopstr[1]){
      # read crop frac file for weighing
      crop.fracs <- read.table("/media/jan/AddData/Data/Stefan_data/cropland_hurtt_mirca_2000_corrsum.out", head=T)
      crop.fracs[, 1:2] <- crop.fracs[, 1:2] + 0.25
      head(crop.fracs)
      cfracs.sp <- SpatialPixelsDataFrame(crop.fracs[, 1:2], as.data.frame(crop.fracs[, 3:8]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      cfracs.stack <- stack(cfracs.sp)
      frac.mask <- crop(cfracs.stack, out.stack[[1]])
      frac.mask <- mask(frac.mask, out.stack[[1]])
      frac.mask <- frac.mask[[c(2, 1, 3, 5, 4, 6)]]
    }
    
    out.stack <- out.stack * 10 / (0.85 * 2.0 * 0.446) # convert to tonnes per hectar and wet matter
    
    yield.wt.tew <- ((out.stack[[1]] * frac.mask[[1]]) + (out.stack[[2]] * frac.mask[[2]]) + (out.stack[[4]] * frac.mask[[4]]) + (out.stack[[5]] * frac.mask[[5]])) / (sum(frac.mask[[c(1,2,4,5)]]))
    yield.wt.teco <- (out.stack[[3]] * frac.mask[[3]] + out.stack[[6]] * frac.mask[[6]]) / (sum(frac.mask[[c(3, 6)]]))
    
    crop.stack <- stack(yield.wt.tew, yield.wt.teco)
    
    idx.j <- which(file.out == j)
    gcm.list[[j]] <- crop.stack
    
  } # end loop GCMs
  
  # calculate mean of different GCms
  all.gcms <- stack(gcm.list)
  ### Save data
  # Save dfs in list
  idx <- which(loopstr == i)
  
  df.list.wheat[[idx]] <- mean(all.gcms[[seq(1,10, by=2)]]) # equal weighted mean!!!
  df.list.maize[[idx]] <- mean(all.gcms[[seq(2,10, by=2)]]) # equal weighted mean!!!
  
  cat("iteration:", i, " \n")
  flush.console()
  
} # end loop combinations

### Store everything in one data frame
df.out.wheat <- lapply(df.list.wheat, as.data.frame) # convert to df
df.out.maize <- lapply(df.list.maize, as.data.frame) # convert to df
df.out.wheat <- do.call(cbind, df.out.wheat)
df.out.maize <- do.call(cbind, df.out.maize)
col.names <- loopstr
names(df.out.wheat) <- col.names  
names(df.out.maize) <- col.names  

#setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/sensitivity/wheat_yield_2040/")
#getwd()
#save(df.out, file = "df_out.RData")

### ANOVA SENSITIVITY 
# prepare for ANOVA
df.out.wheat <- na.omit(df.out.wheat)
df.out.maize <- na.omit(df.out.maize)

n.list.wheat <- list(); n.list.maize <- list()
co2.list.wheat <- list(); co2.list.maize <- list()
clim.list.wheat <- list(); clim.list.maize <- list()

## WHEAT
n <- dim(df.out.wheat)[1]
for (i in 1:n){
  df.temp <- df.out.wheat[i, ]
  df.temp <- t(as.vector(df.temp))
  df <- data.frame(df.temp[, 1])
  df$N <- c(rep("OFF", 4), rep("ON", 4))
  df$co2 <- rep(c("OFF", "OFF", "ON", "ON"), 2)
  df$clim <- rep(c("OFF", "ON"), 4)
  names(df)[1] <- "site"
  
  # ANOVA variance decomposition
  Fit <- summary(aov(site ~ N * co2 * clim, data=df))
  print(Fit)
  SumSq <- Fit[[1]][, 2]
  Total <- (2 ^ 3 - 1) * var(df$site)
  Indices <- 100 * SumSq / Total
  #TabIndices <- cbind(Fit[[1]], Indices)#[order(Indices, decreasing=T),]
  totalIndex.N <- 100 * ((SumSq[1] + SumSq[4] + SumSq[5] + SumSq[7])  / Total)
  totalIndex.co2 <- 100 * ((SumSq[2] + SumSq[4] + SumSq[6] + SumSq[7])  / Total)
  totalIndex.clim <- 100 * ((SumSq[3] + SumSq[5] + SumSq[6] + SumSq[7])  / Total)
  
  # 1st: N, 2nd: co2, 3rd: clim
  n.list.wheat[[i]] <- c(Indices[1], totalIndex.N)
  co2.list.wheat[[i]] <- c(Indices[2], totalIndex.co2)
  clim.list.wheat[[i]] <- c(Indices[3], totalIndex.clim)
  #n.co2.list.wheat[i] <- Indices[4]
  #n.clim.list.wheat[i] <- Indices[5]
  #co2.clim.list.wheat[i] <- Indices[6]
  #n.co2.clim.list.wheat[i] <- Indices[7]
  
  cat("iteration:", i, " \n")
  flush.console()
}

## MAIZE
n <- dim(df.out.maize)[1]
for (i in 1:n){
  df.temp <- df.out.maize[i, ]
  df.temp <- t(as.vector(df.temp))
  df <- data.frame(df.temp[, 1])
  df$N <- c(rep("OFF", 4), rep("ON", 4))
  df$co2 <- rep(c("OFF", "OFF", "ON", "ON"), 2)
  df$clim <- rep(c("OFF", "ON"), 4)
  names(df)[1] <- "site"
  
  # ANOVA variance decomposition
  Fit <- summary(aov(site ~ N * co2 * clim, data=df))
  print(Fit)
  SumSq <- Fit[[1]][, 2]
  Total <- (2 ^ 3 - 1) * var(df$site)
  Indices <- 100 * SumSq / Total
  TabIndices <- cbind(Fit[[1]], Indices)[order(Indices, decreasing=T),]
  totalIndex.N <- 100 * ((SumSq[1] + SumSq[4] + SumSq[5] + SumSq[7])  / Total)
  totalIndex.co2 <- 100 * ((SumSq[2] + SumSq[4] + SumSq[6] + SumSq[7])  / Total)
  totalIndex.clim <- 100 * ((SumSq[3] + SumSq[5] + SumSq[6] + SumSq[7])  / Total)
  
  # 1st: N, 2nd: co2, 3rd: clim
  n.list.maize[[i]] <- c(Indices[1], totalIndex.N)
  co2.list.maize[[i]] <- c(Indices[2], totalIndex.co2)
  clim.list.maize[[i]] <- c(Indices[3], totalIndex.clim)
  #n.co2.list.maize[i] <- Indices[4]
  #n.clim.list.maize[i]<- Indices[5]
  #co2.clim.list.maize[i]<- Indices[6]
  #n.co2.clim.list.maize[i]<- Indices[7]
  
  cat("iteration:", i, " \n")
  flush.console()
}

## rasters
dummy.r.m <- stack(df.list.maize[[1]], df.list.maize[[1]]); dummy.r.w <- stack(df.list.wheat[[1]], df.list.wheat[[1]])
idx.m <- which(!is.na(dummy.r.m[])); idx.w <- which(!is.na(dummy.r.w[]))
# set all values to NA for safety
dummy.r.w[] <- NA
dummy.r.m[] <- NA
n.r.w <- co2.r.w <- clim.r.w <- test.r.w <- dummy.r.w
n.r.m <- co2.r.m <- clim.r.m <- test.r.m <- dummy.r.m

n.r.w[[1]][idx.w] <- sapply(n.list.wheat, "[", 1) #main effect
n.r.w[[2]][idx.w] <- sapply(n.list.wheat, "[", 2) #total effect
co2.r.w[[1]][idx.w] <- sapply(co2.list.wheat, "[", 1) #main effect
co2.r.w[[2]][idx.w] <- sapply(co2.list.wheat, "[", 2) #total effect
clim.r.w[[1]][idx.w] <- sapply(clim.list.wheat, "[", 1) #main effect
clim.r.w[[2]][idx.w] <- sapply(clim.list.wheat, "[", 2) #total effect

n.r.m[[1]][idx.m] <- sapply(n.list.maize, "[", 1) #main effect
n.r.m[[2]][idx.m] <- sapply(n.list.maize, "[", 2) #total effect
co2.r.m[[1]][idx.m] <- sapply(co2.list.maize, "[", 1) #main effect
co2.r.m[[2]][idx.m] <- sapply(co2.list.maize, "[", 2) #total effect
clim.r.m[[1]][idx.m] <- sapply(clim.list.maize, "[", 1) #main effect
clim.r.m[[2]][idx.m] <- sapply(clim.list.maize, "[", 2) #total effect

## Calculate zonal mean (9 zones)
if (calc.mean == "zonal"){
  setwd("/home/jan/GIS_data/Agroclimatic_zones")
  ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")
  # load zone rasters
  files.r <- list.files(pattern="*.tif$")[1:9]
  zones.stack <- stack(files.r)
  zones.rast <- mean(zones.stack, na.rm=T) 
  plot(zones.rast) 
  
  n.mean.w <- as.data.frame(zonal(n.r.w[[2]], zones.rast, fun='mean'))
  co2.mean.w <- as.data.frame(zonal(co2.r.w[[2]], zones.rast, fun='mean'))
  clim.mean.w <- as.data.frame(zonal(clim.r.w[[2]], zones.rast, fun='mean'))
  n.mean.m <- as.data.frame(zonal(n.r.m[[2]], zones.rast, fun='mean'))
  co2.mean.m <- as.data.frame(zonal(co2.r.m[[2]], zones.rast, fun='mean'))
  clim.mean.m <- as.data.frame(zonal(clim.r.m[[2]], zones.rast, fun='mean'))
}

data.df.mean <- as.data.frame(cbind(n.mean.m, co2.mean.m, clim.mean.m, n.mean.w, co2.mean.w, clim.mean.w))
data.df.mean <- data.df.mean[,-c(1,3,5,7,9,11)]
#data.df.mean <- data.df.mean[-10,]

colnames(data.df.mean)[1:6] <- c("n.m","co2.m","clim.m","n.w","co2.w","clim.w")
#rownames(data.df.mean) <- c("mainEff", "totalEff")
#data.df.sd <- as.data.frame(cbind(n.sd.m, co2.sd.m, clim.sd.m, n.sd.w, co2.sd.w, clim.sd.w))
#rownames(data.df.sd) <- c("mainEff", "totalEff")

# save data
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/3_sensitivity/bugfix_yield_2040/")
saveRDS(data.df.mean, file = paste("df_zmean_5gcms_", years, "_", rcp, ".rds", sep=""))
#saveRDS(data.df.sd, file = paste("df_sd_", years, "_", rcp, ".rds", sep=""))


### calculate ANOVA sensitivity for leaching (save as .rds). do each RCP individually! -----------------------------------------------------------

## Settings
output <- "nflux" 
years <- 2040 # 
rcp <- "rcp85"
calc.mean <- "zonal" #zonal

## Initialize
df.list.leach <- list()

loopstr <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/factorial_design", pattern="clim")

for (i in loopstr){ # loop through 8 combinations 
  
  lsf <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/factorial_design/nOFF_co2OFF_climOFF/output")
  file.out <- grep(paste("nflux", ".out", sep=""), lsf, value=T)
  file.out <- grep(rcp, file.out, value=T)
  
  gcm.list <- list()
  
  for (j in file.out){ # loop through GCM
    
    setwd(paste("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/factorial_design/", i, "/output/", sep=""));
    
    file.to.read <- list.files(pattern = paste(rcp, "_nflux", sep=""))
    out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
    
    ### Convert to spatial objects
    # Subset according to years
    #years.tmp <- years - 1850 + 500
    
    out.df <- subset(out.df, out.df[3] >= min(years) & out.df[3] <= max(years))
    
    out.df <- SpatialPixelsDataFrame(out.df[, 1:2], as.data.frame(out.df[, 8]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    out.stack <- stack(out.df) 
    
    idx.j <- which(file.out == j)
    gcm.list[[j]] <- out.stack
    
  } # end loop GCMs
  
  # calculate mean of different GCms
  all.gcms <- stack(gcm.list)
  ### Save data
  # Save dfs in list
  idx <- which(loopstr == i)
  
  df.list.leach[[idx]] <- mean(all.gcms) # equal weighted mean!!!
  
  cat("iteration:", i, " \n")
  flush.console()
  
} # end loop combinations

### Store everything in one data frame
df.out.leach <- lapply(df.list.leach, as.data.frame) # convert to df
df.out.leach <- do.call(cbind, df.out.leach)
col.names <- loopstr
names(df.out.leach) <- col.names  

#setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/sensitivity/wheat_yield_2040/")
#getwd()
#save(df.out, file = "df_out.RData")

### ANOVA SENSITIVITY 
# prepare for ANOVA
df.out.leach <- na.omit(df.out.leach)

n.list.leach <- list()
co2.list.leach <- list()
clim.list.leach <- list()

## LEACHING
n <- dim(df.out.leach)[1]
for (i in 1:n){
  df.temp <- df.out.leach[i, ]
  df.temp <- t(as.vector(df.temp))
  df <- data.frame(df.temp[, 1])
  df$N <- c(rep("OFF", 4), rep("ON", 4))
  df$co2 <- rep(c("OFF", "OFF", "ON", "ON"), 2)
  df$clim <- rep(c("OFF", "ON"), 4)
  names(df)[1] <- "site"
  
  # ANOVA variance decomposition
  Fit <- summary(aov(site ~ N * co2 * clim, data=df))
  print(Fit)
  SumSq <- Fit[[1]][, 2]
  Total <- (2 ^ 3 - 1) * var(df$site)
  Indices <- 100 * SumSq / Total
  #TabIndices <- cbind(Fit[[1]], Indices)#[order(Indices, decreasing=T),]
  totalIndex.N <- 100 * ((SumSq[1] + SumSq[4] + SumSq[5] + SumSq[7])  / Total)
  totalIndex.co2 <- 100 * ((SumSq[2] + SumSq[4] + SumSq[6] + SumSq[7])  / Total)
  totalIndex.clim <- 100 * ((SumSq[3] + SumSq[5] + SumSq[6] + SumSq[7])  / Total)
  
  # 1st: N, 2nd: co2, 3rd: clim
  n.list.leach[[i]] <- c(Indices[1], totalIndex.N)
  co2.list.leach[[i]] <- c(Indices[2], totalIndex.co2)
  clim.list.leach[[i]] <- c(Indices[3], totalIndex.clim)
  
  cat("iteration:", i, " \n")
  flush.console()
}

## rasters
dummy.r <- stack(df.list.leach[[1]], df.list.leach[[1]])
idx <- which(!is.na(dummy.r[]))
# set all values to NA for safety
dummy.r[] <- NA
n.r <- co2.r<- clim.r <- test.r <- dummy.r

n.r[[1]][idx] <- sapply(n.list.leach, "[", 1) #main effect
n.r[[2]][idx] <- sapply(n.list.leach, "[", 2) #total effect
co2.r[[1]][idx] <- sapply(co2.list.leach, "[", 1) #main effect
co2.r[[2]][idx] <- sapply(co2.list.leach, "[", 2) #total effect
clim.r[[1]][idx] <- sapply(clim.list.leach, "[", 1) #main effect
clim.r[[2]][idx] <- sapply(clim.list.leach, "[", 2) #total effect


## Calculate zonal mean (9 zones)
if (calc.mean == "zonal"){
  setwd("/home/jan/GIS_data/Agroclimatic_zones")
  ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")
  # load zone rasters
  files.r <- list.files(pattern="*.tif$")[1:9]
  zones.stack <- stack(files.r)
  zones.rast <- mean(zones.stack, na.rm=T) 
  plot(zones.rast) 
  
  n.mean <- as.data.frame(zonal(n.r[[2]], zones.rast, fun='mean'))
  co2.mean <- as.data.frame(zonal(co2.r[[2]], zones.rast, fun='mean'))
  clim.mean <- as.data.frame(zonal(clim.r[[2]], zones.rast, fun='mean'))
  n.mean <- as.data.frame(zonal(n.r[[2]], zones.rast, fun='mean'))
  co2.mean <- as.data.frame(zonal(co2.r[[2]], zones.rast, fun='mean'))
  clim.mean <- as.data.frame(zonal(clim.r[[2]], zones.rast, fun='mean'))
}

data.df.mean <- as.data.frame(cbind(n.mean, co2.mean, clim.mean))
data.df.mean <- data.df.mean[,-c(3,5)]

colnames(data.df.mean) <- c("zone","n","co2","clim")

# save data
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/3_sensitivity/bugfix_yield_2040/")
saveRDS(data.df.mean, file = paste("nleach_zmean_5gcms_", years, "_", rcp, ".rds", sep=""))
#saveRDS(data.df.sd, file = paste("df_sd_", years, "_", rcp, ".rds", sep=""))


### bargraph of yield total effect for 2040 and zones (load .rds) -----------------------------------------------------------------------------
setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/3_sensitivity/bugfix_yield_2040/")

mean.effects <- list.files(pattern="df_zmean_5gcms_2040")
meanEff <- list()

for (i in mean.effects) {
 idx <- which(i == mean.effects)
  meanEff[[idx]] <- readRDS(i)
} 

data.df.mean <- meanEff[[1]] # first column rcp45, 2nd column rcp85
#data.df.mean <- data.df.mean[, - c(1)]
data.df.mean[10,] <- colMeans(data.df.mean[1:9,])

yield.df.mean <- data.df.mean
#yield.df.mean <- yield.df.mean[, c(1,3,5,7,)]
data.mean.m <- melt(yield.df.mean)

ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South", "European_Total")
ac.zones.reordered <- c(ac.zones[10], ac.zones[1:9])

data.mean.m$zones <- rep(ac.zones, 6)
data.mean.m$groups <- factor(data.mean.m$zones, levels=ac.zones.reordered)
data.mean.m$crop <- rep(rep(c("Maize", "Wheat"), each=30), 1)
data.mean.m$variable <- rep(rep(c("N", "CO2", "Climate"), each=10), 2)
data.mean.m$rcp <- "8.5"

cols <- c("white", "grey70", "grey30")

(cxc <- ggplot(data.mean.m, aes(x=groups, y=value, fill=variable, width=0.9)) + geom_bar(stat = "identity", size=0.5, col="black", alpha = 1) + theme_bw()  + ggtitle("Driver importance") + ylab("Proportion of sum of squares (%)")  + theme(panel.grid.major = element_line(colour = "white", size = 0.2)) + theme(strip.background = element_rect(fill="grey90")) + facet_grid(~ crop) + scale_fill_manual(values=cols, name="Driver") + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

ggsave(cxc, file="/home/jan/Results/Paper_2_intensity/driver_importance_2040_c12_5gcms_stacked.pdf", width=12, height=8, units="in")



### Plot Nleaching + yield --------------------------------------------------------------------------------------------------

setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/3_sensitivity/bugfix_yield_2040/")

yieldEff <- readRDS("df_zmean_5gcms_2040_rcp85.rds")
leachEff <- readRDS("nleach_zmean_5gcms_2040_rcp85.rds")

data.df.mean <- leachEff
data.df.mean[10,] <- colMeans(data.df.mean)
data.df.mean$zone[10] <- 10
data.df.mean <- data.df.mean[,-1]
leach.df.mean <- data.df.mean

data.df.mean <- yieldEff
data.df.mean[10,] <- colMeans(data.df.mean[1:9, ])
yield.df.mean <- data.df.mean

# combine leaching and yield
dim(leach.df.mean); dim(yield.df.mean)
comb.df.mean <- cbind(yield.df.mean, leach.df.mean)

# melt
data.mean.m <- melt(comb.df.mean)

ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South", "European_Total")
ac.zones.reordered <- c(ac.zones[10], ac.zones[1:9])

data.mean.m$zones <- rep(ac.zones, 9)
data.mean.m$groups <- factor(data.mean.m$zones, levels=ac.zones.reordered)
data.mean.m$variable <- rep(rep(c("N", "CO2", "Climate"), each=10), 3)
data.mean.m$output <- rep(c("Maize", "Wheat", "Leaching"), each=30)
data.mean.m$output <- factor(data.mean.m$output, levels=c("Wheat", "Maize", "Leaching"))

cols <- c("white", "grey70", "grey30")

cxc <- ggplot(data.mean.m, aes(x=groups, y=value, fill=variable, width=0.9)) + 
  geom_bar(stat = "identity", size=0.5, col="black", alpha = 1) + 
  theme_bw() + 
  ggtitle("Driver importance") + 
  ylab("Proportion of sum of squares (%)") + xlab("") +
  theme(panel.grid.major = element_line(colour = "white", size = 0.2)) + 
  theme(strip.background = element_rect(fill = "grey40"), strip.text = element_text(colour = 'white', face="bold")) +
  scale_fill_manual(values=cols, name="Driver") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ output, ncol=1) + 
  theme(legend.position = "top")
cxc

ggsave(cxc, file="/home/jan/Dropbox/Paper_2_nitro/Figures/driver_importance.pdf", width=6, height=10, units="in")
#write.csv(data.mean.m, "/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/3_sensitivity/bugfix_yield_2040/yield_leaching_2040_effects.csv")


### stacked plot of total effects for all years (for appendix) ----------------------------------------------------------------------

setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/3_sensitivity/effects_ts")

mean.effects <- list.files(pattern="mean")
mean.effects.list <- list()
cnt <- 1

for (i in mean.effects){
  mean.effects.list[[cnt]] <- readRDS(i)
  cnt <- cnt + 1 
}

mean.ts.melt <- data.frame(unlist(sapply(mean.effects.list, FUN=function(x) x[2, ]))); colnames(mean.ts.melt) <- "value"

mean.ts.melt$year <- rep(c(2020, 2030, 2040), each=12)
mean.ts.melt$crop <- rep(rep(c("maize", "wheat"), each=3), 6) 
mean.ts.melt$driver <- rep(rep(c("n", "co2", "clim"), 2), 6)
mean.ts.melt$rcp <- rep(rep(c("4.5", "8.5"), 3), each=6)

library(ggplot2)
lap <- ggplot(mean.ts.melt, aes(x=year, y=value, fill=driver)) + geom_area(color='black', size=.2, alpha=.6,  position = 'stack') + theme_bw() + facet_wrap(~ crop) + theme(panel.grid.major = element_line(colour = "white", size = 0.2)) + scale_fill_manual(values=rev(brewer.pal(n=3, 'Blues'))) + facet_grid(rcp ~ crop, scales="free_x")  + theme(panel.grid.major = element_line(colour = "white", size = 0.2))
lap

ggsave(lap, file="/home/jan/Results/Paper_2_intensity/total_effect_ts.pdf", width=8, height=6, units="in")


