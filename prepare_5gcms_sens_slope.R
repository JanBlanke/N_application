#****************************************************
#* Name: prepare_5gcms_sens_slope.R                 *
#* Author: Jan Blanke                               *
#* Description: Process LPJ-GUESS sims for analysis *
#* Dependecies: standalone, factors simulations     *
#* Comments: Run each RCP individually!             *
#****************************************************

library(sp)
library(raster)
library(rasterVis)
library(plyr)

### Settings - change accordingly 
output <- "leaching" # "yield", "leaching"    


### Process wheat and maize yield -----------------------------------------------------------------------------------------------
if (output == "yield"){

  crops <- c("TeWW", "TeCo") # TeCo, TeWW
  years.ts <- seq(2008, 2040, by=1) 
  years.avg <- 0 # year +- 1 
  firstyear <- 1850
  output <- "yield" # lai, cpool, cflux, cmass
  years <- 2040 # 
  rcp <- "rcp85"
  
  ## Initialize
  df.list.wheat <- list()
  df.list.maize <- list()
  
  loopstr <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/")[1:4]
  
  for (i in loopstr){ # loop through 4 combinations 
    
    lsf <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/allDynamic/output/")
    file.out <- grep(paste("yield", ".out", sep=""), lsf, value=T)
    file.out <- grep(rcp, file.out, value=T)
    
    gcm.list <- list()
    
    for (j in file.out){ # loop through GCM
      
      setwd(paste("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/", i, "/output/", sep=""));
      
      file.to.read <- list.files(pattern = paste(rcp, "_yield", sep=""))
      
      out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
      
      ## Convert to spatial objects
      
      time.list <- list()
      for (k in years.ts){
        years.tmp <- k 
        out.df.sub <- subset(out.df, out.df[3] >= min(years.tmp) & out.df[3] <= max(years.tmp))
        
        out.df.sub <- SpatialPixelsDataFrame(out.df.sub[, 1:2], as.data.frame(out.df.sub[, 4:9]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        out.stack <- stack(out.df.sub) 
        
        if (i == loopstr[1]){
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
        
        out.stack <- out.stack * 10 / (0.85 * 2.0 * 0.446) # convert to tonnes per hectar and wet matter
        
        yield.wt.tew <- ((out.stack[[1]] * frac.mask[[1]]) + (out.stack[[2]] * frac.mask[[2]]) + (out.stack[[3]] * frac.mask[[3]]) + (out.stack[[4]] * frac.mask[[4]])) / (sum(frac.mask[[1:4]]))
        yield.wt.teco <- (out.stack[[5]] * frac.mask[[5]] + out.stack[[6]] * frac.mask[[6]]) / (sum(frac.mask[[5:6]]))
        
        crop.stack <- stack(yield.wt.tew, yield.wt.teco)
        idx.k <- which(years.ts == k)
        time.list[[idx.k]] <- crop.stack
        
      }
      
      
      all.years <- stack(time.list)
      idx.j <- which(file.out == j)
      gcm.list[[j]] <- all.years
      
      
    } # end loop GCMs
    
    # calculate mean of different GCms
    all.gcms <- stack(gcm.list)
    
    all.wheat <- all.gcms[[seq(1, 330, by=2)]]
    all.maize <- all.gcms[[seq(2, 330, by=2)]]
    
    # average over gcms with equal weights
    all.wheat <- (all.wheat[[1:length(years.ts)]] + all.wheat[[(length(years.ts) + 1):(length(years.ts) * 2)]] + all.wheat[[(1 + (length(years.ts) * 2)):(length(years.ts) * 3)]] + all.wheat[[(1 + (length(years.ts) * 3)):(length(years.ts) * 4)]] + all.wheat[[(1 + (length(years.ts) * 4)):(length(years.ts) * 5)]]) / 5
    
    all.maize <- (all.maize[[1:length(years.ts)]] + all.maize[[(length(years.ts) + 1):(length(years.ts) * 2)]] + all.maize[[(1 + (length(years.ts) * 2)):(length(years.ts) * 3)]] + all.maize[[(1 + (length(years.ts) * 3)):(length(years.ts) * 4)]] + all.maize[[(1 + (length(years.ts) * 4)):(length(years.ts) * 5)]]) / 5
    
    ## Save data
    # Save dfs in list
    idx <- which(loopstr == i)
    
    df.list.wheat[[idx]] <- all.wheat
    df.list.maize[[idx]] <- all.maize
    
    cat("iteration:", i, " \n")
    flush.console()
    
    # write out
    setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/5GCMs_netcdf/bugfix_33y")
    
    file.name.wheat <- paste("yield_wheat_5gcms_", rcp, "_", i, ".nc", sep="") 
    file.name.maize <- paste("yield_maize_5gcms_", rcp, "_", i, ".nc", sep="") 
    
    writeRaster(all.wheat, file=file.name.wheat)
    writeRaster(all.maize, file=file.name.maize)
    
  } # end loop combinations

  
###### Process N leaching ----------------------------------------------------------------------------------------------------------
} else if (output == "leaching"){
  
  crops <- c("leach") 
  years.ts <- seq(2008, 2040, by=1)  #c(2008, 2020, 2030, 2040)
  years.avg <- 0 # year +- 1 
  firstyear <- 1850
  output <- "nflux" # lai, cpool, cflux, cmass
  years <- 2040 # 
  rcp <- "rcp85"
  
  
  ## Initialize
  df.list.leach <- list()
  
  loopstr <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/")[1:4]
  
  for (i in loopstr){ # loop through 8 combinations 
    
    lsf <- list.files("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/allDynamic/output/")
    file.out <- grep(paste(output, ".out", sep=""), lsf, value=T)
    file.out <- grep(rcp, file.out, value=T)
    
    gcm.list <- list()
    
    for (j in file.out){ # loop through GCM
      
      setwd(paste("/home/jan/LPJ_GUESS/Paper_2/bugfix_run/ensemble/", i, "/output/", sep=""));
      
      file.to.read <- list.files(pattern = paste(rcp, "_nflux", sep=""))
      
      out.df <- read.table(j, header=T); out.df[, 1:2] <- out.df[, 1:2] + 0.25 
      
      ## Convert to spatial objects
      # Subset according to years
      
      time.list <- list()
      for (k in years.ts){
        years.tmp <- k 
        
        out.df.sub <- subset(out.df, out.df[3] >= min(years.tmp) & out.df[3] <= max(years.tmp))
        
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
    all.gcms <- stack(gcm.list)
    all.leach <- all.gcms
    
    # average over gcms with equal weights
    all.leach <- (all.leach[[1:length(years.ts)]] + all.leach[[(length(years.ts) + 1):(length(years.ts) * 2)]] + all.leach[[(1 + (length(years.ts) * 2)):(length(years.ts) * 3)]] + all.leach[[(1 + (length(years.ts) * 3)):(length(years.ts) * 4)]] + all.leach[[(1 + (length(years.ts) * 4)):(length(years.ts) * 5)]]) / 5
    
    ## Save data
    # Save dfs in list
    idx <- which(loopstr == i)
    
    df.list.leach[[idx]] <- all.leach
    
    
    cat("iteration:", i, " \n")
    flush.console()
    
    # write out
    setwd("/media/jan/AddData/Simulations_processed/Paper_2_nitrogen/4_sens_slope/5GCMs_netcdf/bugfix_33y")
    file.name.leach <- paste("leach_5gcms_", rcp, "_", i, ".nc", sep="") 
    writeRaster(all.leach, file=file.name.leach)
    
  }
  
} # end leaching



