###################################
### Plot CAPRI data, e.g. trend ###
###################################

library(rgdal)
library(maptools)
library(rasterVis)
library(gdata)

## Read NUTS2 shapefile
eu.nuts <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/", "NUTS_ETRS_1989_LAEA")
eu.nuts <- spTransform(eu.nuts, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
eu.nuts <- eu.nuts[!is.na(eu.nuts$CAPRI_NUTS), ]

#alternative?
#eu.2010 <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/NUTS_2_shapefiles/", "NUTS2_eurostat_2010")
#eu.nuts <- eu.2010[!is.na(eu.2010$NUTS_ID), ]
#names(eu.nuts)[1] <- "CAPRI_NUTS"

eu.ext <- readOGR("/home/jan/GIS_data", "world_borders")
eu.ext <- crop(eu.ext, extent(-20, 40, 30, 70))
map_nuts2 <- readOGR("/home/jan/Dropbox/Paper_2_nitro/Data/", "Europe_country_borders")

source("/home/jan/Dropbox/R/functions/color_ramps.R")

## Read nitrogen data from Julia
load("/home/jan/Dropbox/Paper_2_nitro/Data/CAPRI/ToClueA2_no_header_sheet2.Rdata")
nappl.a2 <- nappl
load("/home/jan/Dropbox/Paper_2_nitro/Data/CAPRI/ToClue_no_header_sheet3.Rdata")
nappl.hist <- nappl; rm(nappl)
#load("/home/jan/Dropbox/Paper_2_nitro/Data/CAPRI/ToClueB1_no_header_sheet2.Rdata")
#nappl.b1 <- nappl

colnames(nappl.a2)[1:2] <- c("CAPRI_NUTS", "Year") 
#colnames(nappl.hist)[1:2] <- c("CAPRI_NUTS", "Year")
#colnames(nappl.b1)[1:2] <- c("CAPRI_NUTS", "Year")

## Relative change maps
#wheat.list <- list()
#maize.list <- list()

# Scenario subsets
#years <- unique(nappl.a2$Year)

#for (i in years){
  #i <- years[[1]] 

nappl.a2 <- subset(nappl.a2, nappl.a2$Year == "A2_2040Y1")

nappl.to.merge <- nappl.a2[, c(1:2, 10, 15)]
merged <- merge(eu.nuts, nappl.to.merge, by="CAPRI_NUTS")
  
 # idx <- which(i == years)
 # wheat.list[[idx]] <- merged[, "SWHE"]
 # maize.list[[idx]] <- merged[, "MAIZ"] 
#}


# temp test #
#spplot(merged[, "SWHE"])
# temp test #

#wheat.list.df <- lapply(wheat.list, as.data.frame)
#maize.list.df <- lapply(maize.list, as.data.frame)
#wheat.df.scen <- data.frame(matrix(unlist(wheat.list.df), nrow=253, byrow=T), stringsAsFactors=FALSE)
#maize.df.scen <- data.frame(matrix(unlist(maize.list.df), nrow=253, byrow=T), stringsAsFactors=FALSE)

# Historical subset
nappl.sub <- subset(nappl.hist, nappl.hist$Year == 2008)[, 1:10]  
# join NUTS2 polygons and nitrogen data
merged <- merge(eu.nuts, nappl.sub, by="CAPRI_NUTS")
wheat.df.2008 <- as.data.frame(merged[, "SWHE"])
maize.df.2008 <- as.data.frame(merged[, "MAIZ"])

# Calculate absolute change for wheat
wheat.ac.2020 <- wheat.df.scen[, 1] - wheat.df.2008
wheat.ac.2030 <- wheat.df.scen[, 2] - wheat.df.2008
wheat.ac.2040 <- wheat.df.scen[, 3] - wheat.df.2008

library(matrixStats)
wheat.ac.all <- as.data.frame(cbind(wheat.ac.2020, wheat.ac.2030, wheat.ac.2040))
wheat.ac.mean <- rowMeans(wheat.ac.all, na.rm=T)
wheat.ac.median <- rowMedians(as.matrix(wheat.ac.all), na.rm=T)
merged$Wheat <- wheat.ac.median

# Calculate absolute change for maize
maize.ac.2020 <- maize.df.scen[, 1] - maize.df.2008
maize.ac.2030 <- maize.df.scen[, 2] - maize.df.2008
maize.ac.2040 <- maize.df.scen[, 3] - maize.df.2008

maize.ac.all <- as.data.frame(cbind(maize.ac.2020, maize.ac.2030, maize.ac.2040))
maize.ac.mean <- rowMeans(maize.ac.all, na.rm=T)
maize.ac.median <- rowMedians(as.matrix(maize.ac.all), na.rm=T)
merged$Maize<- maize.ac.median

rng.w <- max(range(as.data.frame(merged["Wheat"]), na.rm = T))
rng.m <- max(range(as.data.frame(merged["Maize"]), na.rm = T))


#--- Make maps for paper ------------------------------------------------------------------------------------------------------------

### plot actual nmin
names(merged)[c(5, 6)] <- c("Wheat", "Maize")

cols <- c(viridis(7, alpha = 0.9, begin = 0.0, end = 0.95, option = "viridis"))
cols <- colorRampPalette(c(brewer.pal(9,"RdYlGn")[c(3,4,5,6,7,8,9)], "darkgreen"))(7)
cols <- c(brewer.pal(7,"Reds"))

library(classInt)
library(RColorBrewer)
#pal = viridis(7, alpha = 0.9, begin = 0.0, end = 0.95, option = "viridis")


brks.qt = classIntervals(merged$Maize, n = 7, style = "quantile")
brks.jk = classIntervals(merged$Maize, n = 7, style = "jenks"); brks.jk$brks[8] <- brks.jk$brks[8] + 1

my.settings <- list(strip.background=list(col="grey40"), background=list(col="white"))

png("/home/jan/Dropbox/Paper_2_nitro/Figures/nappl_a2_2008_viridis.png", width = 20, height = 13, units = 'cm', res = 600)
spplot(merged[, c("Wheat", "Maize")], 
       col.regions = cols, 
       at=brks.jk$brks,
       col = "black", 
       main=NULL, 
       sp.layout=list(eu.ext, fill="grey90", col="black", lwd=0.25), 
       colorkey = list(space = "bottom", width = 2, height = 1),  
       xlim = c(-11, 32), ylim = c(34, 70), 
       lwd=0.07, 
       #at=seq(0,669, length=15), 
       do.log=T,
       par.settings=my.settings,
       par.strip.text = list(col = "white", font=2)
)
dev.off()


### plot mean 
#library(Cairo)
#CairoPDF(file="/home/jan/Results/Paper_2_intensity/nappl_a2_mean_abs_change.pdf", width=20/2.54, height=13/2.54, family="Helvetica", pointsize=11)
#pdf("/home/jan/Results/Paper_2_intensity/nappl_a2_mean_abs_change.pdf", width=20/2.54 , height=13/2.54)
#png("/home/jan/Results/Paper_2_intensity/nappl_a2_mean_abs_change.png", width = 18, height = 12, units = 'cm', res = 600)
spplot(merged[, c("Wheat", "Maize")], col.regions=rev(rdbuPal), col = "grey25", main=NULL, sp.layout=list(eu.ext, fill="grey90", col="grey55", lwd=0.01), at=seq(-rng.m, rng.m, length=20), colorkey = list(space = "bottom", width = 2, height = 1),  xlim = c(-11, 32), ylim = c(34, 70), lwd=0.3, par.settings = my.settings)
#dev.off()












