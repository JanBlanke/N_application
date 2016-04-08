# *******************************
# title: "Agro-climatic zone map"
# author: "Jan H. Blanke"
# date: "February 20, 2016"
# ******************************

library(rasterVis)
library(rgdal)
library(raster)
library(rasterVis)

eu.countries.mask <- readOGR("/home/jan/GIS_data", "EU_match_clue")
eu.countries.underlay <- readOGR("/home/jan/GIS_data", "world_borders")

ac.zones <- c("Alpine", "Boreal2", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")
setwd("/home/jan/GIS_data/Agroclimatic_zones")
alpine <- readOGR(".", ac.zones[1])
boreal <- readOGR(".", ac.zones[2])
contnorth <- readOGR(".", ac.zones[3])
contsouth <- readOGR(".", ac.zones[4])
medsouth <- readOGR(".", ac.zones[5])
mednorth <- readOGR(".", ac.zones[6])
atcen <- readOGR(".", ac.zones[7])
atnorth <- readOGR(".", ac.zones[8])
atsouth <- readOGR(".", ac.zones[9])

zone.labels <- ac.zones; zone.labels[2] <- "Boreal"

# load zone rasters
files.r <- list.files(pattern="*.tif$")[1:9]
zones.stack <- stack(files.r)
zones.rast <- mean(zones.stack, na.rm=T) 
plot(zones.rast) 

idx <- which(zones.rast[] == 7.5)
zones.rast[idx] <- 6

zones.rast <- ratify(zones.rast)
rat <- levels(zones.rast)[[1]]
rat$zones <- zone.labels
levels(zones.rast) <- rat
zones.rast 

# spcify colors (same as P1 ANOVA results)
#cols <- brewer.pal(9,'Pastel1'); cols[9] <- "white"
cols.tab <- (tableau_color_pal(palette = "tableau10medium")(9))
cols.tab[7] <- "white"
cols.tab[3] <- "#74c476"
cols.tab[5] <- "#807dba"
cols.tab[9] <- "#66c2a5" # brown

cols <- cols.tab[c(8,2,4,9,5,6,1,7,3)]


specTheme <- rasterTheme(region=cols)

png("/home/jan/Results/Paper_2_intensity/zones_raster.png", width = 15 , height = 22, units = 'cm', res = 600)
levelplot(zones.rast, par.settings=specTheme) + layer(sp.polygons(alpine, lwd=.5)) + layer(sp.polygons(boreal, lwd=.25)) + layer(sp.polygons(contnorth, lwd=.25)) + layer(sp.polygons(contsouth, lwd=.25)) + layer(sp.polygons(medsouth, lwd=.25)) + layer(sp.polygons(mednorth, lwd=.25)) + layer(sp.polygons(atcen, lwd=.25)) + layer(sp.polygons(atnorth, lwd=.25)) + layer(sp.polygons(atsouth, lwd=.25)) + layer(sp.polygons(eu.countries.mask, lwd=1, alpha=1, col="black")) + layer(sp.polygons(eu.countries.underlay, lwd=.2, alpha=1, fil="gray90", col="gray55"), under=T)
dev.off()

