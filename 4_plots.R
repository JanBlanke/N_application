#######################################
### Name: 4_plots.R                 ###
### Author: Jan Blanke              ###
### Description: Analysis and plots ### 
#######################################

library(raster)
library(reshape2)
library(ggplot2)
library(scales)
library(ggthemes)

source("/home/jan/Dropbox/R/functions/colorramps.R")
eu.countries.mask <- readOGR("/home/jan/GIS_data", "EU_match_clue")

setwd("/home/jan/Results/Paper_2_intensity/yield_rasters")


### Load raster files -------------------------------------------------------
list.files(pattern="mean_allDynamic_rcp85")

## Read as rasters
# all dynamic WHEAT
ref.stack.45.wheat <- stack("TeWW_yield_ensemble_mean_allDynamic_rcp45.nc")
ref.stack.85.wheat <- stack("TeWW_yield_ensemble_mean_allDynamic_rcp85.nc")
ref.stack.45.corn <- stack("TeCo_yield_ensemble_mean_allDynamic_rcp45.nc")
ref.stack.85.corn <- stack("TeCo_yield_ensemble_mean_allDynamic_rcp85.nc")
# N constant
Noff.stack.45.wheat <- stack("TeWW_yield_ensemble_mean_constantN_rcp45.nc")
Noff.stack.85.wheat <- stack("TeWW_yield_ensemble_mean_constantN_rcp85.nc")
Noff.stack.45.corn <- stack("TeCo_yield_ensemble_mean_constantN_rcp45.nc")
Noff.stack.85.corn <- stack("TeCo_yield_ensemble_mean_constantN_rcp85.nc")
# CO2 constant
CO2off.stack.45.wheat <- stack("TeWW_yield_ensemble_mean_constantCO2_rcp45.nc")
CO2off.stack.85.wheat <- stack("TeWW_yield_ensemble_mean_constantCO2_rcp85.nc")
CO2off.stack.45.corn <- stack("TeCo_yield_ensemble_mean_constantCO2_rcp45.nc")
CO2off.stack.85.corn <- stack("TeCo_yield_ensemble_mean_constantCO2_rcp85.nc")
# N and CO2 constant
CO2Noff.stack.45.wheat <- stack("TeWW_yield_ensemble_mean_constantCO2N_rcp45.nc")
CO2Noff.stack.85.wheat <- stack("TeWW_yield_ensemble_mean_constantCO2N_rcp85.nc")
CO2Noff.stack.45.corn <- stack("TeCo_yield_ensemble_mean_constantCO2N_rcp45.nc")
CO2Noff.stack.85.corn <- stack("TeCo_yield_ensemble_mean_constantCO2N_rcp85.nc")

## Read as data frame
# all dynamic WHEAT
ref.df.45.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_allDynamic_rcp45.nc"))
ref.df.85.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_allDynamic_rcp85.nc"))
ref.df.45.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_allDynamic_rcp45.nc"))
ref.df.85.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_allDynamic_rcp85.nc"))
# N constant
Noff.df.45.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_constantN_rcp45.nc"))
Noff.df.85.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_constantN_rcp85.nc"))
Noff.df.45.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_constantN_rcp45.nc"))
Noff.df.85.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_constantN_rcp85.nc"))
# CO2 constant
CO2off.df.45.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_constantCO2_rcp45.nc"))
CO2off.df.85.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_constantCO2_rcp85.nc"))
CO2off.df.45.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_constantCO2_rcp45.nc"))
CO2off.df.85.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_constantCO2_rcp85.nc"))
# N and CO2 constant
CO2Noff.df.45.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_constantCO2N_rcp45.nc"))
CO2Noff.df.85.wheat <- as.data.frame(stack("TeWW_yield_ensemble_mean_constantCO2N_rcp85.nc"))
CO2Noff.df.45.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_constantCO2N_rcp45.nc"))
CO2Noff.df.85.corn <- as.data.frame(stack("TeCo_yield_ensemble_mean_constantCO2N_rcp85.nc"))

all.df <- cbind(ref.df.45.wheat, ref.df.85.wheat, ref.df.45.corn, ref.df.85.corn, Noff.df.45.wheat, Noff.df.85.wheat, Noff.df.45.corn, Noff.df.85.corn, CO2off.df.45.wheat, CO2off.df.85.wheat, CO2off.df.45.corn, CO2off.df.85.corn, CO2Noff.df.45.wheat, CO2Noff.df.85.wheat, CO2Noff.df.45.corn, CO2Noff.df.85.corn)
all.df <- na.omit(all.df) # nrow = 1918

all.names <- c("ref.df.45.wheat", "ref.df.85.wheat", "ref.df.45.corn", "ref.df.85.corn", "Noff.df.45.wheat", "Noff.df.85.wheat", "Noff.df.45.corn", "Noff.df.85.corn", "CO2off.df.45.wheat", "CO2off.df.85.wheat", "CO2off.df.45.corn", "CO2off.df.85.corn", "CO2Noff.df.45.wheat", "CO2Noff.df.85.wheat", "CO2Noff.df.45.corn", "CO2Noff.df.85.corn")

all.names <- rep(all.names, each=4)
for (i in seq(1, 64, 4)) all.names[i] <- paste(all.names[i], "_2008", sep="")
for (i in seq(2, 64, 4)) all.names[i] <- paste(all.names[i], "_2020", sep="")
for (i in seq(3, 64, 4)) all.names[i] <- paste(all.names[i], "_2030", sep="")
for (i in seq(4, 64, 4)) all.names[i] <- paste(all.names[i], "_2040", sep="")
colnames(all.df) <- all.names


### 1 Boxplots --------------------------------------------------------------
# seperate dataset again into wheat and maize
wheat.df <- all.df[, grep("wheat", all.names)]
corn.df <- all.df[, -grep("wheat", all.names)]

# wheat
melt.df <- reshape2::melt(corn.df, id.vars = NULL) 
melt.df$year <- rep(rep(c("2008", "2020", "2030", "2040"), each=1914), times = 8)
melt.df$rcp <- rep(rep(rep(c("4.5", "8.5"), each=1914), each=4), times=4)
 
# change order
#melt.df$variable <- factor(melt.df$variable, levels = c('ref.df.45.wheat_2008','refrast.85', 'co2off.45', 'co2off.85', 'noff.45', 'noff.85'), ordered = TRUE)

ggplot(melt.df, aes(x=variable, y=value)) + geom_violin(aes(fill=rcp)) + geom_boxplot(width=0.3) + theme_bw(15) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap( rcp ~ year, scales="free_x") + stat_summary(fun.y=mean, colour="black", geom="point", size=3)
#ggsave(bxpl, file="/home/jan/Results/Paper_2_intensity/boxplot_wheat_2040.pdf", width=11, height=6, units="in")


### 2 Correlation bar plots ---------------------------------------------------
corr1 <- cor(rast.df[, 2], rast.df[, 4], method="pearson")
corr2 <- cor(rast.df[, 1], rast.df[, 3], method="pearson")
corr3 <- cor(rast.df[, 2], rast.df[, 6], method="pearson")
corr4 <- cor(rast.df[, 1], rast.df[, 5], method="pearson")

corr5 <- cor(rast.df[, 8], rast.df[, 10], method="pearson")
corr6 <- cor(rast.df[, 7], rast.df[, 9], method="pearson")
corr7 <- cor(rast.df[, 8], rast.df[, 12], method="pearson")
corr8 <- cor(rast.df[, 7], rast.df[, 11], method="pearson")

corr.df <- data.frame(variable = c("N.wheat.85", "N.wheat.45", "CO2.wheat.85", "CO2.wheat.45", "N.maize.85", "N.maize.45", "CO2.maize.85", "CO2.maize.45"), value = c(corr1, corr2, corr3, corr4, corr5, corr6, corr7, corr8))
corr.df$group <- rep(c("wheat", "maize"), each=4)

cols <- brewer.pal(8,"Paired")
#(tableau_color_pal(palette = "tableau20")(4))

scale.labs <- data.frame(x = rep(1, 10), y = seq(0.1, 1, by=0.1),labels = as.character(seq(0.1, 1, by=0.1)))

(crrplt <- ggplot(corr.df, aes(x=factor(variable), value, fill=variable)) + geom_bar(stat="identity",  colour = "black", width = 1, alpha = 0.6)  + theme_bw() + ggtitle("correlation for 2040")  + xlab("factor") + ylab("pearson correlation coefficient") + scale_fill_manual(values=cols, name= "Output") + coord_cartesian(ylim=c(0.8,1.1))) 

ggsave(crrplt, file="/home/jan/Results/Paper_2_intensity/corrplot_2040.pdf", width=11, height=6, units="in")


### 3 Plot yield predictions (raster maps) ---------------------------------------------------
library(rgdal)
library(rasterVis)
eu.countries <- readOGR("/home/jan/GIS_data", "EU_27_final")
eu.countries.mask <- readOGR("/home/jan/GIS_data", "EU_match_clue")
eu.countries.underlay <- readOGR("/home/jan/GIS_data", "EU_27_underlay")

wheat.stack <- stack(ref.stack.45.wheat[[1]], ref.stack.45.wheat[[4]], ref.stack.85.wheat[[4]])
names(wheat.stack) <- c("wheat.2008", "wheat.2040.rcp45", "wheat.2040.rcp85")

corn.stack <- stack(ref.stack.45.corn[[1]], ref.stack.45.corn[[4]], ref.stack.85.corn[[4]])
names(corn.stack) <- c("maize.2008", "maize.2040.rcp45", "maize.2040.rcp85")

modTheme <- modifyList(RdYlGnTheme, list(layout.heights=list(top.padding=0, bottom.padding=0), strip.background = list(col = 'gray85'), strip.border = list(col = 'black', lwd=0.5), strip.text=list(cex=0.5)))

maxval.wheat <- max(range(wheat.stack)@data@max)
(wheat.plot <- levelplot(wheat.stack, par.settings=modTheme, main="Wheat yield predictions", at=seq(0, maxval.wheat, length=20), layout=c(3,1), scales=list(y=list(draw=1), x=list(draw=0)), xlab=NULL) + layer(sp.polygons(eu.countries.underlay, lwd=.3, alpha=1, fil="gray90", col="gray55"), under=T) + layer(sp.polygons(eu.countries.mask, lwd=1, alpha=1, col="black")))

maxval.corn <- max(range(corn.stack)@data@max)
(corn.plot <- levelplot(corn.stack, par.settings=modTheme, main="Maize yield predictions", at=seq(0, maxval.corn, length=20), layout=c(3,1), scales=list(y=list(draw=1), x=list(draw=1))) + layer(sp.polygons(eu.countries.underlay, lwd=.3, alpha=1, fil="gray90", col="gray55"), under=T) + layer(sp.polygons(eu.countries.mask, lwd=1, alpha=1, col="black")))

pw = list(x=5.7, units="cm") 
pl = list(x=7, units="cm") 
png("/home/jan/Results/Paper_2_intensity/yield_predictions.png", width = 27 , height = 21  , units = 'cm', res = 600)
print(wheat.plot, split=c(1,1,1,2), panel.width=pw, panel.height=pl,  more=T) 
print(corn.plot, split=c(1,2,1,2), panel.width=pw, panel.height=pl) 
dev.off()


### 4 Difference rasters ---------------------------------------------------
# yield
alldyn.diff.85 <- ref.stack.85.wheat[[4]] - ref.stack.85.wheat[[1]]
constN.diff.85 <- Noff.stack.85.wheat[[4]] - Noff.stack.85.wheat[[1]]
alldyn.diff.45 <- ref.stack.45.wheat[[4]] - ref.stack.45.wheat[[1]]
constN.diff.45 <- Noff.stack.45.wheat[[4]] - Noff.stack.45.wheat[[1]]

## nleaching
nleach85 <- stack("/home/jan/Results/Paper_2_intensity/nleach_rasters/TeWW_nflux_ensemble_mean_allDynamic_rcp45.nc")
nleach45 <- stack("/home/jan/Results/Paper_2_intensity/nleach_rasters/TeWW_nflux_ensemble_mean_allDynamic_rcp85.nc")
nleach.alldyn.diff.85 <- nleach85[[4]] - nleach85[[1]]
nleach.alldyn.diff.45 <- nleach45[[4]] - nleach45[[1]]
diff.stack <- stack(nleach.alldyn.diff.45, nleach.alldyn.diff.85); names(diff.stack) <- c("rcp45", "rcp85")

# n leaching alldyn
divTheme <- rasterTheme(region=divRampFun(diff.stack))
modTheme <- modifyList(divTheme, list(layout.heights=list(top.padding=0, bottom.padding=-1), strip.background = list(col = 'gray90'), strip.border = list(col = 'black'), strip.text=list(cex=0.5)))

png("/home/jan/Results/Paper_2_intensity/nleach_diff.png", width = 27 , height = 20  , units = 'cm', res = 600)
levelplot(diff.stack, par.settings=modTheme, contour=T) + layer(sp.polygons(eu.countries.underlay, lwd=.3, alpha=1, fil="gray90", col="gray55"), under=T) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black")) 
dev.off()

## alldynamic
# yield alldyn 85
divTheme85 <- rasterTheme(region=divRampFun(alldyn.diff.85))
levelplot(alldyn.diff.85, par.settings=divThemeOpt, at=seq(min(alldyn.diff.85[], na.rm = T), max(alldyn.diff.85[], na.rm = T), length=100)) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black")) 

# yield alldyn 45
divTheme45 <- rasterTheme(region=divRampFun(alldyn.diff.45))
levelplot(alldyn.diff.45, par.settings=divThemeOpt, at=seq(min(alldyn.diff.45[], na.rm = T), max(alldyn.diff.45[], na.rm = T), length=100)) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black")) 

## N constant
# yield constN 85
divTheme85 <- rasterTheme(region=divRampFun(constN.diff.85))
levelplot(constN.diff.85, par.settings=divTheme85, at=seq(min(constN.diff.85[], na.rm = T), max(constN.diff.85[], na.rm = T), length=100)) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black")) 

# yield constN 45
divTheme45 <- rasterTheme(region=divRampFun(constN.diff.45))
levelplot(constN.diff.45, par.settings=divTheme45, at=seq(min(constN.diff.45[], na.rm = T), max(constN.diff.45[], na.rm = T), length=100)) + layer(sp.polygons(eu.countries.mask, lwd=1.5, alpha=1, col="black")) 


### 5 Line graphs of yield based on zonal statistics ---------------------------------------------------
setwd("/home/jan/GIS_data/Agroclimatic_zones")

ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")

# load zone rasters
files.r <- list.files(pattern="*.tif$")[1:9]
zones.stack <- stack(files.r)
zones.rast <- mean(zones.stack, na.rm=T) 
plot(zones.rast) 

setwd("/home/jan/Results/Paper_2_intensity/yield_rasters//") 
list.files()
wheat85 <- stack("TeWW_yield_ensemble_mean_allDynamic_rcp85.nc")
wheat45 <- stack("TeWW_yield_ensemble_mean_allDynamic_rcp45.nc")
#corn85 <- stack("TeCo_yield_ensemble_mean_allDynamic_rcp85.nc")
#corn45 <- stack("TeCo_yield_ensemble_mean_allDynamic_rcp45.nc")
#nleach85 <- stack("/home/jan/Results/Paper_2_intensity/nleach_rasters/TeWW_nflux_ensemble_mean_allDynamic_rcp45.nc")
#nleach45 <- stack("/home/jan/Results/Paper_2_intensity/nleach_rasters/TeWW_nflux_ensemble_mean_allDynamic_rcp85.nc")

crop45 <- wheat45
crop85 <- wheat85
# generalized from here!

# calculate zonal statistics: mean, sd, sem
meanval85 <- as.data.frame(zonal(crop85, zones.rast, fun='mean'))
sdval85 <- zonal(crop85, zones.rast, fun='sd')
count85 <- zonal(crop85, zones.rast, fun='count')
semval85 <- sdval85[, -1] / sqrt(count85[, -1])

meanval45 <- as.data.frame(zonal(crop45, zones.rast, fun='mean'))
sdval45 <- zonal(crop45, zones.rast, fun='sd')
count45 <- zonal(crop45, zones.rast, fun='count')
semval45 <- sdval45[, -1] / sqrt(count45[, -1])

# combine statistics of interest
meanval <- cbind(meanval45, meanval85); meanval <- meanval[,-6]
semval <- cbind(semval45, semval85); colnames(semval) <- c("X1", "X2", "X3", "X4", "X1.2", "X2.2", "X3.2", "X4.2")

# make graphs
tab10 <- tableau_color_pal(palette = "tableau10")(10)
tab10.light <- tableau_color_pal(palette = "tableau10light")(10)

cols <- c(tab10.light[4], tab10[4])
#cols.reds <- c("#fdae61","#D62728")
#cols.blues <- c("#1F77B4")
#colz2 <- brewer.pal(11, "RdBu")[c(2, 4)]
#colz <- (colorRampPalette(c("grey80", "black"))(9)) # for black and white version

mean.melt <- melt(meanval[,-1]) 
sem.melt <- melt(semval) 

mean.melt$year <- rep(rep(c("2008", "2020", "2030", "2040"), each=9), times=2)
mean.melt$zones <- rep(ac.zones, times=8)
mean.melt$rcps <- rep(c("rcp4.5", "rcp8.5"), each=36)
mean.melt$sem <- sem.melt[,3]

# v1: facet wrap for rcps
lnplt <- ggplot(data=mean.melt, aes(x=year, y=value, color=zones)) + geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.07) + geom_line(aes(group=zones), lwd=1) + geom_point(size=2.2, shape=21, fill="white") + theme_bw() + ggtitle("nitrogen leaching for agro-climatic zones") + ylab(expression("nitrogen leached" ~ (kg ~  m^{-2}))) + scale_color_manual(values=colz, name= "Zones") + facet_wrap(~rcps) 

# v2: facet wrap for zones
(lnplt <- ggplot(data=mean.melt, aes(x=year, y=value, color=rcps)) + geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.1) + geom_line(aes(group=rcps), lwd=1) + geom_point(size=2.2, shape=21, fill="white") + theme_bw() + ggtitle("Yield for agro-climatic zones") + ylab("Yield (tonnes/hectare)") + scale_color_manual(values=cols, name= "RCPs") + facet_wrap(~zones) + theme(strip.background = element_rect(fill="grey90")))

#expression("nitrogen leached" ~ (kg ~  m^{-2}))

ggsave(lnplt, file="/home/jan/Results/Paper_2_intensity/lineplot_wheat_zonesplit.pdf", width=12, height=8, units="in")

# same graph as above but add grey lines of other rcp for comparison
dat2 <- mean.melt
dat2$rcps <- dat2$rcps[72:1]

lnplt2 <-ggplot(data=mean.melt, aes(x=year, y=value, color=zones)) + geom_line(data=dat2, aes(x=year, y=value, group=zones), colour="grey80", linetype="solid") + geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.07) + geom_line(aes(group=zones)) + geom_point(size=2.2, shape=21, fill="white") + theme_bw() + ggtitle("nitrogen leaching for agro-climatic zones") + ylab(expression("nitrogen leached" ~ (kg ~  m^{-2}))) + scale_color_manual(values=colz, name= "Zones") + facet_wrap(~rcps) 

ggsave(lnplt2, file="/home/jan/Results/Paper_2_intensity/lineplot_wheat_rcpcomp.pdf", width=12, height=5, units="in")




### 6 Contribution to change based on zonal statistics ---------------------------------------------------
library(raster)
library(reshape2)

setwd("/home/jan/GIS_data/Agroclimatic_zones")
ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")

# load zone rasters
files.r <- list.files(pattern="*.tif$")[1:9]
zones.stack <- stack(files.r)
zones.rast <- mean(zones.stack, na.rm=T) 
plot(zones.rast) 

setwd("/home/jan/Results/Paper_2_intensity/yield_rasters/") 
list.files()
allDyn85 <- stack("TeWW_yield_ensemble_mean_allDynamic_rcp85.nc")
offN85 <- stack("TeWW_yield_ensemble_mean_constantN_rcp85.nc")
offCO285 <- stack("TeWW_yield_ensemble_mean_constantCO2_rcp85.nc")
allDyn45 <- stack("TeWW_yield_ensemble_mean_allDynamic_rcp45.nc")
offN45 <- stack("TeWW_yield_ensemble_mean_constantN_rcp45.nc")
offCO245 <- stack("TeWW_yield_ensemble_mean_constantCO2_rcp45.nc")

mean.tot.85 <- as.data.frame(zonal(allDyn85, zones.rast, fun='mean'))
mean.N.85 <- as.data.frame(zonal(offN85, zones.rast, fun='mean'))
mean.CO2.85 <- as.data.frame(zonal(offCO285, zones.rast, fun='mean'))
mean.tot.45 <- as.data.frame(zonal(allDyn45, zones.rast, fun='mean'))
mean.N.45 <- as.data.frame(zonal(offN45, zones.rast, fun='mean'))
mean.CO2.45 <- as.data.frame(zonal(offCO245, zones.rast, fun='mean'))

delta.tot.85 <- (mean.tot.85[, 5] - mean.tot.85[, 2]) / mean.tot.85[, 2] * 100
delta.N.85 <- (mean.N.85[, 5] - mean.N.85[, 2]) / mean.N.85[, 2] * 100
delta.CO2.85 <- (mean.CO2.85[, 5] - mean.CO2.85[, 2]) / mean.CO2.85[, 2] * 100
delta.tot.45 <- (mean.tot.45[, 5] - mean.tot.45[, 2]) / mean.tot.45[, 2] * 100
delta.N.45 <- (mean.N.45[, 5] - mean.N.45[, 2]) / mean.N.45[, 2] * 100
delta.CO2.45 <- (mean.CO2.45[, 5] - mean.CO2.45[, 2]) / mean.CO2.45[, 2] * 100

df <- data.frame(tot85=delta.tot.85, noff85=delta.tot.85 - delta.N.85, co2off85=delta.tot.85 - delta.CO2.85,
                 tot45=delta.tot.45, noff45=delta.tot.45 - delta.N.45, co2off45=delta.tot.45 - delta.CO2.45)

## stacked bar graphs
# rcp85
barplot(t(df[,1:3]), main="change", xlab="groups", beside=F, args.legend = list(x = "topright", bty = "n"))
# rcp 45
barplot(t(df[,4:6]), main="change", xlab="groups", beside=F, args.legend = list(x = "topright", bty = "n"))

# ggplot2
y.melt <- melt(df)
y.melt$zones <- rep(ac.zones, 6)
y.melt$rcp <- rep(c("rcp85", "rcp45"), each=27)
y.melt$colcode <- rep(rep(c("total", "n", "co2"), each=9),2)

colz <- rev(c("grey20", "grey70", "white"))
cntrbnplt <- ggplot(y.melt, aes(x = zones, y = value, fill=colcode)) + geom_bar(position="stack", stat="identity", col="black")  + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("change (%)") + facet_wrap(~rcp, scales="free_x") + scale_fill_manual(values=colz, name= "Effect")
  
ggsave(cntrbnplt, file="/home/jan/Results/Paper_2_intensity/contribution_wheat.pdf", width=12, height=8, units="in")

# overlapping barplot
cntrbnplt2 <- ggplot(y.melt, aes(x = zones, y = value, fill=colcode)) + 
  geom_bar(position = 'dodge', alpha = .99, stat='identity', col="black") + 
  facet_wrap(~rcp, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("change (%)")  + scale_fill_manual(values=colz, name= "Effect") 
ggsave(cntrbnplt2, file="/home/jan/Results/Paper_2_intensity/contribution_wheat_overlap.pdf", width=12, height=8, units="in")



### 7 Ago-climatic zones map ------------------------------------------------
library(rasterVis)
library(rgdal)

ac.zones <- c("Alpine", "Boreal", "Continental_North", "Continental_South", "Mediterranean_South", "Mediterranean_North", "Atlantic_Central", "Atlantic_North", "Atlantic_South")
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

# load zone rasters
files.r <- list.files(pattern="*.tif$")[1:9]
zones.stack <- stack(files.r)
zones.rast <- mean(zones.stack, na.rm=T) 
plot(zones.rast) 

idx <- which(zones.rast[] == 7.5)
zones.rast[idx] <- 6

zones.rast <- ratify(zones.rast)
rat <- levels(zones.rast)[[1]]
rat$zones <- ac.zones
levels(zones.rast) <- rat
zones.rast 

#specTheme <- rasterTheme(region=sample(tableau_color_pal(palette = "tableau20light")(9)))
cols <- alpha(tableau_color_pal(palette = "tableau20")(9), 0.6)
#specTheme <- rasterTheme(region=(brewer.pal(9,"Set3"))[c(2,7,3,4,5,6,1,9,8)])
specTheme <- rasterTheme(region=cols)

png("/home/jan/Results/Paper_2_intensity/zones_raster.png", width = 20 , height = 27  , units = 'cm', res = 600)
levelplot(zones.rast, par.settings=specTheme) + layer(sp.polygons(alpine, col="black")) + layer(sp.polygons(boreal)) + layer(sp.polygons(contnorth)) + layer(sp.polygons(contsouth)) + layer(sp.polygons(medsouth)) + layer(sp.polygons(mednorth)) + layer(sp.polygons(atcen)) + layer(sp.polygons(atnorth)) + layer(sp.polygons(atsouth)) + layer(sp.polygons(eu.countries.mask, lwd=2, alpha=1, col="black"))
dev.off()



