# load sites geopackage
library(terra)
library(sp)
library(plotrix)

sites <- vect("data/sites.gpkg")
sites
plot(sites)

# county shapefile
cb.dir <- "K:/_gis/_cb"
in.name <- "cb_2019_us_county_500k.shp"
in.path <- paste(cb.dir, in.name, sep="/")
county <- vect(in.path)
county <- county[which(county$STATEFP == "13"),]
plot(county)

in.name <- "cb_2023_us_state_5m.shp"
in.path <- paste(cb.dir, in.name, sep="/")
usa <- vect(in.path)
conus <- c("AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", 
           "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", 
           "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", 
           "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",
           "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", 
           "WV", "WY")
usa <- usa[which(usa$STUSPS %in% conus),]

# get HMI
# see notes at bottom of file for origin
human <- rast("K:/mammal_ga/human_modification_2022_90m.tif")

crs(human)
crs(sites)

# project to match
sitesp <- project(sites, crs(human))
countyp <- project(county, crs(human))

# 1000 m buffer to crop raster for speed
buff1000 <- buffer(sitesp, 1000)

humanc <- crop(human, ext(buff1000))

par(mfrow=c(1,1), mar=c(0.1, 0.1, 0.1, 10.1))
plot(buff1000)
plot(humanc, add=TRUE)
plot(sites, col="red", add=TRUE)

humanc <- round(humanc, 2)
quantile(
    values(humanc),
    probs = c(seq(0, .9, .1), .925, .95, .975, 0.98, .99, 1),
    na.rm = TRUE
)

mean(values(humanc) == 1, na.rm = TRUE)


humanx <- humanc
th <- table(values(humanx))
th <- c(0, 0, th)
names(th)[1:2] <- c("0", "0.01")

plot(th, type="h", lwd=5,
     col=c("black", rep(rainbow(10), each=10)))

brks <- c(
    0, 0.2, 0.4, 0.6, 0.8,
    0.9,
    0.92, 0.94, 0.96, 0.98, 1
)
brks
cols <- hcl.colors(
    length(brks) - 1,
    palette = "YlOrRd"
)
plot(
    humanx,
    breaks = brks,
    col = cols
)


cols <- hcl.colors(
    length(brks) - 1,
    palette = "Inferno"
)
plot(
    humanx,
    breaks = brks,
    col = cols
)

human98 <- humanx >= 0.98
plot(human98, col=c("grey90", "red"), legend=FALSE)

human96 <- humanx >= 0.96
plot(human96, col=c("grey90", "red"), legend=FALSE)

human80 <- humanx >= 0.8
plot(human80, col=c("grey90", "red"), legend=FALSE)

brks <- c(seq(0, 1, by=0.1))
cols <- hcl.colors(length(brks) - 1, palette = "Temps")
plot(humanx, breaks=brks, col=cols, legend=FALSE)

cents <- centroids(sitesp)
use.names <- c("beach", "cb", "clark",
               "disc", "expo", "fern",
               "field", "fraz",  "gate",
               "hb", "liz", "main",
               "mars", "nance", "noon",
               "pass", "pegg", "pickle",
               "pond", "rags", "rot",
               "terr", "zion")

cents <- cents[which(cents$name %in% use.names),]

eusa <- ext(usa)

jpeg("rv-figure-01.jpg", width=14, height=11,
     units="in", res=600)
layout(matrix(c(1,2), nrow=1), widths=c(7,1))
par(mar=c(0, 0, 0, 0), lend=1, las=1, cex.axis=1.3)
plot(buff1000, axes=FALSE, border=NA,
     mar=c(1, 1, 1, 0), legend=FALSE)
plot(humanx, add=TRUE, col=cols, breaks=brks,
     axes=FALSE, legend=FALSE)
plot(countyp, add=TRUE, col="transparent", lwd=3)
points(cents, pch=21, #add=TRUE, 
     cex=3, bg="white", lwd=4)
boxed.labels(-84.864, 34.175, "Bartow County",
             bg="white", cex=2, ypad=1.8)
boxed.labels(-84.627, 33.882, "Cobb County",
             bg="white", cex=2, ypad=1.8)
boxed.labels(-84.520, 33.704, "Fulton County",
             bg="white", cex=2, ypad=1.8)
arrows(-84.411, 34.005, -84.333, 33.920,
       lwd=6, length=0.2, col="white")
boxed.labels(-84.410, 34.010, "DeKalb County",
             bg="white", cex=2, ypad=1.8)
rect(-84.38, 34.04, -84.316, 34.13, #north
     border=NA, col="#FFFFFF80")
rect(-84.566, 34.13,-84.316, 34.18, #sbar
     border=NA, col="#FFFFFF80")
north(xy=c(-84.35, 34.07), type=2,
      cex=3)

sbar(20, xy=c(-84.55, 34.15), type="bar",
     divs=4, below="km", lwd=3, cex=1.5)
##
par(mar=c(0,0,2,4))
plot(NA, xlim=c(0, 1), ylim=c(0, 10),
     type="n", axes=FALSE, xlab="", ylab="")
rect(-0.8, 0:9, 0, 1:10, col=cols,
     border=NA, xpd=NA)
axis(side=4, at=0:10, labels=c(round(brks, 2)),
     las=1, line=-3.5, cex.axis=2)
par(las=0)
mtext("Human Modification index",
      side=4, line=1.5, cex=3)
#arrows(6, 0, 6, 10, xpd=NA,
#       length=0.2, lwd=3, code=3)
#mtext("Less modified", side=4, line=8.5,
#      adj=0.05, cex=2, font=3, xpd=NA)
#mtext("More modified", side=4, line=8.75,
#      adj=0.95, cex=2, font=3, xpd=NA)
par(fig=c(0.09, 0.45, 0.05, 0.33),
    new=TRUE, mar=c(0, 0, 0, 0))
plot(usa, col=ifelse(usa$STUSPS == "GA", "grey", "white"),
     mar=c(0, 0, 0, 0), axes=FALSE, box=TRUE,
     background="white", lwd=2)
lines(ext(buff1000), col="red", lwd=6)
arrows(-98.5, 41.4, -86.5, 34.39, lwd=6, col="red",
       length=0.15)
dev.off()

###########################################


tiff("figures/figure-01-revision.tif", width=14, height=11,
     units="in", res=500, compression="lzw")
layout(matrix(c(1,2), nrow=1), widths=c(7,1))
par(mar=c(0, 0, 0, 0), lend=1, las=1, cex.axis=1.3)
plot(buff1000, axes=FALSE, border=NA,
     mar=c(1, 1, 1, 0), legend=FALSE)
plot(humanx, add=TRUE, col=cols, breaks=brks,
     axes=FALSE, legend=FALSE)
plot(countyp, add=TRUE, col="transparent", lwd=3)
points(cents, pch=21, #add=TRUE, 
       cex=3, bg="white", lwd=4)
boxed.labels(-84.864, 34.175, "Bartow County",
             bg="white", cex=2, ypad=1.8)
boxed.labels(-84.627, 33.882, "Cobb County",
             bg="white", cex=2, ypad=1.8)
boxed.labels(-84.520, 33.704, "Fulton County",
             bg="white", cex=2, ypad=1.8)
arrows(-84.411, 34.005, -84.333, 33.920,
       lwd=6, length=0.2, col="white")
boxed.labels(-84.410, 34.010, "DeKalb County",
             bg="white", cex=2, ypad=1.8)
rect(-84.38, 34.04, -84.316, 34.13, #north
     border=NA, col="#FFFFFF80")
rect(-84.566, 34.13,-84.316, 34.18, #sbar
     border=NA, col="#FFFFFF80")
north(xy=c(-84.35, 34.07), type=2,
      cex=3)

sbar(20, xy=c(-84.55, 34.15), type="bar",
     divs=4, below="km", lwd=3, cex=1.5)
##
par(mar=c(0,0,2,4))
plot(NA, xlim=c(0, 1), ylim=c(0, 10),
     type="n", axes=FALSE, xlab="", ylab="")
rect(-0.8, 0:9, 0, 1:10, col=cols,
     border=NA, xpd=NA)
axis(side=4, at=0:10, labels=c(round(brks, 2)),
     las=1, line=-3.5, cex.axis=2)
par(las=0)
mtext("Human Modification index",
      side=4, line=1.5, cex=3)
#arrows(6, 0, 6, 10, xpd=NA,
#       length=0.2, lwd=3, code=3)
#mtext("Less modified", side=4, line=8.5,
#      adj=0.05, cex=2, font=3, xpd=NA)
#mtext("More modified", side=4, line=8.75,
#      adj=0.95, cex=2, font=3, xpd=NA)
par(fig=c(0.09, 0.45, 0.05, 0.33),
    new=TRUE, mar=c(0, 0, 0, 0))
plot(usa, col=ifelse(usa$STUSPS == "GA", "grey", "white"),
     mar=c(0, 0, 0, 0), axes=FALSE, box=TRUE,
     background="white", lwd=2)
lines(ext(buff1000), col="red", lwd=6)
arrows(-98.5, 41.4, -86.5, 34.39, lwd=6, col="red",
       length=0.15)
dev.off()

# end script!

## notes on origin of HMI
# extract Global human modification datasets of terrestrial ecosystems for 2022
# raster to site buffers
# https://zenodo.org/records/14502573
# https://www.nature.com/articles/sdata201667

# Produced in a Google Earth Engine script:
# /* Visualize the Global Human Modification v3 datasets, December 20, 2024
# Citations:
#     Theobald, D.M., Oakleaf, J., Moncrieff, G., and Kennedy, C.M. 2024. Global human modification 
# datasets of terrestrial ecosystems from 1990 to 2020. https://zenodo.org/14449495.
# Theobald, D.M., Oakleaf, J., Moncrieff, G., and Kennedy, C.M. 2024. Global human modification 
# datasets of terrestrial ecosystems for 2022. https://zenodo.org/14502573.
# Contacts: 
#     Dave Theobald, davidtheobald8@gmail.com
# Christina Kennedy, The Nature Conservancy, ckennedy@tnc.org
# 
# For 2022 and for each 5-year step from 1990 to 2020, 9 assets are provided (300 m). Values are floating-point numbers, 
# ranging from 0.0 - 1.0 (no modification, full modified).
# 
# The naming convention is as follows: 
#     - “c” signifies a data consistent for the change datasets for 1990-2020; 
# - “s” signifies 'static' data for 2022 that should not be analyzed against the change data series.
# - the IUCN taxonomy threat classes: 
#     - AA = all threats combined, 
# - AG = agricultural, 
# - BU =  residential, commercial and recreation areas, 
# - EX = energy production and mining; 
# - FR = biological resource use; 
# - HA = human accessibility; 
# - NS = natural systems modification; 
# - PO = pollution; and 
# - TI = transportation and service corridors. 
# */
#     
#     var paletteHM = ['4c6100', 'adda25', 'e2ff9b', 'ffff73', 'ffe629', 'ffd37f', 'ffaa00', 'e69808', 'e60000', 'a80000', '730000']
# 
# /////////////////// 2022 assets
# // the overall, cumulative HM, 300 m
# var HMoverall = ee.Image('projects/hm-30x30/assets/output/v20240801/HMv20240801_2022s_AA_300')
# Map.addLayer(HMoverall, {min:0, max:1.0, palette: paletteHM}, 'GHM-v3 2022 overall (300 m)')
# 
# // HM for the threat groups
# var lstThreats = ['AG', 'BU', 'EX', 'FR', 'HI', 'NS', 'PO', 'TI']
# for (var i=0; i<lstThreats.length; i++) {
#     var HM = ee.Image('projects/hm-30x30/assets/output/v20240801/HMv20240801_2022s_' + lstThreats[i] + '_300')
#     Map.addLayer(HM, {min:0, max:1.0, palette: paletteHM}, 'GHM 2022 ' + lstThreats[i] + ' (300 m)',false)
# }
# 
# // the overall, cumulative HM, 90 m
# var HMoverall = ee.Image('projects/hm-30x30/assets/output/v20240801/HMv20240801_2022s_AA_90')
# Map.addLayer(HMoverall, {min:0, max:1.0, palette: paletteHM}, 'GHM-v3 2022 overall (90 m)',false)
# 
# // HM for the threat groups
# var lstThreats = ['AG', 'BU', 'EX', 'FR', 'HI', 'NS', 'PO', 'TI']
# for (var i=0; i<lstThreats.length; i++) {
#     var HM = ee.Image('projects/hm-30x30/assets/output/v20240801/HMv20240801_2022s_' + lstThreats[i] + '_90')
#     Map.addLayer(HM, {min:0, max:1.0, palette: paletteHM}, 'GHM 2022 ' + lstThreats[i] + ' (90 m)',false)
# }
# 
# /////////////////// 1990-2020 assets
# // the overall, cumulative HM
# for (var y=1990; y<=2020; y+=5) {
#     var HMoverall = ee.Image('projects/hm-30x30/assets/output/v20240801/HMv20240801_' + y + 'c_AA_300')
#     Map.addLayer(HMoverall, {min:0, max:1.0, palette: paletteHM}, 'GHM-v3 ' + y + ' overall (300 m)',false)
# }
# for (var y=1990; y<=2020; y+=5) {
#     // HM for the threat groups
#     for (var i=0; i<lstThreats.length; i++) {
#         var HM = ee.Image('projects/hm-30x30/assets/output/v20240801/HMv20240801_' + y + 'c_' + lstThreats[i] + '_300')
#         Map.addLayer(HM, {min:0, max:1.0, palette: paletteHM}, 'GHM ' + y + ' ' + lstThreats[i] + ' (300 m)',false)
#     }
# }