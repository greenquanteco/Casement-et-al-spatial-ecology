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

library(terra)

# get sites
sites <- vect("spatial/sites.gpkg")

buf <- buffer(sites, 100)

plot(buf, col="blue")
plot(sites, add=TRUE, col="white", border=NA)

human <- rast("K:/mammal_ga/human_modification_2022_90m.tif")

crs(human)
crs(sites)

sitesp <- project(sites, crs(human))

buff1000 <- buffer(sitesp, 1000)
buff100 <- buffer(sitesp, 100)
donut100 <- erase(buff100, sitesp)

humanc <- crop(human, ext(buff1000))

plot(humanc)
plot(buff100, add=TRUE, col="red")

plot(humanc, xlim=c(-84.703, -84.645),
     ylim=c(34.014, 34.068))
plot(donut100, add=TRUE,
     border="red", col="transparent", lwd=4)

# extract HM to donuts
hm.mean <- extract(humanc, donut100,
                   fun=mean, weights=TRUE,
                   na.rm=TRUE)
hm.mean

hm.mean$site <- sitesp$name

# Zion returned NaN

plot(sitesp, col=ifelse(sitesp$name=="zion", "red", "blue"))

plot(humanc, xlim=c(-84.406, -84.362),
     ylim=c(33.657, 33.682))
plot(donut100, col="transparent", border="red", lwd=2, add=TRUE)
#plot(sitesp[which(sitesp$name=="zion"),], col="green", add=TRUE)



z <- extract(
    humanc,
    donut100[22, ],
    weights = TRUE
)

z
nrow(z)
sum(!is.na(z[, 2]))
sum(z$weight, na.rm = TRUE)
# Is the polygon geometrically valid?
is.valid(donut100[22, ])

# Does it have nonzero area?
expanse(donut100[22, ])

# Does its extent overlap the raster?
ext(donut100[22, ])
ext(humanc)

# Extract without weights
extract(hm, donut.hm[i, ])

extract(humanc, donut100[22,], exact=TRUE)

# redo them all in a loop with exact=TRUE,
# avoids issues with polygons being read weird
hm <- numeric(nrow(donut100))
cels <- numeric(nrow(donut100))

for (i in 1:nrow(donut100)) {
    
    z <- extract(
        humanc,
        donut100[i, ],
        exact = TRUE
    )
    
    hm[i] <- weighted.mean(
        z$All_threats_combined,
        z$fraction,
        na.rm = TRUE
    )
    cels[i] <- nrow(z)
}

hm
cels
# make an output
out <- data.frame(site=sitesp$name, human=hm, cells=cels)

write.csv(out, "dat-human-modification-2022.csv", row.names=FALSE)                  
