# script to convert the old master shapefile
# that lives on K drive to a GeoPackage to 
# facilitate version control

library(terra)

# read master site boundaries from campus storage
sites <- vect("K:/mammal_ga/gis/master_sites_2025-10-21.shp")


# Take a quick look before copying
sites
crs(sites)
plot(sites)

# fix field station name
data.frame(sites)
sites$name[1] <- "field"

# write a GeoPackage
writeVector(
    sites,
    "D:/_research/mammal_ga/Casement-et-al-spatial-ecology/spatial/sites.gpkg",
    filetype = "GPKG",
    overwrite = TRUE
)


sites2 <- vect("D:/_research/mammal_ga/Casement-et-al-spatial-ecology/spatial/sites.gpkg")

# check it
plot(sites2, col="blue")
data.frame(sites2)
nrow(sites)
nrow(sites2)
crs(sites)
crs(sites2)
crs(sites)==crs(sites2)

round(ext(sites2),2)
# end script!