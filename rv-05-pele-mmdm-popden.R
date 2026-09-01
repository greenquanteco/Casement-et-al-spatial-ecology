# need to recalculate population density

# MMDM for each species
# individuals added back in, using script "06 mark recap set up 2025-01-20.r",
# were not included in MMDM calculations.
mm <- read.csv("dat-spp-mmdm.csv")

# transect coordinates and lengths
tr <- read.csv("dat-trans-coords.csv")

# biodata 
## use capture history because we need MNKA, 
## NOT number of captures
bio <- read.csv("dat-cap-history.csv")

# don't need coordinates
keep.cols <- c("sitena", "utran", "length")
tr <- tr[,keep.cols]
tr$length <- round(tr$length, 2)

# calculate species-specific areas using MMDM as
# radii. Equation is:
# A = 2rL + pi*(mmdm/2)^2
pm <- mm$mmdm[which(mm$spp == "pele")]

tr$aream2 <- (tr$length * pm) + (pi * (pm/2)^2)
tr$areaha <- tr$aream2 / 1e4
tr <- tr[order(tr$utran),]

# aggregate tr to site
## first, total area sampled for each species.
## combine this with total number captured at a site.
## later, can compare this to "mean of transects" values.
site.sum <- aggregate(areaha~sitena, data=tr, sum)

# merge with biodata 
## site level MNKA
## counts individuals that split their time among 2 transects as "half"
## an individual at each
### how many animals each capture history counts as
### most are 1, but animals counted as half on each transect if seen on 2
### so they are not double counted

bio <- bio[which(bio$spp != "divi"),]
bio$indmax <- apply(bio[,6:ncol(bio)], 1, max)
agg1 <- aggregate(indmax~spp+site, data=bio, sum)

# calculate pop density using whole site data
den1 <- data.frame(site=site.sum$sitena)

# matrix of sites x species
spps <- sort(unique(bio$spp))
nspp <- length (spps)

ma <- matrix(0, nrow=nrow(site.sum), ncol=nspp)
colnames(ma) <- spps
rownames(ma) <- site.sum$sitena
ma <- as.data.frame(ma)

# make some vectors of names and their lengths
sites <- rownames(ma)
nsites <- length(sites)

# must be TRUE
all(colnames(ma) %in% mm$spp)

# fill in ma in a loop
for(i in 1:nsites){
    for(j in 1:nspp){
        flag1 <- which(agg1$site == sites[i])
        flag2 <- which(agg1$spp == spps[j])
        flag3 <- intersect(flag1, flag2)
        if(length(flag3) == 1){ma[sites[i],spps[j]] <- agg1$indmax[flag3]}
        if(length(flag3) > 1){break}    
    }#j
}#i

# must be TRUE
all(site.sum$sitena == sites)
all(den1$site == rownames(ma))

# calculate density
rownames(ma)==site.sum$sitena

popden <- ma
for(i in 1:nspp){
    popden[,i] <- round(ma[,i] / site.sum$areaha, 3)
}#i

# combine with site for site-level pop density
den1 <- cbind(den1, popden)
rownames(den1) <- NULL

# write out population densities
# units are individuals / ha
# aka       MNKA / ha
den1$site[which(den1$site == "field station")] <- "field"
write.csv(den1, "dat-popden.csv", row.names=FALSE)

# end script!