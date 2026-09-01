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
# A = 2rL + pi*r^2
# basically assumes a rectangle, with a half-circle at each end
ma <- matrix(NA, nrow=nrow(tr), ncol=nrow(mm))
colnames(ma) <- paste("ar", mm$spp, sep="_")
for(i in 1:ncol(ma)){
    ma[,i] <- (2*mm$mmdm[i]*tr$length) + (pi * mm$mmdm[i] * mm$mmdm[i])
}#i
ma <- round(ma, 1)
tr <- cbind(tr, as.data.frame(ma))
tr <- tr[order(tr$utran),]

# aggregate tr to site
## first, total area sampled for each species.
## combine this with total number captured at a site.
## later, can compare this to "mean of transects" values.
site.sum <- aggregate(ar_mipi~sitena, data=tr, sum)
site.sum$ar_mumu <- aggregate(ar_mumu~sitena, data=tr, sum)$ar_mumu
site.sum$ar_pego <- aggregate(ar_pego~sitena, data=tr, sum)$ar_pego
site.sum$ar_pele <- aggregate(ar_pele~sitena, data=tr, sum)$ar_pele
site.sum$ar_rano <- aggregate(ar_rano~sitena, data=tr, sum)$ar_rano
site.sum$ar_rehu <- aggregate(ar_rehu~sitena, data=tr, sum)$ar_rehu
site.sum$ar_sihi <- aggregate(ar_sihi~sitena, data=tr, sum)$ar_sihi
site.sum$ar_tast <- aggregate(ar_tast~sitena, data=tr, sum)$ar_tast
site.sum$ar_blca <- aggregate(ar_blca~sitena, data=tr, sum)$ar_blca
site.sum$ar_glvo <- aggregate(ar_glvo~sitena, data=tr, sum)$ar_glvo
site.sum$ar_ocnu <- aggregate(ar_ocnu~sitena, data=tr, sum)$ar_ocnu
site.sum$ar_scca <- aggregate(ar_scca~sitena, data=tr, sum)$ar_scca

# convert square meters to ha
site.sum$ar_mipi <- site.sum$ar_mipi / 1e4
site.sum$ar_mumu <- site.sum$ar_mumu / 1e4
site.sum$ar_pego <- site.sum$ar_pego / 1e4
site.sum$ar_pele <- site.sum$ar_pele / 1e4
site.sum$ar_rano <- site.sum$ar_rano / 1e4
site.sum$ar_rehu <- site.sum$ar_rehu / 1e4
site.sum$ar_sihi <- site.sum$ar_sihi / 1e4
site.sum$ar_tast <- site.sum$ar_tast / 1e4
site.sum$ar_blca <- site.sum$ar_blca / 1e4
site.sum$ar_glvo <- site.sum$ar_glvo / 1e4
site.sum$ar_ocnu <- site.sum$ar_ocnu / 1e4
site.sum$ar_scca <- site.sum$ar_scca / 1e4

# check
## should be 0s
apply(site.sum[,-1], 2, function(x){length(which(is.na(x)))})

# merge with biodata 
## site level MNKA
## counts individuals that split their time among 2 transects as "half"
## an individual at each
### how many animals each capture history counts as
### most are 1, but animals counted as half on each transect if seen on 2
### so they are not double counted
bio$indmax <- apply(bio[,6:ncol(bio)], 1, max)
agg1 <- aggregate(indmax~spp+site, data=bio, sum)

## transect level MNKA
agg2 <- aggregate(indmax~spp+tran, data=bio, sum)

# calculate pop density using whole site data
den1 <- data.frame(site=site.sum$sitena)
ma <- matrix(0, nrow=nrow(site.sum), ncol=ncol(site.sum)-1)
colnames(ma) <- sapply(strsplit(names(site.sum)[-1], "_"), "[", 2)
ma <- as.data.frame(ma)

# make some vectors of names and their lengths
spp <- colnames(ma)
nspp <- length(spp)
sites <- den1$site
nsites <- length(sites)

# must be TRUE
all(colnames(ma) %in% mm$spp)

rownames(ma) <- sites

# fill in ma in a loop
for(i in 1:nsites){
    for(j in 1:nspp){
        flag1 <- which(agg1$site == sites[i])
        flag2 <- which(agg1$spp == spp[j])
        flag3 <- intersect(flag1, flag2)
        if(length(flag3) == 1){ma[i,j] <- agg1$indmax[flag3]}
        if(length(flag3) > 1){break}    
    }#j
}#i

# must be TRUE
all(site.sum$sitena == sites)
all(rownames(ma) == site.sum$sitena)

# calculate density
for(i in 1:nsites){
    ma[i,] <- ma[i,] / site.sum[i,-1]
}#i

# combine with site for site-level pop density
den1 <- cbind(den1, ma)
rownames(den1) <- NULL

# calculate transect level density
den2 <- data.frame(utran=tr$utran)

utrans <- tr$utran
ntrans <- length(utrans)

ma <- matrix(0, nrow=ntrans, ncol=nspp)
colnames(ma) <- spp
rownames(ma) <- utrans

# can be FALSE because not all transects caught an animal;
# utrans contains all transects but agg2 contains only 
# those that had an animal.
all(utrans %in% agg2$tran)

# fill in ma in a loop
for(i in 1:ntrans){
    flag1 <- which(agg2$tran == utrans[i])
    if(length(flag1) == 0){next}
    for(j in 1:nspp){
        flag2 <- which(agg2$spp == spp[j])
        flag3 <- intersect(flag1, flag2)
        if(length(flag3) == 0){next} 
        
        ma[i,j] <- agg2$indmax[flag3]
    }#j
}#i

# get columns with areas from tr data frame
ar.cols <- grep("ar_", names(tr))
for(i in 1:ntrans){
    flag1 <- which(tr$utran == utrans[i])
    # pull species-specific areas and convert to ha
    i.areas <- as.numeric(tr[flag1, ar.cols]) / 1e4
    i.mnka <- ma[i,]
    ma[i,] <- i.mnka / i.areas
}

den2 <- cbind(den2, as.data.frame(ma))
rownames(den2) <- NULL

# write out population densities
# units are individuals / ha
# aka       MNKA / ha
out.name <- "dat_rodent_popden_site_2025-02-25.csv"
out.path <- paste(use.dir, out.name, sep="/")
write.csv(den1, out.path, row.names=FALSE)

out.name <- "dat_rodent_popden_transect_2025-02-25.csv"
out.path <- paste(use.dir, out.name, sep="/")
write.csv(den2, out.path, row.names=FALSE)

# end script!