# script to investigate potential autocorrelation in
# response variables before revision 1 reanalysis

library(terra)
library(sp)
library(spdep)

# site polygons
sites <- vect("spatial/sites.gpkg")
plot(sites)

# sites in analysis
dx <- read.csv("dat-expl-100.csv")

# site mammal data
dy <- read.csv("dat-popden.csv")


# do we have all the right site names?
asites <- dx$site
ssites <- sites$name
## must be TRUE
all(asites %in% ssites)

# just need sites in asites
sites <- sites[which(sites$name %in% asites),]

# must be all TRUE
sort(sites$name) == sort(asites)

# reorder to match
sites <- sites[match(dx$site, sites$name),]
## check
dx$site == sites$name

# check mammal data
## TRUE
all(dy$site == dx$site)
 
rownames(dy) <- dy$site
dy$site <- NULL

rownames(dx) <- dx$site

# check coordinates and get centroids
crs(sites)
is.lonlat(sites)

# project to utm
sites.utm <- project(sites, "EPSG:32616")
is.lonlat(sites.utm)

# get centroids
cent <- centroids(sites.utm)

# distance matrix
xy <- crds(cent)
d <- as.matrix(dist(xy)) # in meters
d.km <- d/1000           # in km

summary(d.km[upper.tri(d.km)])

diag(d.km) <- NA

# nearest neighbors
nearest <- apply(d.km, 1, min, na.rm = TRUE)

summary(nearest)
sort(nearest)
hist(d.km[upper.tri(d.km)],
     xlab = "Distance between site centroids (km)",
     main = "")

# find neighbors in 10 km distance class
nb10 <- dnearneigh(xy, d1 = 0, d2 = 10000)
card(nb10)

# inspect progressively larger thresholds
thresholds <- seq(10, 30, by = 2)

ncomp <- sapply(thresholds, function(x) {
    nb <- dnearneigh(xy, d1 = 0, d2 = x * 1000)
    n.comp.nb(nb)$nc
})

data.frame(
    distance_km = thresholds,
    components = ncomp
)

# refine based on cut off of 14-16 km
# in previous result
thresholds <- seq(14, 16, by = 0.1)

ncomp <- sapply(thresholds, function(x) {
    nb <- dnearneigh(xy, d1 = 0, d2 = x * 1000)
    n.comp.nb(nb)$nc
})

thresholds[which(ncomp == 1)[1]]

# build the neighborhood of sites using
# the threshold of 15.7 km
nb <- dnearneigh(xy, d1 = 0, d2 = 15700)
card(nb)
summary(card(nb))

# nb is essentially a list defining the neighborhood:
# site 1 is connected to sites x, y, z; site 2 is connected
# to sites w, x; and so on.

# now to convert that neighborhood to spatial weights
# for calculating Moran's I. Simplest choice is 
# row-standardized weights
lw <- nb2listw(nb, style="W")

# now calculate and test Moran's I vs. null expectation
# under spatial randomness. Use permutation test because
# n is only 23
dy$rich <- apply(dy, 1, function(x){sum(x > 0)})

set.seed(123)
moran.mc(dy$rich, lw, nsim=9999, alternative="two.sided")
moran.mc(log(dy$pele+1), lw, nsim=9999, alternative="two.sided")
moran.mc(as.numeric(dy$sihi > 0), lw, nsim=9999, alternative="two.sided")
moran.mc(as.numeric(dy$tast > 0), lw, nsim=9999, alternative="two.sided")

# profile Moran's I across different non-overlapping distance
# classes

# check how many unique pairs are in each distance band
breaks <- seq(0, ceiling(max(d.km, na.rm = TRUE) / 10) * 10, by = 10)
d.pairs <- d.km[upper.tri(d.km)]
table(cut(d.pairs,
          breaks = breaks,
          include.lowest = TRUE,
          right = FALSE))

# compute the Moran's I profile at different bands
lower <- seq(0, 70, by = 10)
upper <- seq(10, 80, by = 10)

rich.corr <- data.frame(
    lower = lower,
    upper = upper,
    midpoint = (lower + upper) / 2,
    I = NA,
    p = NA
)

set.seed(42)

for (i in seq_along(lower)) {
    
    nb.i <- dnearneigh(
        xy,
        d1 = lower[i] * 1000,
        d2 = upper[i] * 1000
    )
    
    lw.i <- nb2listw(
        nb.i,
        style = "W",
        zero.policy = TRUE
    )
    
    test.i <- moran.mc(
        dy$rich,
        lw.i,
        nsim = 9999,
        alternative = "two.sided",
        zero.policy = TRUE
    )
    
    rich.corr$I[i] <- test.i$statistic
    rich.corr$p[i] <- test.i$p.value
}

rich.corr

# plot the autocorrelation correlogram:
plot(
    rich.corr$midpoint,
    rich.corr$I,
    type = "b",
    pch = 16,
    xlab = "Distance (km)",
    ylab = "Moran's I"
)

abline(h = -1 / (nrow(dy) - 1), lty = 2)
text(35, -0.058, "Null expectation with n = 23", adj=0)
text(35, -0.088, expression(italic(I)[null]==-0.0455), adj=0)

# automate across the other variables

## global Moran's I
vars <- c("rich", "simp", "pele", "sihi", "tast")

set.seed(100)
global <- data.frame(
    variable=vars,
    I=NA, p=NA
)

for (i in seq_along(vars)) {
    
    test.i <- moran.mc(
        dy[[vars[i]]],
        lw,
        nsim = 9999,
        alternative = "two.sided"
    )
    
    global$I[i] <- test.i$statistic
    global$p[i] <- test.i$p.value
    cat("variable", vars[i], "complete!\n")
    flush.console()
}
global

# profile each response
corr <- data.frame()

set.seed(123)

for (v in vars) {
    
    for (i in seq_along(lower)) {
        
        nb.i <- dnearneigh(
            xy,
            d1 = lower[i] * 1000,
            d2 = upper[i] * 1000
        )
        
        lw.i <- nb2listw(
            nb.i,
            style = "W",
            zero.policy = TRUE
        )
        
        test.i <- moran.mc(
            dy[[v]],
            lw.i,
            nsim = 9999,
            alternative = "two.sided",
            zero.policy = TRUE
        )
        
        corr <- rbind(
            corr,
            data.frame(
                variable = v,
                lower = lower[i],
                upper = upper[i],
                midpoint = (lower[i] + upper[i]) / 2,
                I = as.numeric(test.i$statistic),
                p = test.i$p.value
            )
        )
    }#i for distance
    cat("variable", vars[v], "complete!\n")
    flush.console()
}#v for vars

corr
corr[corr$p < 0.05,]

# Holm correction for each response
corr$p.holm <- ave(
    corr$p,
    corr$variable,
    FUN = function(x) p.adjust(x, method = "holm")
)

global
#  variable            I      p
#1     rich -0.101545295 0.7344
#2     simp -0.129580506 0.5310
#3     pele  0.187350894 0.0914
#4     sihi  0.009969548 0.4444
#5     tast  0.044573766 0.3406

corr[corr$p < 0.1, ]
#   variable lower upper midpoint          I      p p.holm
#4      rich    30    40       35  0.2661259 0.0626 0.5008
#12     simp    30    40       35  0.2747481 0.0640 0.5120
#17     pele     0    10        5  0.3886125 0.0348 0.2436
#21     pele    40    50       45  0.3548870 0.0186 0.1488
#22     pele    50    60       55  0.2429485 0.0824 0.4944
#28     sihi    30    40       35  0.2878325 0.0628 0.4396
#32     sihi    70    80       75 -0.2224064 0.0210 0.1680
#33     tast     0    10        5  0.3852587 0.0378 0.3024

# Takeaways:
# 1. There was no evidence of global spatial autocorrelation
#    in any of the 5 response variables. 
# 2. PELE exhibited some potential positive AC at the closest
#    and some farther distance bands.
# 3. SIHI, Rich, and Simp shows some potential positive
#    autocorrelation at the 30-40 band
# 4. TAST had some potential AC at closest band (0-10)

# now try to model the responses and check for autocorrelation
# in the residuals

# Matrix of pairwise centroid-to-centroid distances among sites, in km
# Class: matrix; requires base R only.
saveRDS(d.km, "spatial/site_distance_km.rds")

# Row-standardized spatial weights based on the 15.7-km neighbor network
# Class: listw; requires package spdep.
saveRDS(lw, "spatial/weights_15.7km.rds")

# Global Moran's I permutation-test results for the five focal response variables
# Class: data.frame; requires base R only.
saveRDS(global, "results/moran_global.rds")

# Moran's I correlogram results for the five focal responses in 10-km distance classes
# Class: data.frame; requires base R only.
saveRDS(corr, "results/moran_correlograms.rds")



