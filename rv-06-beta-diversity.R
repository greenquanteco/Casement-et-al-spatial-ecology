dx <- read.csv("dat-popden.csv")

# get explanatory variables
dz <- read.csv("dat-human-modification-2022.csv")
rownames(dz) <- dz$site
dz <- dz[dx$site,]
all(dx$site == dz$site)

rownames(dx) <- dx$site
dx$site <- NULL

# Convert abundance/density data to presence-absence
pa <- (dx > 0) * 1

# Pairwise beta diversity
library(betapart)

beta <- beta.pair(pa, index.family = "sorensen")

# Components
bsor <- beta$beta.sor   # total beta diversity
bsim <- beta$beta.sim   # turnover
bsne <- beta$beta.sne   # nestedness

dhuman <- dist(dz$human)

identical(rownames(dx), rownames(dz))

par(mfrow = c(1, 3))

plot(dhuman, bsor,
     xlab = "Difference in human modification",
     ylab = expression(beta[SOR]))

plot(dhuman, bsim,
     xlab = "Difference in human modification",
     ylab = expression(beta[SIM]))

plot(dhuman, bsne,
     xlab = "Difference in human modification",
     ylab = expression(beta[SNE]))

cor.test(dhuman, bsor)
cor.test(dhuman, bsim)
cor.test(dhuman, bsne)


library(vegan)
set.seed(1)
mantel(bsor, dhuman,
       method = "pearson",
       permutations = 9999)

mantel(bsim, dhuman,
       method = "pearson",
       permutations = 9999)

mantel(bsne, dhuman,
       method = "pearson",
       permutations = 9999)


library(ecodist)

site_dist_km <- readRDS("spatial/site_distance_km.rds")
diag(site_dist_km) <- 0
rownames(site_dist_km) <- rownames(dx)
colnames(site_dist_km) <- rownames(dx)
dhuman <- dist(setNames(dz$human, rownames(dz)))

dgeo <- as.dist(site_dist_km)
str(dgeo)
length(dgeo)
str(bsor)


str(bsor)
str(dhuman)
str(dgeo)

mrm.dat <- data.frame(
    bsor   = as.vector(bsor),
    bsim   = as.vector(bsim),
    bsne   = as.vector(bsne),
    dhuman = as.vector(dhuman),
    dgeo   = as.vector(dgeo)
)

set.seed(123)

MRM(bsor ~ dhuman + dgeo,
    data = mrm.dat,
    nperm = 9999)

MRM(bsim ~ dhuman + dgeo,
    data = mrm.dat,
    nperm = 9999)

MRM(bsne ~ dhuman + dgeo,
    data = mrm.dat,
    nperm = 9999)


cor(dhuman, dgeo)

vegan::mantel(dhuman, dgeo,
       method = "pearson",
       permutations = 9999)









set.seed(42)
MRM(bsor ~ dhuman + dgeo,
    nperm = 9999)
str(bsor)
str(dhuman)
str(site_dist_km)

MRM(bsim ~ dhuman + site_dist_km,
    nperm = 9999)

MRM(bsne ~ dhuman + site_dist_km,
    nperm = 9999)
