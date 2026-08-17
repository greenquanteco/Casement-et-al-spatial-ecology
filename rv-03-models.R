# packages

library(MuMIn)
library(spdep)
library(lmtest)

# get response variables

dx <- read.csv("dat-mamm.csv")

# get explanatory variables
dy <- read.csv("dat-expl-100.csv")
dz <- read.csv("dat-human-modification-2022.csv")

# make sure in same order
all(dx$site == dy$site)
all(dy$site == dz$site)
rownames(dz) <- dz$site
dz <- dz[dx$site,]
all(dx$site == dz$site)
all(dy$site == dz$site)

nrow(dx)
nrow(dy)
nrow(dz)

# general explanatory matrix
names(dy)[which(names(dy) == "pct_imp_perim")] <- "perim"
ex.cols <- c("site", "areaha", "shape", "age", "perim", "island",
             "imp", "forest", "dev", "tree", "open",
             "popden", "povrate")
de <- dy[,ex.cols]
de$humanmod <- dz$human

# de is the matrix of explanatory variables
# to be used with all response variables

# check correlations between the explanatory variables
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(de[,-1], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0)

cormat <- round(cor(de[,-1]),3)

apply(de[,-1], 2, var)

apply(de[,-1], 2, function(x){length(unique(x))})

# make a data frame for each response
dx.rich <- data.frame(y=dx$rich, de[,-1])
dx.pele <- data.frame(y=dx$pele, de[,-1])
dx.sihi <- data.frame(y=dx$sihi, de[,-1])
dx.tast <- data.frame(y=dx$tast, de[,-1])

rownames(dx.rich) <- de$site
rownames(dx.pele) <- de$site
rownames(dx.sihi) <- de$site
rownames(dx.tast) <- de$site

dx.sihi$y <- ifelse(dx.sihi$y > 0, 1, 0)
dx.tast$y <- ifelse(dx.tast$y > 0, 1, 0)

# set up lists to hold model objects model dx.rich
list.rich <- vector("list", length=ncol(dx.rich))
names(list.rich) <- c("null", names(dx.rich)[-1])
list.pele <- list.rich
list.sihi <- list.rich
list.tast <- list.rich

# RICHNESS models
list.rich[[01]] <- glm(y~1,        data=dx.rich, family=poisson)
list.rich[[02]] <- glm(y~areaha,   data=dx.rich, family=poisson)
list.rich[[03]] <- glm(y~shape,    data=dx.rich, family=poisson)
list.rich[[04]] <- glm(y~age,      data=dx.rich, family=poisson)
list.rich[[05]] <- glm(y~perim,    data=dx.rich, family=poisson)
list.rich[[06]] <- glm(y~island,   data=dx.rich, family=poisson)
list.rich[[07]] <- glm(y~imp,      data=dx.rich, family=poisson)
list.rich[[08]] <- glm(y~forest,   data=dx.rich, family=poisson)
list.rich[[09]] <- glm(y~dev,      data=dx.rich, family=poisson)
list.rich[[10]] <- glm(y~tree,     data=dx.rich, family=poisson)
list.rich[[11]] <- glm(y~open,     data=dx.rich, family=poisson)
list.rich[[12]] <- glm(y~popden,   data=dx.rich, family=poisson)
list.rich[[13]] <- glm(y~povrate,  data=dx.rich, family=poisson)
list.rich[[14]] <- glm(y~humanmod, data=dx.rich, family=poisson)

aic.rich <- data.frame(mod=1:length(list.rich))
aic.rich$pred <- names(list.rich)
aic.rich$aicc <- sapply(list.rich, MuMIn::AICc)

aic.rich$delta <- aic.rich$aicc - min(aic.rich$aicc)
aic.rich$wt <- exp(-0.5*aic.rich$delta)
aic.rich$wt <- aic.rich$wt/sum(aic.rich$wt)

# order by descending AIC weight
aic.rich <- aic.rich[order(-aic.rich$wt),]

# calculate evidence ratio
aic.rich$ER <- max(aic.rich$wt) / aic.rich$wt

aic.rich

# get residuals
rich.res <- sapply(list.rich, residuals, type="pearson")

# get spatial weights
w <- readRDS("spatial/weights_15.7km.rds")
class(w)

# run permutation-based Moran's I across models
set.seed(123)
moran.rich <- t(apply(rich.res, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))

moran.rich <- as.data.frame(moran.rich)
moran.rich
#                   I      p
#null     -0.10154530 0.6371
#areaha   -0.05622163 0.4800
#age      -0.11827005 0.7016
#perim    -0.09831333 0.6348
#island   -0.09516099 0.6210
#imp      -0.07732384 0.5492
#forest   -0.08584304 0.5789
#dev      -0.06213514 0.4964
#tree     -0.09940847 0.6435
#open     -0.10934176 0.6746
#popden   -0.06074462 0.4925
#povrate  -0.10220549 0.6492
#humanmod -0.08576378 0.5843

# check Poisson dispersion and the 
# Pearson goodness of fit p-value
disp.rich <- t(sapply(list.rich, function(m) {
    
    X2 <- sum(residuals(m, type = "pearson")^2)
    rdf <- df.residual(m)
    
    c(X2 = X2,
      df = rdf,
      dispersion = X2 / rdf,
      p = pchisq(X2, df = rdf, lower.tail = FALSE))
}))

disp.rich <- round(disp.rich, 3)
disp.rich
#             X2 df dispersion     p
#null     20.812 22      0.946 0.532
#areaha   19.584 21      0.933 0.548
#shape    20.136 21      0.959 0.513
#age      20.180 21      0.961 0.510
#perim    18.334 21      0.873 0.628
#island   19.495 21      0.928 0.553
#imp      20.399 21      0.971 0.496
#forest   18.889 21      0.899 0.592
#dev      19.426 21      0.925 0.558
#tree     20.931 21      0.997 0.463
#open     20.826 21      0.992 0.470
#popden   19.411 21      0.924 0.559
#povrate  20.809 21      0.991 0.471
#humanmod 20.202 21      0.962 0.508

########################################################################

# PELE population density
head(dx.pele)
range(dx.pele$y)
sort(dx.pele$y)
mean(dx.pele$y)
sd(dx.pele$y)

list.pele[[01]] <- glm(y~1,        data=dx.pele)
list.pele[[02]] <- glm(y~areaha,   data=dx.pele)
list.pele[[03]] <- glm(y~shape,    data=dx.pele)
list.pele[[04]] <- glm(y~age,      data=dx.pele)
list.pele[[05]] <- glm(y~perim,    data=dx.pele)
list.pele[[06]] <- glm(y~island,   data=dx.pele)
list.pele[[07]] <- glm(y~imp,      data=dx.pele)
list.pele[[08]] <- glm(y~forest,   data=dx.pele)
list.pele[[09]] <- glm(y~dev,      data=dx.pele)
list.pele[[10]] <- glm(y~tree,     data=dx.pele)
list.pele[[11]] <- glm(y~open,     data=dx.pele)
list.pele[[12]] <- glm(y~popden,   data=dx.pele)
list.pele[[13]] <- glm(y~povrate,  data=dx.pele)
list.pele[[14]] <- glm(y~humanmod, data=dx.pele)

aic.pele <- data.frame(mod=1:length(list.pele))
aic.pele$pred <- names(list.pele)
aic.pele$aicc <- sapply(list.pele, MuMIn::AICc)

aic.pele$delta <- aic.pele$aicc - min(aic.pele$aicc)
aic.pele$wt <- exp(-0.5*aic.pele$delta)
aic.pele$wt <- aic.pele$wt/sum(aic.pele$wt)

# order by descending AIC weight
aic.pele <- aic.pele[order(-aic.pele$wt),]

# calculate evidence ratio
aic.pele$ER <- max(aic.pele$wt) / aic.pele$wt

aic.pele

# get residuals
pele.res <- sapply(list.pele, residuals, type="pearson")

# run permutation-based Moran's I across models
set.seed(123)
moran.pele <- t(apply(pele.res, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))

moran.pele <- as.data.frame(moran.pele)
moran.pele

par(mfrow = c(2, 2))
plot(list.pele[["dev"]])
plot(list.pele[["perim"]])

bptest(list.pele[["dev"]])
bptest(list.pele[["perim"]])

###############################################################

## try log(y+1) transform on pele
dx.pele$z <- log(dx.pele$y + 1)

list.pele2 <- list.pele

list.pele2[[01]] <- glm(z~1,        data=dx.pele)
list.pele2[[02]] <- glm(z~areaha,   data=dx.pele)
list.pele2[[03]] <- glm(z~shape,    data=dx.pele)
list.pele2[[04]] <- glm(z~age,      data=dx.pele)
list.pele2[[05]] <- glm(z~perim,    data=dx.pele)
list.pele2[[06]] <- glm(z~island,   data=dx.pele)
list.pele2[[07]] <- glm(z~imp,      data=dx.pele)
list.pele2[[08]] <- glm(z~forest,   data=dx.pele)
list.pele2[[09]] <- glm(z~dev,      data=dx.pele)
list.pele2[[10]] <- glm(z~tree,     data=dx.pele)
list.pele2[[11]] <- glm(z~open,     data=dx.pele)
list.pele2[[12]] <- glm(z~popden,   data=dx.pele)
list.pele2[[13]] <- glm(z~povrate,  data=dx.pele)
list.pele2[[14]] <- glm(z~humanmod, data=dx.pele)

aic.pele2 <- data.frame(mod=1:length(list.pele2))
aic.pele2$pred <- names(list.pele2)
aic.pele2$aicc <- sapply(list.pele2, MuMIn::AICc)

aic.pele2$delta <- aic.pele2$aicc - min(aic.pele2$aicc)
aic.pele2$wt <- exp(-0.5*aic.pele2$delta)
aic.pele2$wt <- aic.pele2$wt/sum(aic.pele2$wt)

# order by descending AIC weight
aic.pele2 <- aic.pele2[order(-aic.pele2$wt),]

# calculate evidence ratio
aic.pele2$ER <- max(aic.pele2$wt) / aic.pele2$wt

aic.pele2

# get residuals
pele.res2 <- sapply(list.pele2, residuals, type="pearson")

# run permutation-based Moran's I across models
set.seed(123)
moran.pele2 <- t(apply(pele.res2, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))

moran.pele2 <- as.data.frame(moran.pele2)
moran.pele2

par(mfrow = c(2, 2))
plot(list.pele2[["dev"]])
plot(list.pele2[["perim"]])

bptest(list.pele2[["dev"]])
bptest(list.pele2[["perim"]])

# conclusion: use log(y+1) version


# SIHI presence/absence models
list.sihi[[01]] <- glm(y~1,        data=dx.sihi, family=binomial)
list.sihi[[02]] <- glm(y~areaha,   data=dx.sihi, family=binomial)
list.sihi[[03]] <- glm(y~shape,    data=dx.sihi, family=binomial)
list.sihi[[04]] <- glm(y~age,      data=dx.sihi, family=binomial)
list.sihi[[05]] <- glm(y~perim,    data=dx.sihi, family=binomial)
list.sihi[[06]] <- glm(y~island,   data=dx.sihi, family=binomial)
list.sihi[[07]] <- glm(y~imp,      data=dx.sihi, family=binomial)
list.sihi[[08]] <- glm(y~forest,   data=dx.sihi, family=binomial)
list.sihi[[09]] <- glm(y~dev,      data=dx.sihi, family=binomial)
list.sihi[[10]] <- glm(y~tree,     data=dx.sihi, family=binomial)
list.sihi[[11]] <- glm(y~open,     data=dx.sihi, family=binomial)
list.sihi[[12]] <- glm(y~popden,   data=dx.sihi, family=binomial)
list.sihi[[13]] <- glm(y~povrate,  data=dx.sihi, family=binomial)
list.sihi[[14]] <- glm(y~humanmod, data=dx.sihi, family=binomial)

aic.sihi <- data.frame(mod=1:length(list.sihi))
aic.sihi$pred <- names(list.sihi)
aic.sihi$aicc <- sapply(list.sihi, MuMIn::AICc)

aic.sihi$delta <- aic.sihi$aicc - min(aic.sihi$aicc)
aic.sihi$wt <- exp(-0.5*aic.sihi$delta)
aic.sihi$wt <- aic.sihi$wt/sum(aic.sihi$wt)

# order by descending AIC weight
aic.sihi <- aic.sihi[order(-aic.sihi$wt),]

# calculate evidence ratio
aic.sihi$ER <- max(aic.sihi$wt) / aic.sihi$wt

aic.sihi

# get residuals
rich.res <- sapply(list.rich, residuals, type="pearson")

# get spatial weights
w <- readRDS("spatial/weights_15.7km.rds")
class(w)

# run permutation-based Moran's I across models
set.seed(123)
moran.rich <- t(apply(rich.res, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))

moran.rich <- as.data.frame(moran.rich)
moran.rich
#                   I      p
#null     -0.10154530 0.6371
#areaha   -0.05622163 0.4800
#age      -0.11827005 0.7016
#perim    -0.09831333 0.6348
#island   -0.09516099 0.6210
#imp      -0.07732384 0.5492
#forest   -0.08584304 0.5789
#dev      -0.06213514 0.4964
#tree     -0.09940847 0.6435
#open     -0.10934176 0.6746
#popden   -0.06074462 0.4925
#povrate  -0.10220549 0.6492
#humanmod -0.08576378 0.5843

# check Poisson dispersion and the 
# Pearson goodness of fit p-value
disp.rich <- t(sapply(list.rich, function(m) {
    
    X2 <- sum(residuals(m, type = "pearson")^2)
    rdf <- df.residual(m)
    
    c(X2 = X2,
      df = rdf,
      dispersion = X2 / rdf,
      p = pchisq(X2, df = rdf, lower.tail = FALSE))
}))

disp.rich <- round(disp.rich, 3)
disp.rich
