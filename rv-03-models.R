# packages

library(MuMIn)
library(spdep)
library(lmtest)

# get response variables

dx <- read.csv("dat-popden.csv")

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
dx$rich <- apply(dx[,-1], 1, function(x){sum(x > 0)})

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
#    mod     pred     aicc    delta         wt       ER
# 1    1     null 85.88121 0.000000 0.15826404 1.000000
# 5    5    perim 87.05608 1.174868 0.08795546 1.799366
# 2    2   areaha 87.14798 1.266770 0.08400530 1.883977
# 3    3    shape 87.28309 1.401881 0.07851770 2.015648
# 9    9      dev 87.39578 1.514572 0.07421594 2.132481
# 12  12   popden 87.40780 1.526588 0.07377136 2.145332
# 8    8   forest 87.53137 1.650162 0.06935121 2.282066
# 6    6   island 87.82717 1.945957 0.05981677 2.645814
# 4    4      age 87.90561 2.024399 0.05751612 2.751647
# 11  11     open 88.04983 2.168620 0.05351462 2.957398
# 7    7      imp 88.06457 2.183358 0.05312171 2.979272
# 14  14 humanmod 88.06601 2.184797 0.05308351 2.981416
# 13  13  povrate 88.23971 2.358502 0.04866762 3.251937
# 10  10     tree 88.25908 2.377868 0.04819865 3.283578


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
#                    I      p
# null     -0.10154530 0.6371
# areaha   -0.05622163 0.4800
# shape    -0.07920240 0.5502
# age      -0.11827005 0.7016
# perim    -0.09831333 0.6348
# island   -0.09516099 0.6210
# imp      -0.07732384 0.5492
# forest   -0.08584304 0.5789
# dev      -0.06213514 0.4964
# tree     -0.09940847 0.6435
# open     -0.10934176 0.6746
# popden   -0.06074462 0.4925
# povrate  -0.10220549 0.6492
# humanmod -0.08576378 0.5843

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
# use log(y + 1) transform
head(dx.pele)
range(dx.pele$y)
sort(dx.pele$y)
mean(dx.pele$y)
sd(dx.pele$y)

dx.pele$z <- log(dx.pele$y + 1)

list.pele[[01]] <- lm(z~1,        data=dx.pele)
list.pele[[02]] <- lm(z~areaha,   data=dx.pele)
list.pele[[03]] <- lm(z~shape,    data=dx.pele)
list.pele[[04]] <- lm(z~age,      data=dx.pele)
list.pele[[05]] <- lm(z~perim,    data=dx.pele)
list.pele[[06]] <- lm(z~island,   data=dx.pele)
list.pele[[07]] <- lm(z~imp,      data=dx.pele)
list.pele[[08]] <- lm(z~forest,   data=dx.pele)
list.pele[[09]] <- lm(z~dev,      data=dx.pele)
list.pele[[10]] <- lm(z~tree,     data=dx.pele)
list.pele[[11]] <- lm(z~open,     data=dx.pele)
list.pele[[12]] <- lm(z~popden,   data=dx.pele)
list.pele[[13]] <- lm(z~povrate,  data=dx.pele)
list.pele[[14]] <- lm(z~humanmod, data=dx.pele)

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

# round
aic.pele[,3:6] <- round(aic.pele[,3:6], 3)

aic.pele
#    mod     pred   aicc  delta    wt      ER
# 5    5    perim 64.282  0.000 0.623   1.000
# 9    9      dev 66.069  1.787 0.255   2.444
# 11  11     open 69.513  5.231 0.046  13.675
# 14  14 humanmod 71.409  7.127 0.018  35.290
# 8    8   forest 71.750  7.468 0.015  41.849
# 3    3    shape 72.680  8.398 0.009  66.623
# 7    7      imp 73.173  8.891 0.007  85.250
# 1    1     null 73.246  8.964 0.007  88.432
# 12  12   popden 73.361  9.079 0.007  93.661
# 4    4      age 73.765  9.483 0.005 114.608
# 10  10     tree 75.118 10.836 0.003 225.425
# 13  13  povrate 75.639 11.357 0.002 292.479
# 2    2   areaha 75.753 11.471 0.002 309.684
# 6    6   island 75.834 11.552 0.002 322.492


# get residuals
pele.res <- sapply(list.pele, residuals, type="pearson")

# run permutation-based Moran's I across models
set.seed(123)
moran.pele <- t(apply(pele.res, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))
moran.pele <- round(moran.pele, 3)
moran.pele <- as.data.frame(moran.pele)
moran.pele <- moran.pele[aic.pele$pred,]
moran.pele
#               I     p
# perim     0.059 0.191
# dev       0.052 0.198
# open     -0.073 0.546
# humanmod  0.069 0.166
# forest    0.234 0.028
# shape     0.081 0.143
# imp       0.100 0.128
# null      0.204 0.039
# popden    0.086 0.138
# age       0.108 0.120
# tree      0.181 0.053
# povrate   0.170 0.058
# areaha    0.206 0.038
# island    0.231 0.026

par(mfrow = c(2, 2))
plot(list.pele[["dev"]])
plot(list.pele[["perim"]])

summary(list.pele$dev)

summary(list.pele$perim)


bptest(list.pele[["perim"]])
# 
# studentized Breusch-Pagan test
# 
# data:  list.pele[["perim"]]
# BP = 2.9944, df = 1, p-value = 0.08355

bptest(list.pele[["dev"]])
# studentized Breusch-Pagan test
# 
# data:  list.pele[["dev"]]
# BP = 0.090189, df = 1, p-value = 0.7639

par(mfrow=c(1,2))
plot(dx.pele$perim, dx.pele$y)
plot(dx.pele$dev, dx.pele$y)

##############################################################################

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

# check aic table
aic.sihi[,3:6] <- round(aic.sihi[,3:6], 3)
aic.sihi
#    mod     pred   aicc  delta    wt      ER
# 4    4      age 22.114  0.000 0.815   1.000
# 11  11     open 25.736  3.622 0.133   6.115
# 13  13  povrate 31.078  8.964 0.009  88.391
# 1    1     null 31.911  9.796 0.006 134.053
# 5    5    perim 32.109  9.995 0.006 148.043
# 10  10     tree 32.170 10.056 0.005 152.634
# 14  14 humanmod 32.334 10.219 0.005 165.625
# 12  12   popden 32.502 10.388 0.005 180.185
# 9    9      dev 32.906 10.792 0.004 220.481
# 7    7      imp 32.985 10.871 0.004 229.353
# 6    6   island 33.567 11.453 0.003 306.858
# 3    3    shape 33.799 11.685 0.002 344.583
# 2    2   areaha 34.036 11.922 0.002 388.036
# 8    8   forest 34.318 12.203 0.002 446.601


# get residuals
sihi.res <- sapply(list.sihi, residuals, type="pearson")

# run permutation-based Moran's I across models
set.seed(123)
moran.sihi <- t(apply(sihi.res, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))

moran.sihi <- round(moran.sihi, 3)
moran.sihi <- as.data.frame(moran.sihi)
moran.sihi <- moran.sihi[aic.sihi$pred,]
moran.sihi
#               I     p
# age      -0.092 0.619
# open     -0.052 0.469
# povrate   0.122 0.100
# null      0.293 0.014
# perim     0.237 0.028
# tree      0.241 0.026
# humanmod  0.121 0.095
# popden    0.157 0.069
# dev       0.151 0.074
# imp       0.148 0.074
# island    0.427 0.003
# shape     0.252 0.023
# areaha    0.263 0.024
# forest    0.290 0.017

# check coefficients
summary(list.sihi$age)$coefficients

summary(list.sihi$open)$coefficients

##############################################################################

# TAST presence/absence models

list.tast[[01]] <- glm(y~1,        data=dx.tast, family=binomial)
list.tast[[02]] <- glm(y~areaha,   data=dx.tast, family=binomial)
list.tast[[03]] <- glm(y~shape,    data=dx.tast, family=binomial)
list.tast[[04]] <- glm(y~age,      data=dx.tast, family=binomial)
list.tast[[05]] <- glm(y~perim,    data=dx.tast, family=binomial)
list.tast[[06]] <- glm(y~island,   data=dx.tast, family=binomial)
list.tast[[07]] <- glm(y~imp,      data=dx.tast, family=binomial)
list.tast[[08]] <- glm(y~forest,   data=dx.tast, family=binomial)
list.tast[[09]] <- glm(y~dev,      data=dx.tast, family=binomial)
list.tast[[10]] <- glm(y~tree,     data=dx.tast, family=binomial)
list.tast[[11]] <- glm(y~open,     data=dx.tast, family=binomial)
list.tast[[12]] <- glm(y~popden,   data=dx.tast, family=binomial)
list.tast[[13]] <- glm(y~povrate,  data=dx.tast, family=binomial)
list.tast[[14]] <- glm(y~humanmod, data=dx.tast, family=binomial)

aic.tast <- data.frame(mod=1:length(list.tast))
aic.tast$pred <- names(list.tast)
aic.tast$aicc <- sapply(list.tast, MuMIn::AICc)

aic.tast$delta <- aic.tast$aicc - min(aic.tast$aicc)
aic.tast$wt <- exp(-0.5*aic.tast$delta)
aic.tast$wt <- aic.tast$wt/sum(aic.tast$wt)

# order by descending AIC weight
aic.tast <- aic.tast[order(-aic.tast$wt),]

# calculate evidence ratio
aic.tast$ER <- max(aic.tast$wt) / aic.tast$wt

# check aic table
aic.tast[,3:6] <- round(aic.tast[,3:6], 3)
aic.tast
#    mod     pred   aicc  delta    wt      ER
# 9    9      dev 25.797  0.000 0.436   1.000
# 5    5    perim 28.007  2.210 0.144   3.019
# 14  14 humanmod 28.158  2.361 0.134   3.256
# 11  11     open 28.161  2.363 0.134   3.260
# 12  12   popden 29.654  3.857 0.063   6.879
# 7    7      imp 30.467  4.670 0.042  10.329
# 8    8   forest 32.800  7.002 0.013  33.155
# 6    6   island 33.536  7.738 0.009  47.901
# 1    1     null 33.683  7.885 0.008  51.559
# 2    2   areaha 34.796  8.999 0.005  89.951
# 4    4      age 35.438  9.640 0.004 123.985
# 13  13  povrate 36.030 10.233 0.003 166.760
# 10  10     tree 36.036 10.239 0.003 167.221
# 3    3    shape 36.060 10.263 0.003 169.283

# open threw a warning, but looks like it's 
# because there are no y=1 observations above open
# = 10, so fitted probabilities at higher open
# are effectively 0.
plot(dx.tast$open, jitter(dx.tast$y))
coef(summary(list.tast$open))
range(fitted(list.tast$open))

# get residuals
tast.res <- sapply(list.tast, residuals, type="pearson")

# run permutation-based Moran's I across models
set.seed(123)
moran.tast <- t(apply(tast.res, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))

moran.tast <- round(moran.tast, 3)
moran.tast <- as.data.frame(moran.tast)
moran.tast <- moran.tast[aic.tast$pred,]
moran.tast
#               I     p
# dev      -0.172 0.878
# perim    -0.113 0.699
# humanmod -0.151 0.807
# open     -0.198 0.919
# popden   -0.083 0.598
# imp      -0.078 0.549
# forest   -0.037 0.413
# island   -0.021 0.368
# null      0.098 0.128
# areaha    0.059 0.190
# age       0.058 0.188
# povrate   0.116 0.118
# tree      0.097 0.126
# shape     0.089 0.146

# check coefficients
summary(list.tast$dev)$coefficients
summary(list.tast$perim)$coefficients
summary(list.tast$humanmod)$coefficients
summary(list.tast$open)$coefficients

##############################################################################


##############################################################################

# make some figures
# rich: none
# simp: none
# pele: perim, dev
# sihi: age, open
# tast: dev, perim, humanmod, open

# figure 2: correlation plot
# lower triangle is |r|, scaled by magnitude
# and colored by direction.
# show any correlation with |r| > 0.3
# upper triangle labels "blocks" of thematically
# associated predictors
dmat <- de[,-1]
rownames(dmat) <- de$site
rmat <- cor(dmat, use="pairwise.complete.obs")
n <- ncol(rmat)

vars <- colnames(rmat)
vars2 <- c(
    "Area", "Shape", "Age", "%Imp.\nperim.",
    "Isol.", "Imperv.\ncover",
    "Forest\ncover","Dev.\ncover",
    "Tree\ncover","Open","Pop.\nDens.",
    "Poverty\nrate","Human\nMod.")

jpeg("rv-figure-02.jpeg", width=9, height=6.5,
     units="in", res=500)
par(mar=c(0.1, 0.1, 0.1, 0.1), lend=1,
    las=1, bty="n")
plot(NA, xlim=c(0.5, n+0.5),
     ylim=c(n+0.5, 0.5),
     xaxt="n", yaxt="n",
     bty="n", asp=(6.5/9),
     xlab="", ylab="")
for(i in 1:n){
    for(j in 1:n){
        if(i > j){
            rr <- rmat[i,j]
            text(j, i,
                 sprintf("%.2f", rr),
                 cex = 0.35 + 1.6 * abs(rr),
                 col=ifelse(rr > 0.3, "blue",
                            ifelse(rr < -0.3, "red", "grey70")))
        }#if
    }#j
}#i
xleft <- 0.5+(0):(n-1)
xright <- 1:n + 0.5
ybottom <- xleft
ytop <- xright
rect(xleft, ybottom, xright, ytop, lwd=2,
     col="grey65")
text(1:n, 1:n, vars2, font=3)
segments(xleft, ytop, xleft, n+0.5, col="grey70")
segments(0.5, seq(2.5, n+0.5, by=1),
         1:12+0.5, seq(2.5, n+0.5, by=1), col="grey70")
segments(1.5, 0.3, 5.5, 4.3, lwd=2)
segments(6.5, 5.3, 10.5, 9.3, lwd=2)
segments(11.5, 10.3, 13.5, 12.3, lwd=2)
SRT <- -36
text(4.2, 1.6, "Island\nbiogeography-like\npredictors",
     cex=1.5, srt=SRT)
text(8.8, 7, "Land cover", cex=1.6, srt=SRT)
text(12.8, 10.7, "Socioeconomic\nfactors", 
     cex=1.5, srt=SRT)
dev.off()


################################################

# figure 3: scatterplots of best supported models
# 3 rows x 2 columns
# pele vs. perim
# pele vs. dev
# sihi vs. age
# sihi vs. open
# tast vs. dev
# tast vs. perim or humanmod

n <- 200
px1 <- seq(min(de$perim, na.rm=TRUE), max(de$perim, na.rm=TRUE), length=n)
px2 <- seq(min(de$dev, na.rm=TRUE), max(de$dev, na.rm=TRUE), length=n)
px3 <- seq(min(de$age, na.rm=TRUE), max(de$age, na.rm=TRUE), length=n)
px4 <- seq(min(de$open, na.rm=TRUE), max(de$open, na.rm=TRUE), length=n)
px5 <- px2
px6 <- seq(min(de$humanmod, na.rm=TRUE), max(de$humanmod, na.rm=TRUE), length=n)

prx1 <- data.frame(perim=px1)
prx2 <- data.frame(dev=px2)
prx3 <- data.frame(age=px3)
prx4 <- data.frame(open=px4)
prx5 <- data.frame(dev=px5)
prx6 <- data.frame(humanmod=px6)

pred1 <- predict(list.pele$perim,    newdata=prx1,
                 se.fit=TRUE, interval="confidence")
pred2 <- predict(list.pele$dev,      newdata=prx2,
                 se.fit=TRUE, interval="confidence")
pred3 <- predict(list.sihi$age,      newdata=prx3,
                 se.fit=TRUE, type="link")
pred4 <- predict(list.sihi$open,     newdata=prx4,
                 se.fit=TRUE, type="link")
pred5 <- predict(list.tast$dev,      newdata=prx5,
                 se.fit=TRUE, type="link")
pred6 <- predict(list.tast$humanmod, newdata=prx6,
                 se.fit=TRUE, type="link")

prx1$lo <- exp(pred1$fit[,2])-1
prx1$mn <- exp(pred1$fit[,1])-1
prx1$up <- exp(pred1$fit[,3])-1

prx2$lo <- exp(pred2$fit[,2])-1
prx2$mn <- exp(pred2$fit[,1])-1
prx2$up <- exp(pred2$fit[,3])-1

prx3$lo <- plogis(qnorm(0.025, pred3$fit, pred3$se.fit))
prx3$mn <- plogis(qnorm(0.5, pred3$fit, pred3$se.fit))
prx3$up <- plogis(qnorm(0.975, pred3$fit, pred3$se.fit))

prx4$lo <- plogis(qnorm(0.025, pred4$fit, pred4$se.fit))
prx4$mn <- plogis(qnorm(0.5, pred4$fit, pred4$se.fit))
prx4$up <- plogis(qnorm(0.975, pred4$fit, pred4$se.fit))

prx5$lo <- plogis(qnorm(0.025, pred5$fit, pred5$se.fit))
prx5$mn <- plogis(qnorm(0.5, pred5$fit, pred5$se.fit))
prx5$up <- plogis(qnorm(0.975, pred5$fit, pred5$se.fit))

prx6$lo <- plogis(qnorm(0.025, pred6$fit, pred6$se.fit))
prx6$mn <- plogis(qnorm(0.5, pred6$fit, pred6$se.fit))
prx6$up <- plogis(qnorm(0.975, pred6$fit, pred6$se.fit))

poly.col <- "grey80"
use.pch <- 16
pcex <- 1.2

# silhouettes:
# https://www.phylopic.org/images/78d30905-8878-41ee-a612-700c1bf09ae9/peromyscus-leucopus
# https://www.phylopic.org/images/6239499e-114f-4828-bb0c-643351c70b9c/tamias-striatus
# https://www.phylopic.org/images/81930c02-5f26-43f7-9c19-e9831e780e53/sigmodon-hispidus

library(rsvg)
library(png)

rsvg_png("images/pele.svg", "images/pele.png")
rsvg_png("images/sihi.svg", "images/sihi.png")
rsvg_png("images/tamias.svg", "images/tamias.png")

pele.img   <- readPNG("images/pele.png")
sihi.img   <- readPNG("images/sihi.png")
tamias.img <- readPNG("images/tamias.png")

add.silhouette <- function(img, x, y, width = 0.12,
                           adj = c(0.5, 0.5)) {
    
    # Plot region in user coordinates
    usr <- par("usr")
    x.range <- usr[2] - usr[1]
    y.range <- usr[4] - usr[3]
    
    # Plot region in inches
    pin <- par("pin")
    
    # Image aspect ratio: width / height
    img.aspect <- dim(img)[2] / dim(img)[1]
    
    # Desired width in user coordinates
    w.user <- width * x.range
    
    # Convert that width to physical inches
    w.in <- w.user / x.range * pin[1]
    
    # Height in physical inches preserving image aspect
    h.in <- w.in / img.aspect
    
    # Convert height back to user coordinates
    h.user <- h.in / pin[2] * y.range
    
    # Position according to adj
    xleft   <- x - adj[1] * w.user
    xright  <- xleft + w.user
    ybottom <- y - adj[2] * h.user
    ytop    <- ybottom + h.user
    
    rasterImage(img,
                xleft = xleft,
                ybottom = ybottom,
                xright = xright,
                ytop = ytop)
}


yline <- 3.75

jpeg("rv-figure-03.jpg", width=9.4, height=10.8,
     units="in", res=600)
par(mfrow=c(3,2), mar=c(5.1, 6.1, 1.1, 1.1), 
    bty="n", lend=1, las=1,
    oma=c(0, 1, 2, 0),
    cex.axis=2.1, cex.lab=2.1,
    xpd=NA)
plot(de$perim, dx.pele$y,
     xlim=c(0, 50), ylim=c(0, 40),
     xlab="Site perimeter imperviousness (%)",
     ylab="")
title(main="A", adj=0, font.main=2, cex.main=2)
polygon(x=c(px1, rev(px1)), y=c(prx1$lo, rev(prx1$up)),
        border=NA, col=poly.col)
points(px1, prx1$mn, type="l", lwd=3)
points(de$perim, dx.pele$y, pch=use.pch, cex=pcex)
add.silhouette(pele.img, x=5, y=40, width=0.3, adj=c(0,1))
title(ylab=expression(italic("P. leucopus")~pop.~density~(n/ha)), line=yline)

plot(de$dev, dx.pele$y,
     xlim=c(0, 100), ylim=c(0, 40),
     xlab="Developed cover (%)",
     ylab="")
title(main="B", adj=0, font.main=2, cex.main=2)
polygon(x=c(px2, rev(px2)), y=c(prx2$lo, rev(prx2$up)),
        border=NA, col=poly.col)
points(px2, prx2$mn, type="l", lwd=3)
points(de$dev, dx.pele$y, pch=use.pch, cex=pcex)
add.silhouette(pele.img, x=0, y=40, width=0.3, adj=c(0,1))
title(ylab=expression(italic("P. leucopus")~pop.~density~(n/ha)), line=yline)

plot(de$age, jitter(dx.sihi$y, amount=0.02),
     xlim=c(0, 90), type="n",
     xlab="Site age (years)",
     ylab="")
title(main="C", adj=0, font.main=2, cex.main=2)
polygon(x=c(px3, rev(px3)), y=c(prx3$lo, rev(prx3$up)),
        border=NA, col=poly.col)
points(px3, prx3$mn, type="l", lwd=3)
points(de$age, jitter(dx.sihi$y, amount=0.02), pch=use.pch, cex=pcex)
add.silhouette(sihi.img, x=90, y=1, width=0.3, adj=c(1,1))
title(ylab=expression(italic("S. hispidus")~detection~probability), line=yline)

plot(de$open, jitter(dx.sihi$y, amount=0.02),
     xlim=c(0, 55), type="n",
     xlab="Open land cover (%)",
     ylab="")
title(main="D", adj=0, font.main=2, cex.main=2)
polygon(x=c(px4, rev(px4)), y=c(prx4$lo, rev(prx4$up)),
        border=NA, col=poly.col)
points(px4, prx4$mn, type="l", lwd=3)
points(de$open, jitter(dx.sihi$y, amount=0.02), pch=use.pch, cex=pcex)
add.silhouette(sihi.img, x=55, y=0.95, width=0.3, adj=c(1,1))
title(ylab=expression(italic("S. hispidus")~detection~probability), line=yline)

plot(de$dev, jitter(dx.tast$y, amount=0.02),
     xlim=c(0, 100), type="n",
     xlab="Developed cover (%)",
     ylab="")
title(main="E", adj=0, font.main=2, cex.main=2)
polygon(x=c(px5, rev(px5)), y=c(prx5$lo, rev(prx5$up)),
        border=NA, col=poly.col)
points(px5, prx5$mn, type="l", lwd=3)
points(de$dev, jitter(dx.tast$y, amount=0.02), pch=use.pch, cex=pcex)
add.silhouette(tamias.img, x=0, y=1, width=0.25, adj=c(0,1))
title(ylab=expression(italic("T. striatus")~detection~probability), line=yline)

plot(de$humanmod, jitter(dx.tast$y, amount=0.02),
     xlim=c(0, 1), type="n",
     xlab="Human modification index (unitless)",
     ylab="")
title(main="F", adj=0, font.main=2, cex.main=2)
polygon(x=c(px6, rev(px6)), y=c(prx6$lo, rev(prx6$up)),
        border=NA, col=poly.col)
points(px6, prx6$mn, type="l", lwd=3)
points(de$humanmod, jitter(dx.tast$y, amount=0.02), pch=use.pch, cex=pcex)
add.silhouette(tamias.img, x=0, y=1, width=0.25, adj=c(0,1))
title(ylab=expression(italic("T. striatus")~detection~probability), line=yline)
dev.off()
