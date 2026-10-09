# script to support manuscript by B. Casement, L. Lopez,
# K. Nestor, and N. Green entitled, "Urban small mammal
# communities are different, not depauperate: 
# species-specific responses along an urban-to-rural
# gradient" submitted to Journal of Mammalogy

# requires data files
# (in folder /data)
# "dat-mammal.csv"
# "dat-explanatory.csv"
# "weights_15.7km.rds"

# requires image files
# for figure 3
# (in folder /images)
# "pele.png"
# "sihi.png"
# "tamias.png"

# produces outputs:
## "rv-figure-02.jpg"
## "rv-figure-03.jpg"
## "rv-figure-04.jpg"

# requires packages:
library(MuMIn)
library(spdep)
library(png)
library(betapart)
library(vegan)

# get response variables
dy <- read.csv("data/dat-mammal.csv")

# get explanatory variables
dx <- read.csv("data/dat-explanatory.csv")

# get spatial weights
w <- readRDS("data/weights_15.7km.rds")

###########################################################
#                                                         #
# preliminary steps                                       #
#                                                         #
########################################################### 

# must be TRUE
all(dx$site == dy$site)

# sites -> rownames
rownames(dy) <- dy$site
rownames(dx) <- dx$site
dy$site <- NULL
dx$site <- NULL

# make a data frame for each response
dy$rich <- apply(dy, 1, function(x){sum(x > 0)})
dx.rich <- data.frame(y=dy$rich, dx)
dx.pele <- data.frame(y=dy$pele, dx)
dx.sihi <- data.frame(y=dy$sihi, dx)
dx.tast <- data.frame(y=dy$tast, dx)

# set up lists to hold model objects
list.rich <- vector("list", length=ncol(dx.rich))
names(list.rich) <- c("null", names(dx.rich)[-1])
list.pele <- list.rich
list.sihi <- list.rich
list.tast <- list.rich

###########################################################
#                                                         #
# make figure 2                                           #
#                                                         #
########################################################### 

rmat <- cor(dx, use="pairwise.complete.obs")
n <- ncol(rmat)

vars <- colnames(rmat)
vars2 <- c(
    "Area", "Shape", "Age", "Perim.\nImperv.",
    "Isol.", "Imperv.\ncover",
    "Forest\ncover","Dev.\ncover",
    "Tree\ncover","Open","Pop.\nDens.",
    "Poverty\nrate","Human\nMod.")

jpeg("figures/figure-02-revision.jpg",
     width=9, height=6.5,
     units="in", res=800)
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
text(3.9, 1.75, "Colonization\nand extinction",
     cex=1.5, srt=SRT)
text(8.8, 7, "Land cover", cex=1.6, srt=SRT)
text(12.8, 10.7, "Socioeconomic\nfactors", 
     cex=1.5, srt=SRT)
dev.off()

# TIFF version
tiff("figures/figure-02-revision.tif",
     width=9, height=6.5,
     units="in", res=800,
     compression="lzw")
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
text(3.9, 1.75, "Colonization\nand extinction",
     cex=1.5, srt=SRT)
text(8.8, 7, "Land cover", cex=1.6, srt=SRT)
text(12.8, 10.7, "Socioeconomic\nfactors", 
     cex=1.5, srt=SRT)
dev.off()

###########################################################
#                                                         #
# analysis 1: Species richness                            #
#                                                         #
###########################################################

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
list.rich[[14]] <- glm(y~human,    data=dx.rich, family=poisson)

# check convergence
# must all  be TRUE:
sapply(list.rich, function(x){x$converged})

# find AICc and calculate AICc weights
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

# round
aic.rich[,3:ncol(aic.rich)] <- round(aic.rich[,3:ncol(aic.rich)], 3)

# inspect
aic.rich

# get residuals
rich.res <- sapply(list.rich, residuals, type="pearson")

# run permutation-based Moran's I across models
set.seed(123)
moran.rich <- t(apply(rich.res, 2, function(x) {
    z <- moran.mc(x, listw = w, nsim = 9999)
    
    c(I = unname(z$statistic),
      p = z$p.value)
}))
moran.rich <- as.data.frame(round(moran.rich, 3))
moran.rich

###########################################################
#                                                         #
# analysis 2: Peromyscus (PELE) population density        #
#                                                         #
###########################################################

# use log(y + 1) transform
dx.pele$z <- log(dx.pele$y + 1)

# fit models
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
list.pele[[14]] <- lm(z~human,    data=dx.pele)

# AIC model selection
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

# inspect
aic.pele

# aic.pele at hypothesis level
hyp.weights <- function(aic.tab) {
    
    H1 <- c("areaha", "shape", "age", "perim", "island")
    H2 <- c("imp", "forest", "dev", "tree", "open")
    H3 <- c("popden", "povrate", "humanmod")
    
    hyp <- rep(NA_character_, nrow(aic.tab))
    hyp[aic.tab$pred %in% H1] <- "H1"
    hyp[aic.tab$pred %in% H2] <- "H2"
    hyp[aic.tab$pred %in% H3] <- "H3"
    
    rel.like <- exp(-0.5 * aic.tab$delta)
    
    hyp.like <- tapply(rel.like, hyp, mean)
    
    hyp.like / sum(hyp.like)
}
hyp.pele <- hyp.weights(aic.pele)
hyp.pele

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


###########################################################
#                                                         #
# analysis 3: Sigmodon hispidus (SIHI) detection          #
#                                                         #
###########################################################

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
list.sihi[[14]] <- glm(y~human,    data=dx.sihi, family=binomial)

# check convergence
# must all  be TRUE:
sapply(list.sihi, function(x){x$converged})

# AICc model selection
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

# hypothesis level weights
hyp.sihi <- hyp.weights(aic.sihi)
hyp.sihi

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

###########################################################
#                                                         #
# analysis 4: Tamias striatus (TAST) detection            #
#                                                         #
###########################################################

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
list.tast[[14]] <- glm(y~human   , data=dx.tast, family=binomial)

# check convergence
# must all  be TRUE:
sapply(list.tast, function(x){x$converged})

# AICc model selection
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

# model for open threw a warning, but looks like it's
# because there are no y=1 observations above open
# = 10, so fitted probabilities at higher open
# are effectively 0.
plot(dx.tast$open, jitter(dx.tast$y))
coef(summary(list.tast$open))
range(fitted(list.tast$open))

# hypothesis level weights
hyp.tast <- hyp.weights(aic.tast)
hyp.tast
max(hyp.tast) / hyp.tast

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

###########################################################
#                                                         #
# analyses 3-5 summary: hypothesis support table and      #
#   coefficient table (Tables 3 and SD5)                  #
###########################################################

# hypothesis support table
hyp.table <- data.frame(
    hyp=1:3,
    pele=round(hyp.pele, 3), 
    sihi=round(hyp.sihi, 3), 
    tast=round(hyp.tast, 3)
)
hyp.table

# table of coefficients
sum.list <- list(
    list.pele$perim,
    list.pele$dev,
    list.sihi$age,
    list.sihi$open,
    list.tast$dev,
    list.tast$perim)

# check coefficients
lapply(sum.list, function(x){summary(x)$coefficients})

# compile coefficients for Table SD5
coef.tables <- lapply(sum.list, function(m) {
    ci <- round(confint(m),3)
    
    data.frame(
        Parameter = names(coef(m)),
        Estimate = coef(m),
        Lower95 = ci[, 1],
        Upper95 = ci[, 2],
        row.names = NULL
    )
})
do.call(rbind, coef.tables)

###########################################################
#                                                         #
# make figure 3                                           #
#                                                         #
###########################################################

# figure 3: scatterplots of best supported models
# 3 rows x 2 columns
# pele vs. perim
# pele vs. dev
# sihi vs. age
# sihi vs. open
# tast vs. dev
# tast vs. humanmod

# n for smooth curve
n <- 200

# set up x values for prediction
px1 <- seq(min(dx$perim, na.rm=TRUE), max(dx$perim, na.rm=TRUE), length=n)
px2 <- seq(min(dx$dev, na.rm=TRUE), max(dx$dev, na.rm=TRUE), length=n)
px3 <- seq(min(dx$age, na.rm=TRUE), max(dx$age, na.rm=TRUE), length=n)
px4 <- seq(min(dx$open, na.rm=TRUE), max(dx$open, na.rm=TRUE), length=n)
px5 <- px2
px6 <- px1

# put in data frames
prx1 <- data.frame(perim=px1)
prx2 <- data.frame(dev=px2)
prx3 <- data.frame(age=px3)
prx4 <- data.frame(open=px4)
prx5 <- data.frame(dev=px5)
prx6 <- data.frame(perim=px6)

# calculate predictions
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
pred6 <- predict(list.tast$perim, newdata=prx6,
                 se.fit=TRUE, type="link")

# backtransform from link to data scale
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

# silhouettes:
# download the SVG version and convert to png
# using rsvg::rsvg_png() or another utility
# https://www.phylopic.org/images/78d30905-8878-41ee-a612-700c1bf09ae9/peromyscus-leucopus
# https://www.phylopic.org/images/6239499e-114f-4828-bb0c-643351c70b9c/tamias-striatus
# https://www.phylopic.org/images/81930c02-5f26-43f7-9c19-e9831e780e53/sigmodon-hispidus

pele.img   <- readPNG("images/pele.png")
sihi.img   <- readPNG("images/sihi.png")
tamias.img <- readPNG("images/tamias.png")

# helper function to add silhouettes to figure 3
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
}#function

# constant to adjust line of y axis labels
yline <- 3.75

# color for 95% CI polygon
poly.col <- "grey80"

# point character for figure 3
use.pch <- 16

# point size for figure 3
pcex <- 1.2


jpeg("figures/figure-03-revision.jpg",
     width=7.3, height=8.4,
     units="in", res=800)
par(mfrow=c(3,2), mar=c(5.1, 6.1, 1.1, 1.1), 
    bty="n", lend=1, las=1,
    oma=c(0, 1, 2, 0),
    cex.axis=1.7, cex.lab=1.7,
    xpd=NA)
# panel A: pele pop den vs. perimeter imperviousness
plot(dx$perim, dx.pele$y,
     xlim=c(0, 50), ylim=c(0, 40),
     xlab="Perimeter imperviousness (%)",
     ylab="")
title(main="A", adj=0, font.main=2, cex.main=2)
polygon(x=c(px1, rev(px1)), y=c(prx1$lo, rev(prx1$up)),
        border=NA, col=poly.col)
points(px1, prx1$mn, type="l", lwd=3)
points(dx$perim, dx.pele$y, pch=use.pch, cex=pcex)
add.silhouette(pele.img, x=5, y=40, width=0.3, adj=c(0,1))
title(ylab=expression(italic("P. leucopus")~pop.~density~(n/ha)),
      line=yline)
# panel B: pele pop den vs. developed land cover
plot(dx$dev, dx.pele$y,
     xlim=c(0, 100), ylim=c(0, 40),
     xlab="Developed cover (%)",
     ylab="")
title(main="B", adj=0, font.main=2, cex.main=2)
polygon(x=c(px2, rev(px2)), y=c(prx2$lo, rev(prx2$up)),
        border=NA, col=poly.col)
points(px2, prx2$mn, type="l", lwd=3)
points(dx$dev, dx.pele$y, pch=use.pch, cex=pcex)
add.silhouette(pele.img, x=0, y=40, width=0.3, adj=c(0,1))
title(ylab=expression(italic("P. leucopus")~pop.~density~(n/ha)),
      line=yline)
# panel C: sihi detection vs. site age
plot(dx$age, dx.sihi$y, type="n",
     xlim=c(0, 90),
     xlab="Site age (years)",
     ylab="")
title(main="C", adj=0, font.main=2, cex.main=2)
polygon(x=c(px3, rev(px3)), y=c(prx3$lo, rev(prx3$up)),
        border=NA, col=poly.col)
points(px3, prx3$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$age, jitter(dx.sihi$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(sihi.img, x=90, y=1, width=0.3, adj=c(1,1))
title(ylab=expression(italic("S. hispidus")~detection~probability),
      line=yline)
# panel D: sihi detection vs. open land cover
plot(dx$open, dx.sihi$y, type="n",
     xlim=c(0, 55), 
     xlab="Open land cover (%)",
     ylab="")
title(main="D", adj=0, font.main=2, cex.main=2)
polygon(x=c(px4, rev(px4)), y=c(prx4$lo, rev(prx4$up)),
        border=NA, col=poly.col)
points(px4, prx4$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$open, jitter(dx.sihi$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(sihi.img, x=55, y=0.95,
               width=0.3, adj=c(1,1))
title(ylab=expression(italic("S. hispidus")~detection~probability),
      line=yline)
# panel E: tast detection vs. developed land cover
plot(dx$dev, dx.tast$y, type="n",
     xlim=c(0, 100), 
     xlab="Developed cover (%)",
     ylab="")
title(main="E", adj=0, font.main=2, cex.main=2)
polygon(x=c(px5, rev(px5)), y=c(prx5$lo, rev(prx5$up)),
        border=NA, col=poly.col)
points(px5, prx5$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$dev, jitter(dx.tast$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(tamias.img, x=0, y=1, width=0.25, adj=c(0,1))
title(ylab=expression(italic("T. striatus")~detection~probability),
      line=yline)
# panel F: tast detection vs. perim
plot(dx$perim, dx.tast$y, type="n",
     xlim=c(0, 50), 
     xlab="Perimeter imperviousness (%)",
     ylab="")
title(main="F", adj=0, font.main=2, cex.main=2)
polygon(x=c(px6, rev(px6)), y=c(prx6$lo, rev(prx6$up)),
        border=NA, col=poly.col)
points(px6, prx6$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$perim, jitter(dx.tast$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(tamias.img, x=30, y=0.95,
               width=0.25, adj=c(0,1))
title(ylab=expression(italic("T. striatus")~detection~probability),
      line=yline)
dev.off()

# TIFF version for revisions
tiff("figures/figure-03-revision.tif",
     width=7.3, height=8.4,
     units="in", res=800,
     compression="lzw")
par(mfrow=c(3,2), mar=c(5.1, 6.1, 1.1, 1.1), 
    bty="n", lend=1, las=1,
    oma=c(0, 1, 2, 0),
    cex.axis=1.7, cex.lab=1.7,
    xpd=NA)
# panel A: pele pop den vs. perimeter imperviousness
plot(dx$perim, dx.pele$y,
     xlim=c(0, 50), ylim=c(0, 40),
     xlab="Perimeter imperviousness (%)",
     ylab="")
title(main="A", adj=0, font.main=2, cex.main=2)
polygon(x=c(px1, rev(px1)), y=c(prx1$lo, rev(prx1$up)),
        border=NA, col=poly.col)
points(px1, prx1$mn, type="l", lwd=3)
points(dx$perim, dx.pele$y, pch=use.pch, cex=pcex)
add.silhouette(pele.img, x=5, y=40, width=0.3, adj=c(0,1))
title(ylab=expression(italic("P. leucopus")~pop.~density~(n/ha)),
      line=yline)
# panel B: pele pop den vs. developed land cover
plot(dx$dev, dx.pele$y,
     xlim=c(0, 100), ylim=c(0, 40),
     xlab="Developed cover (%)",
     ylab="")
title(main="B", adj=0, font.main=2, cex.main=2)
polygon(x=c(px2, rev(px2)), y=c(prx2$lo, rev(prx2$up)),
        border=NA, col=poly.col)
points(px2, prx2$mn, type="l", lwd=3)
points(dx$dev, dx.pele$y, pch=use.pch, cex=pcex)
add.silhouette(pele.img, x=0, y=40, width=0.3, adj=c(0,1))
title(ylab=expression(italic("P. leucopus")~pop.~density~(n/ha)),
      line=yline)
# panel C: sihi detection vs. site age
plot(dx$age, dx.sihi$y, type="n",
     xlim=c(0, 90),
     xlab="Site age (years)",
     ylab="")
title(main="C", adj=0, font.main=2, cex.main=2)
polygon(x=c(px3, rev(px3)), y=c(prx3$lo, rev(prx3$up)),
        border=NA, col=poly.col)
points(px3, prx3$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$age, jitter(dx.sihi$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(sihi.img, x=90, y=1, width=0.3, adj=c(1,1))
title(ylab=expression(italic("S. hispidus")~detection~probability),
      line=yline)
# panel D: sihi detection vs. open land cover
plot(dx$open, dx.sihi$y, type="n",
     xlim=c(0, 55), 
     xlab="Open land cover (%)",
     ylab="")
title(main="D", adj=0, font.main=2, cex.main=2)
polygon(x=c(px4, rev(px4)), y=c(prx4$lo, rev(prx4$up)),
        border=NA, col=poly.col)
points(px4, prx4$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$open, jitter(dx.sihi$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(sihi.img, x=55, y=0.95,
               width=0.3, adj=c(1,1))
title(ylab=expression(italic("S. hispidus")~detection~probability),
      line=yline)
# panel E: tast detection vs. developed land cover
plot(dx$dev, dx.tast$y, type="n",
     xlim=c(0, 100), 
     xlab="Developed cover (%)",
     ylab="")
title(main="E", adj=0, font.main=2, cex.main=2)
polygon(x=c(px5, rev(px5)), y=c(prx5$lo, rev(prx5$up)),
        border=NA, col=poly.col)
points(px5, prx5$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$dev, jitter(dx.tast$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(tamias.img, x=0, y=1, width=0.25, adj=c(0,1))
title(ylab=expression(italic("T. striatus")~detection~probability),
      line=yline)
# panel F: tast detection vs. perim
plot(dx$perim, dx.tast$y, type="n",
     xlim=c(0, 50), 
     xlab="Perimeter imperviousness (%)",
     ylab="")
title(main="F", adj=0, font.main=2, cex.main=2)
polygon(x=c(px6, rev(px6)), y=c(prx6$lo, rev(prx6$up)),
        border=NA, col=poly.col)
points(px6, prx6$mn, type="l", lwd=3)
## set random number seed for jitter
set.seed(123)
points(dx$perim, jitter(dx.tast$y, amount=0.02),
       pch=use.pch, cex=pcex)
add.silhouette(tamias.img, x=30, y=0.95,
               width=0.25, adj=c(0,1))
title(ylab=expression(italic("T. striatus")~detection~probability),
      line=yline)
dev.off()


###########################################################
#                                                         #
# analysis 5: Beta diversity and species turnover         #
#                                                         #
###########################################################

# set up for beta diversity analysis
# need presence/absence only
db <- dy
db$rich <- NULL
db$pele <- as.numeric(db$pele > 0)

# must be TRUE
all(rownames(db) == rownames(dx))

# Pairwise beta diversity
beta.pw <- beta.pair(db, index.family = "sorensen")

# Components
bsor <- beta.pw$beta.sor   # total beta diversity
bsim <- beta.pw$beta.sim   # turnover
bsne <- beta.pw$beta.sne   # nestedness

# distance matrix based on human modification index
# pairwise differnces
dhuman <- dist(setNames(dx$human, rownames(dx)))

# mantel tests for correlation between
# dhuman and beta diversity components
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


###########################################################
#                                                         #
# make figure 4                                           #
#                                                         #
###########################################################

# make sure in same order
all(rownames(db)==rownames(dx))

# site names
sites <- dx$site

# reorder copies of the datasets
dx4 <- dx[order(dx$human),]
dy4 <- db[rownames(dx4),]

# tranpose to make figure 4
dm <- as.matrix(dy4)
dm <- t(dm)

# make vector of 1s and 0s to put species in 
# order from along HMI
sporder <- apply(dm, 1, paste, collapse="")
sporder <- sort(sporder, decreasing=TRUE)
dm <- dm[names(sporder),]

# data frame to see mean HMI of sites where
# each species is present
sp <- data.frame(spp=rownames(dm), mn=NA, sd=NA,
                 n=NA)
hmi <- dx4$human
names(hmi) <- rownames(dx4)

for(i in 1:nrow(sp)){
    isites <- colnames(dm)[which(dm[i,] == 1)]
    sp$mn[i] <- mean(hmi[isites])
    sp$sd[i] <- sd(hmi[isites])
    sp$n[i] <- length(isites)
}
sp <- sp[order(sp$mn),]
sp
dm2 <- dm[sp$spp,]

# y coordinates for plot
y1 <- 1:nrow(dm2)
y2 <- rev(y1)

# x coordinates for plot (hmi)
h <- hmi[colnames(dm)]
hmi.diff <- apply(dm, 1, function(x) {
    mean(h[x == 1]) - mean(h[x == 0])
})

dm3 <- dm2[names(sort(hmi.diff, decreasing=TRUE)),]

splabs <- expression(
    italic(Peromyscus~gossypinus),
    italic(Reithrodontomys~humulis),
    italic(Sigmodon~hispidus),
    italic(Ochrotomys~nuttalli),
    italic(Microtus~pinetorum),
    italic(Blarina~carolinensis),
    italic(Glaucomys~volans),
    italic(Sciurus~carolinensis),
    italic(Peromyscus~leucopus),
    italic(Tamias~striatus),
    italic(Rattus~norvegicus),
    italic(Mus~musculus)
)


jpeg("figures/figure-04-revision.jpg", width=10.25,
     height=5.3, units="in", res=800)
set.seed(123)
par(mar=c(5.1, 16.1, 1.1, 1.1), bty="n",
    cex.axis=1.3, cex.lab=1.3, lend=1,
    xpd=NA, las=1)
plot(hmi, rep(1,23), type="n",
     xlim=c(0.2, 1), ylim=c(0.5, 12.5),
     yaxt="n", ylab="",
     xlab="Human modification index (unitless)")
segments(0.2, 1:12, 1, 1:12, lty=2, col="grey80")
for(i in 1:12){
    points(hmi, 
           jitter(rep(y2[i], length(hmi)), amount=0.09),
           pch=16,
           col=ifelse(dm3[i,]==0, "transparent", "black"),
           cex=1.5)
}
axis(side=2, at=y1, labels=splabs)
dev.off()


tiff("figures/figure-04-revision.tif", width=10.25,
     height=5.3, units="in", res=800)
set.seed(123)
par(mar=c(5.1, 16.1, 1.1, 1.1), bty="n",
    cex.axis=1.3, cex.lab=1.3, lend=1,
    xpd=NA, las=1)
plot(hmi, rep(1,23), type="n",
     xlim=c(0.2, 1), ylim=c(0.5, 12.5),
     yaxt="n", ylab="",
     xlab="Human modification index (unitless)")
segments(0.2, 1:12, 1, 1:12, lty=2, col="grey80")
for(i in 1:12){
    points(hmi, 
           jitter(rep(y2[i], length(hmi)), amount=0.09),
           pch=16,
           col=ifelse(dm3[i,]==0, "transparent", "black"),
           cex=1.5)
}
axis(side=2, at=y1, labels=splabs)
dev.off()




jpeg("figures/figure-04-revision.jpg", width=10.25,
     height=5.3, units="in", res=800)
set.seed(123)
par(mar=c(5.1, 16.1, 1.1, 1.1), bty="n",
    cex.axis=1.3, cex.lab=1.3, lend=1,
    xpd=NA, las=1)
plot(hmi, rep(1,23), type="n",
     xlim=c(0.2, 1), ylim=c(0.5, 12.5),
     yaxt="n", ylab="",
     xlab="Human modification index (unitless)")
segments(0.2, 1:12, 1, 1:12, lty=2, col="grey80")
for(i in 1:12){
  points(hmi, 
         jitter(rep(y2[i], length(hmi)), amount=0.09),
         pch=16,
         col=ifelse(dm3[i,]==0, "transparent", "black"),
         cex=1.5)
}
axis(side=2, at=y1, labels=splabs)
dev.off()


tiff("figures/figure-04-revision.tif", width=10.25,
     height=5.3, units="in", res=800,
     compression="lzw")
set.seed(123)
par(mar=c(5.1, 16.1, 1.1, 1.1), bty="n",
    cex.axis=1.3, cex.lab=1.3, lend=1,
    xpd=NA, las=1)
plot(hmi, rep(1,23), type="n",
     xlim=c(0.2, 1), ylim=c(0.5, 12.5),
     yaxt="n", ylab="",
     xlab="Human modification index (unitless)")
segments(0.2, 1:12, 1, 1:12, lty=2, col="grey80")
for(i in 1:12){
  points(hmi, 
         jitter(rep(y2[i], length(hmi)), amount=0.09),
         pch=16,
         col=ifelse(dm3[i,]==0, "transparent", "black"),
         cex=1.5)
}
axis(side=2, at=y1, labels=splabs)
dev.off()


# end script!