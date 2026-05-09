# univariate analysis of urban rodent abundance 
# and presence/absence data
library(vegan)
library(dunn.test)
library(pROC)

# final explanatory variable matrix @ 100 m
in.name <- "dat_explanatory_0100_2026-04-27.csv"
ex <- read.csv(in.name)

# final response variable matrix - by species
in.name <- "dat_response_species_2026-04-27.csv"
dx <- read.csv(in.name)

# final response variable matrix - by trait group
in.name <- "dat_response_traits_2026-04-27.csv"
kx <- read.csv(in.name)

# put site names on rownames
rownames(dx) <- dx$site
rownames(kx) <- kx$site
rownames(ex) <- ex$site

# should all be TRUE: 
# means dataframes in same row order by site!
rownames(dx) == rownames(kx)
rownames(ex) == rownames(kx)

# delete site columns
dx$site <- NULL
ex$site <- NULL
kx$site <- NULL

# univariate vs. urbanization
## species richness
kruskal.test(dx$rich~ex$type)

## esg richness
kruskal.test(kx$rich~ex$type)

## pele pop density
kruskal.test(dx$pele~ex$type)
dunn.test(x=dx$pele, g=ex$type)

## sihi pop density
kruskal.test(dx$sihi~ex$type)

## tast pop density
kruskal.test(dx$tast~ex$type)
dunn.test(x=dx$tast, g=ex$type)

## esg2 individual density (specialists)
kruskal.test(kx$k2~ex$type)

## esg6 individual density (generalists)
kruskal.test(kx$k6~ex$type)

## hypotheses 2, 3, 4
## general explanatory matrix
use.cols <- c("areaha", "shape", "age", "island", "pct_imp_perim",
              "tree", "imp", "forest", "dev", "open",
              "popden", "income", "povrate")
dex <- ex[,use.cols]
nmodels <- length(use.cols)+1

#####################################################################

# species richness - Poisson GLM
d01 <- data.frame(y=dx$rich, dex)

list.spprich <- vector("list", length=nmodels)
list.spprich[[1]] <- glm(y~1,             data=d01, family=poisson)
list.spprich[[2]] <- glm(y~areaha,        data=d01, family=poisson)
list.spprich[[3]] <- glm(y~shape,         data=d01, family=poisson)
list.spprich[[4]] <- glm(y~age,           data=d01, family=poisson)
list.spprich[[5]] <- glm(y~island,        data=d01, family=poisson)
list.spprich[[6]] <- glm(y~pct_imp_perim, data=d01, family=poisson)
list.spprich[[7]] <- glm(y~tree,          data=d01, family=poisson)
list.spprich[[8]] <- glm(y~imp,           data=d01, family=poisson)
list.spprich[[9]] <- glm(y~forest,        data=d01, family=poisson)
list.spprich[[10]] <- glm(y~dev,           data=d01, family=poisson)
list.spprich[[11]] <- glm(y~open,          data=d01, family=poisson)
list.spprich[[12]] <- glm(y~popden,        data=d01, family=poisson)
list.spprich[[13]] <- glm(y~income,        data=d01, family=poisson)
list.spprich[[14]] <- glm(y~povrate,       data=d01, family=poisson)

# write omnibus tests out in text file
out.name <- "res_01_species_richness.txt"

# clear output file
cat("", file=out.name)

# fill output file in loop
for(i in 2:nmodels){
    lrt <- anova(list.spprich[[1]], list.spprich[[i]], test="Chisq")
    iname <- use.cols[i-1]
    block <- paste(
        sprintf("=== %s vs Null ===", iname),
        paste(capture.output(lrt), collapse = "\n"),
        sep = "\n"
    )
    cat(block, "\n\n", file = out.name, append = TRUE)
}#i

#####################################################################

# PELE pop density - linear regression
d02 <- data.frame(y=dx$pele, dex)
list.pelepop <- vector("list", length=nmodels)

list.pelepop[[01]] <- lm(y~1,             data=d02)
list.pelepop[[02]] <- lm(y~areaha,        data=d02)
list.pelepop[[03]] <- lm(y~shape,         data=d02)
list.pelepop[[04]] <- lm(y~age,           data=d02)
list.pelepop[[05]] <- lm(y~island,        data=d02)
list.pelepop[[06]] <- lm(y~pct_imp_perim, data=d02)
list.pelepop[[07]] <- lm(y~tree,          data=d02)
list.pelepop[[08]] <- lm(y~imp,           data=d02)
list.pelepop[[09]] <- lm(y~forest,        data=d02)
list.pelepop[[10]] <- lm(y~dev,           data=d02)
list.pelepop[[11]] <- lm(y~open,          data=d02)
list.pelepop[[12]] <- lm(y~popden,        data=d02)
list.pelepop[[13]] <- lm(y~income,        data=d02)
list.pelepop[[14]] <- lm(y~povrate,       data=d02)

# write omnibus tests out in text file
out.name <- "res_02_pele_pop_den.txt"

# clear output file
cat("", file=out.name)

# fill output file in loop

for(i in 2:nmodels){
    cat("===================\n", file = out.name, append=TRUE)
    
    omnibus <- anova(list.pelepop[[1]], list.pelepop[[i]], test="F")
    anovatab <- anova(list.pelepop[[i]])
    
    iname <- use.cols[i-1]
    block <- paste(
        sprintf("=== %s omnibus test ===", iname),
        paste(capture.output(omnibus), collapse = "\n"),
        "",
        sprintf("=== %s ANOVA table ===", iname),
        paste(capture.output(anovatab), collapse="\n"),
        sep="\n"
    )
    cat(block, "\n\n", file = out.name, append = TRUE)
}#i

#####################################################################

# ESG6 pop density - linear regression
d03 <- data.frame(y=kx$k6, dex)
list.k6pop <- vector("list", length=nmodels)

list.k6pop[[01]] <- lm(y~1,             data=d03)
list.k6pop[[02]] <- lm(y~areaha,        data=d03)
list.k6pop[[03]] <- lm(y~shape,         data=d03)
list.k6pop[[04]] <- lm(y~age,           data=d03)
list.k6pop[[05]] <- lm(y~island,        data=d03)
list.k6pop[[06]] <- lm(y~pct_imp_perim, data=d03)
list.k6pop[[07]] <- lm(y~tree,          data=d03)
list.k6pop[[08]] <- lm(y~imp,           data=d03)
list.k6pop[[09]] <- lm(y~forest,        data=d03)
list.k6pop[[10]] <- lm(y~dev,           data=d03)
list.k6pop[[11]] <- lm(y~open,          data=d03)
list.k6pop[[12]] <- lm(y~popden,        data=d03)
list.k6pop[[13]] <- lm(y~income,        data=d03)
list.k6pop[[14]] <- lm(y~povrate,       data=d03)

# write omnibus tests out in text file
out.name <- "res_03_esg6_pop_den.txt"

# clear output file
cat("", file=out.name)

# fill output file in loop

for(i in 2:nmodels){
    cat("===================\n", file = out.name, append=TRUE)
    
    omnibus <- anova(list.k6pop[[1]], list.k6pop[[i]], test="F")
    anovatab <- anova(list.k6pop[[i]])
    
    iname <- use.cols[i-1]
    block <- paste(
        sprintf("=== %s omnibus test ===", iname),
        paste(capture.output(omnibus), collapse = "\n"),
        "",
        sprintf("=== %s ANOVA table ===", iname),
        paste(capture.output(anovatab), collapse="\n"),
        sep="\n"
    )
    cat(block, "\n\n", file = out.name, append = TRUE)
}#i


#####################################################################

# SIHI presence/absence - logistic regression
d04 <- data.frame(y=dx$sihi, dex)
d04$y <- ifelse(d04$y > 0, 1, 0)
list.sihi <- vector("list", length=nmodels)

list.sihi[[01]] <- glm(y~1,             data=d04, family=binomial)
list.sihi[[02]] <- glm(y~areaha,        data=d04, family=binomial)
list.sihi[[03]] <- glm(y~shape,         data=d04, family=binomial)
list.sihi[[04]] <- glm(y~age,           data=d04, family=binomial)
list.sihi[[05]] <- glm(y~island,        data=d04, family=binomial)
list.sihi[[06]] <- glm(y~pct_imp_perim, data=d04, family=binomial)
list.sihi[[07]] <- glm(y~tree,          data=d04, family=binomial)
list.sihi[[08]] <- glm(y~imp,           data=d04, family=binomial)
list.sihi[[09]] <- glm(y~forest,        data=d04, family=binomial)
list.sihi[[10]] <- glm(y~dev,           data=d04, family=binomial)
list.sihi[[11]] <- glm(y~open,          data=d04, family=binomial)
list.sihi[[12]] <- glm(y~popden,        data=d04, family=binomial)
list.sihi[[13]] <- glm(y~income,        data=d04, family=binomial)
list.sihi[[14]] <- glm(y~povrate,       data=d04, family=binomial)

# write omnibus tests out in text file
out.name <- "res_04_sihi_pres.txt"

# clear output file
cat("", file=out.name)

# fill output file in loop

for(i in 2:nmodels){
    cat("===================\n", file = out.name, append=TRUE)
    
    omnibus <- anova(list.sihi[[1]], list.sihi[[i]], test="Chisq")

    p1 <- predict(list.sihi[[i]], type="response")
    roc1 <- roc(d04$y~p1, plot=FALSE)$auc
    
    iname <- use.cols[i-1]
    block <- paste(
        sprintf("=== %s omnibus test ===", iname),
        paste(capture.output(omnibus), collapse = "\n"),
        "",
        sprintf("=== %s AUC ===", iname),
        paste(capture.output(roc1), collapse="\n"),
        sep="\n"
    )
    cat(block, "\n\n", file = out.name, append = TRUE)
}#i

#####################################################################

# ESG2 pop density - linear regression
d05 <- data.frame(y=kx$k2, dex)
list.k2pop <- vector("list", length=nmodels)

list.k2pop[[01]] <- lm(y~1,             data=d05)
list.k2pop[[02]] <- lm(y~areaha,        data=d05)
list.k2pop[[03]] <- lm(y~shape,         data=d05)
list.k2pop[[04]] <- lm(y~age,           data=d05)
list.k2pop[[05]] <- lm(y~island,        data=d05)
list.k2pop[[06]] <- lm(y~pct_imp_perim, data=d05)
list.k2pop[[07]] <- lm(y~tree,          data=d05)
list.k2pop[[08]] <- lm(y~imp,           data=d05)
list.k2pop[[09]] <- lm(y~forest,        data=d05)
list.k2pop[[10]] <- lm(y~dev,           data=d05)
list.k2pop[[11]] <- lm(y~open,          data=d05)
list.k2pop[[12]] <- lm(y~popden,        data=d05)
list.k2pop[[13]] <- lm(y~income,        data=d05)
list.k2pop[[14]] <- lm(y~povrate,       data=d05)

# write omnibus tests out in text file
out.name <- "res_05_esg2_pop_den.txt"

# clear output file
cat("", file=out.name)

# fill output file in loop

for(i in 2:nmodels){
    cat("===================\n", file = out.name, append=TRUE)
    
    omnibus <- anova(list.k2pop[[1]], list.k2pop[[i]], test="F")
    anovatab <- anova(list.k2pop[[i]])
    
    iname <- use.cols[i-1]
    block <- paste(
        sprintf("=== %s omnibus test ===", iname),
        paste(capture.output(omnibus), collapse = "\n"),
        "",
        sprintf("=== %s ANOVA table ===", iname),
        paste(capture.output(anovatab), collapse="\n"),
        sep="\n"
    )
    cat(block, "\n\n", file = out.name, append = TRUE)
}#i

#####################################################################

# TAST presence/absence - logistic regression
d06 <- data.frame(y=dx$tast, dex)
d06$y <- ifelse(d06$y > 0, 1, 0)
list.tast <- vector("list", length=nmodels)

list.tast[[01]] <- glm(y~1,             data=d06, family=binomial)
list.tast[[02]] <- glm(y~areaha,        data=d06, family=binomial)
list.tast[[03]] <- glm(y~shape,         data=d06, family=binomial)
list.tast[[04]] <- glm(y~age,           data=d06, family=binomial)
list.tast[[05]] <- glm(y~island,        data=d06, family=binomial)
list.tast[[06]] <- glm(y~pct_imp_perim, data=d06, family=binomial)
list.tast[[07]] <- glm(y~tree,          data=d06, family=binomial)
list.tast[[08]] <- glm(y~imp,           data=d06, family=binomial)
list.tast[[09]] <- glm(y~forest,        data=d06, family=binomial)
list.tast[[10]] <- glm(y~dev,           data=d06, family=binomial)
list.tast[[11]] <- glm(y~open,          data=d06, family=binomial)
list.tast[[12]] <- glm(y~popden,        data=d06, family=binomial)
list.tast[[13]] <- glm(y~income,        data=d06, family=binomial)
list.tast[[14]] <- glm(y~povrate,       data=d06, family=binomial)

# write omnibus tests out in text file
out.name <- "res_06_tast_pres.txt"

# clear output file
cat("", file=out.name)

# fill output file in loop

for(i in 2:nmodels){
    cat("===================\n", file = out.name, append=TRUE)
    
    omnibus <- anova(list.tast[[1]], list.tast[[i]], test="Chisq")
    
    p1 <- predict(list.tast[[i]], type="response")
    roc1 <- roc(d06$y~p1, plot=FALSE)$auc
    
    iname <- use.cols[i-1]
    block <- paste(
        sprintf("=== %s omnibus test ===", iname),
        paste(capture.output(omnibus), collapse = "\n"),
        "",
        sprintf("=== %s AUC ===", iname),
        paste(capture.output(roc1), collapse="\n"),
        sep="\n"
    )
    cat(block, "\n\n", file = out.name, append = TRUE)
}#i

# generate model output tables for supplement
prep.out <- function(Q, yname){
    z <- lapply(Q, function(x){summary(x)$coefficients})
    z <- do.call(rbind, z)
    #z <- round(z, 3)
    z <- formatC(z, format="f", digits=3)
    z <- gsub("-0.000", "zebra", z)
    z <- gsub("0.000", "<0.001", z)
    z <- gsub("zebra", ">-0.001", z)
    z <- as.data.frame(z)
    z$param <- rownames(z)
    z$param[grep("Inter", z$param)] <- "Intercept"
    z$y <- yname
    rownames(z) <- NULL
    z <- z[,c(6, 5,1:4)]
    
    out <- z[,1:2]
    out$estse <- paste(z[,3], "±", z[,4])
    out$test <- z[,5]
    out$pval <- z[,6]
    return(out)
}#prep.out

out1 <- prep.out(list.spprich, "rich")
out2 <- prep.out(list.pelepop, "pele")
out3 <- prep.out(list.k6pop, "esg6")
out4 <- prep.out(list.sihi, "sihi")
out5 <- prep.out(list.k2pop, "esg2")
out6 <- prep.out(list.tast, "tast")

write.csv(out1, "model_coefs_rich.csv")
write.csv(out2, "model_coefs_pele.csv")
write.csv(out3, "model_coefs_esg6.csv")
write.csv(out4, "model_coefs_sihi.csv")
write.csv(out5, "model_coefs_esg2.csv")
write.csv(out6, "model_coefs_tast.csv")


#####################################################################

# produce figures 

# figure 2: 2x2 boxplots y vs. urbanization
## 2a: rich
## 2b: pele
## 2c: sihi
## 2d: tast
fx <- dx[,c("rich", "pele", "sihi", "tast")]
fx$type <- ex$type
    
cols <- c("olivedrab", "steelblue", "orchid3")

fig.name <- "figure_02.jpg"
jpeg(fig.name, width=8, height=8, units="in", res=600)
par(mfrow=c(2,2), bty="n", las=1, lend=1, mar=c(5.1, 5.1, 1.7, 1.1),
    cex.axis=1.4, cex.lab=1.5, oma=c(0, 0, 0.2, 0), xpd=NA)
boxplot(rich~type, data=fx, col=cols, xaxt="n",
        ylab="Species richness", xlab="Site type", ylim=c(0, 9))
axis(side=1, at=1:3, labels=c("Rural", "Suburban", "Urban"))
text(1:3, c(5.27, 7.34, 5.30), "a", cex=1.5)
title(main="a", adj=0, cex.main=2.5)
boxplot(pele~type, data=fx, col=cols, xaxt="n",
        ylab=expression(italic("P. leucopus")~population~density~(italic(n)/ha)),
        xlab="Site type", ylim=c(0, 20))
title(main="b", adj=0, cex.main=2.5)
axis(side=1, at=1:3, labels=c("Rural", "Suburban", "Urban"))
text(1:3, c(6.85, 17.23, 12.39), c("a", "b", "ab"), cex=1.5)
boxplot(sihi~type, data=fx, col=cols, xaxt="n",
        ylab=expression(italic("S. hispidus")~population~density~(italic(n)/ha)),
        xlab="Site type", ylim=c(0, 2))
title(main="c", adj=0, cex.main=2.5)
axis(side=1, at=1:3, labels=c("Rural", "Suburban", "Urban"))
text(1:3, c(1.19, 1.69, 0.04)+0.05, c("a", "a", "b"), cex=1.5)
boxplot(tast~type, data=fx, col=cols, xaxt="n", ylim=c(0, 3),
        ylab=expression(italic("T. striatus")~population~density~(italic(n)/ha)),
        xlab="Site type", yaxt="n")
axis(side=1, at=1:3, labels=c("Rural", "Suburban", "Urban"))
axis(side=2, at=0:3)
title(main="d", adj=0, cex.main=2.5)
text(1:3, c(0.14, 0.88, 1.42), c("a", "b", "b"), cex=1.5)
dev.off()





# figure 3: nmds (other script)

# figure 4: 2 x 2 scatterplots
## 4a: pele vs. pip; tast vs. pip (right axis)
## 4b: pele vs. dev; tast vs. dev (right axis)
## 4c: pele vs. open or forest
## 4d: sihi vs. age

# predictions
n <- 100
# pele vs. pip
px1 <- seq(min(ex$pct_imp_perim), max(ex$pct_imp_perim), length=n)
prx1 <- data.frame(pct_imp_perim=px1)
pred <- predict(list.pelepop[[6]], newdata=prx1, se.fit=TRUE)
prx1$lo <- qnorm(0.025, pred$fit, pred$se.fit)
prx1$up <- qnorm(0.975, pred$fit, pred$se.fit)
prx1$mm <- qnorm(0.5,   pred$fit, pred$se.fit)

# tast vs. pip
px2 <- px1
prx2 <- data.frame(pct_imp_perim=px2)
pred <- predict(list.tast[[6]], newdata=prx1, se.fit=TRUE, type="link")
prx2$lo <- plogis(qnorm(0.025, pred$fit, pred$se.fit))
prx2$up <- plogis(qnorm(0.975, pred$fit, pred$se.fit))
prx2$mm <- plogis(qnorm(0.5,   pred$fit, pred$se.fit))

# pele vs. dev
px3 <- seq(min(ex$dev), max(ex$dev), length=n)
prx3 <- data.frame(dev=px3)
pred <- predict(list.pelepop[[10]], newdata=prx3, se.fit=TRUE)
prx3$lo <- qnorm(0.025, pred$fit, pred$se.fit)
prx3$up <- qnorm(0.975, pred$fit, pred$se.fit)
prx3$mm <- qnorm(0.5,   pred$fit, pred$se.fit)

# tast vs. dev
px4 <- px3
prx4 <- data.frame(dev=px4)
pred <- predict(list.tast[[10]], newdata=prx4, se.fit=TRUE, type="link")
prx4$lo <- plogis(qnorm(0.025, pred$fit, pred$se.fit))
prx4$up <- plogis(qnorm(0.975, pred$fit, pred$se.fit))
prx4$mm <- plogis(qnorm(0.5,   pred$fit, pred$se.fit))

# pele vs. forest
px5a <- seq(min(ex$forest), max(ex$forest), length=n)
prx5a <- data.frame(forest=px5a)
pred <- predict(list.pelepop[[9]], newdata=prx5a, se.fit=TRUE)
prx5a$lo <- qnorm(0.025, pred$fit, pred$se.fit)
prx5a$up <- qnorm(0.975, pred$fit, pred$se.fit)
prx5a$mm <- qnorm(0.5,   pred$fit, pred$se.fit)

# pele vs. open
px5b <- seq(min(ex$open), max(ex$open), length=n)
prx5b <- data.frame(open=px5b)
pred <- predict(list.pelepop[[11]], newdata=prx5b, se.fit=TRUE)
prx5b$lo <- qnorm(0.025, pred$fit, pred$se.fit)
prx5b$up <- qnorm(0.975, pred$fit, pred$se.fit)
prx5b$mm <- qnorm(0.5,   pred$fit, pred$se.fit)

# sihi vs. age
px6 <- seq(min(ex$age), max(ex$age), length=n)
prx6 <- data.frame(age=px6)
pred <- predict(list.sihi[[4]], newdata=prx6, se.fit=TRUE, type="link")
prx6$lo <- plogis(qnorm(0.025, pred$fit, pred$se.fit))
prx6$up <- plogis(qnorm(0.975, pred$fit, pred$se.fit))
prx6$mm <- plogis(qnorm(0.5,   pred$fit, pred$se.fit))


cols <- c("olivedrab", "steelblue", "orchid3")
names(cols) <- c("rural", "suburban", "urban")
ex$color <- cols[ex$type]
pchs <- c(16, 15, 17)
names(pchs) <- names(cols)
ex$pch <- pchs[ex$type]

ace.up <- 22
bdf.up <- 1.1

jpeg("figure_04.jpg", width=6.5, height=9, units="in",
     res=500)
par(mfrow=c(3,2), mar=c(5.1, 5.1, 2.1, 1.1), bty="n",
    lend=1, las=1, cex.axis=1.4, cex.lab=1.4)
plot(ex$pct_imp_perim, dx$pele, pch=ex$pch, col=ex$color, cex=1.3,
     xaxt="n", yaxt="n", ylim=c(0, 20),
     xlab="Perimeter imperviousness (%)",
     ylab=expression(italic("P. leucopus")~population~density~(italic(n)/ha)))
xx <- prx1$pct_imp_perim
uu <- prx1$up
ll <- prx1$lo
mm <- prx1$mm
polygon(x=c(xx, rev(xx)), y=c(ll, rev(uu)), border=NA, col="grey80")
points(xx, mm, type="l", lwd=4)
axis(side=1, at=0:5*10)
axis(side=2, at=0:4*5)
legend("topright", legend=c("Rural", "Suburban", "Urban"),
       pch=pchs, col=cols, cex=1.3, bty="n")
text(0, ace.up, "A", cex=3, font=2, adj=0, xpd=NA)
points(ex$pct_imp_perim, dx$pele, pch=ex$pch, col=ex$color, cex=1.3)

##############################################################################

plot(ex$pct_imp_perim, d06$y, pch=ex$pch, col=ex$color, cex=1.3,
     ylim=c(0, 1), xlim=c(0, 50), xaxt="n", yaxt="n",
     xlab="Perimeter imperviousness (%)",
     ylab=expression("Prob."~of~italic("T. striatus")~detection))
axis(side=1, at=0:5*10)
axis(side=2, at=0:5*0.2)
text(0, bdf.up, "B", cex=3, font=2, adj=0, xpd=NA)
xx <- prx2$pct_imp_perim
uu <- prx2$up
ll <- prx2$lo
mm <- prx2$mm
polygon(x=c(xx, rev(xx)), y=c(ll, rev(uu)), border=NA, col="grey80")
points(xx, mm, type="l", lwd=4)
points(ex$pct_imp_perim, d06$y, pch=ex$pch, col=ex$color, cex=1.3)

###############################################################################
plot(ex$dev, dx$pele, pch=ex$pch, col=ex$color, cex=1.3,
     xlim=c(0, 100), ylim=c(0, 20), xaxt="n", yaxt="n",
     xlab="Developed cover (%)",
     ylab=expression(italic("P. leucopus")~population~density~(italic(n)/ha)))
axis(side=1, at=0:5*20)
axis(side=2, at=0:4*5)
text(0, ace.up, "C", cex=3, font=2, adj=0, xpd=NA)
xx <- prx3$dev
uu <- prx3$up
ll <- prx3$lo
mm <- prx3$mm
polygon(x=c(xx, rev(xx)), y=c(ll, rev(uu)), border=NA, col="grey80")
points(xx, mm, type="l", lwd=4)
points(ex$dev, dx$pele, pch=ex$pch, col=ex$color, cex=1.3)

###############################################################################
plot(ex$dev, d06$y, pch=ex$pch, col=ex$color, cex=1.3,
     ylim=c(0, 1), xlim=c(0, 100), xaxt="n", yaxt="n",
     xlab="Developed cover (%)",
     ylab=expression("Prob."~of~italic("T. striatus")~detection))
axis(side=1, at=0:5*20)
axis(side=2, at=0:5*0.2)
text(0, bdf.up, "D", cex=3, font=2, adj=0, xpd=NA)
xx <- prx4$dev
uu <- prx4$up
ll <- prx4$lo
mm <- prx4$mm
polygon(x=c(xx, rev(xx)), y=c(ll, rev(uu)), border=NA, col="grey80")
points(xx, mm, type="l", lwd=4)
points(ex$dev, d06$y, pch=ex$pch, col=ex$color, cex=1.3)

###############################################################################
plot(ex$forest, dx$pele, pch=ex$pch, col=ex$color, cex=1.3,
     xlim=c(0, 70), ylim=c(0, 20), xaxt="n", yaxt="n",
     xlab="Forest cover (%)",
     ylab=expression(italic("P. leucopus")~population~density~(italic(n)/ha)))
axis(side=1, at=0:7*10)
axis(side=2, at=0:4*5)
text(0, ace.up, "E", cex=3, font=2, adj=0, xpd=NA)
xx <- prx5a$forest
uu <- prx5a$up
ll <- pmax(0, prx5a$lo)
mm <- pmax(0, prx5a$mm)
polygon(x=c(xx, rev(xx)), y=c(ll, rev(uu)), border=NA, col="grey80")
points(xx, mm, type="l", lwd=4)
points(ex$forest, dx$pele, pch=ex$pch, col=ex$color, cex=1.3)
###############################################################################
# plot(ex$open, dx$pele, pch=ex$pch, col=ex$color, cex=1.3,
#      xlim=c(0, 50), ylim=c(0, 20), xaxt="n", yaxt="n",
#      xlab="Open cover (%)",
#      ylab=expression(italic("P. leucopus")~population~density~(italic(n)/ha)))
# axis(side=1, at=0:5*10)
# axis(side=2, at=0:4*5)
# text(0, 21, "E", cex=3, font=2, adj=0, xpd=NA)
# xx <- prx5b$open
# uu <- prx5b$up
# ll <- pmax(0, prx5b$lo)
# mm <- pmax(0, prx5b$mm)
# polygon(x=c(xx, rev(xx)), y=c(ll, rev(uu)), border=NA, col="grey80")
# points(xx, mm, type="l", lwd=4)
# points(ex$open, dx$pele, pch=ex$pch, col=ex$color, cex=1.3)

###############################################################################
plot(ex$age, d04$y, pch=ex$pch, col=ex$color, cex=1.3,
     xlim=c(0, 90), ylim=c(0, 1), xaxt="n", yaxt="n",
     xlab="Fragment age (y)",
     ylab=expression("Prob."~of~italic("S. hispidus")~detection))
axis(side=1, at=0:3*30)
axis(side=2, at=0:5*0.2)
text(0, bdf.up, "F", cex=3, font=2, adj=0, xpd=NA)
xx <- prx6$age
uu <- prx6$up
ll <- prx6$lo
mm <- prx6$mm
polygon(x=c(xx, rev(xx)), y=c(ll, rev(uu)), border=NA, col="grey80")
points(xx, mm, type="l", lwd=4)
points(ex$age, d04$y, pch=ex$pch, col=ex$color, cex=1.3)
###############################################################################
dev.off()
