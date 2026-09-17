# script to visualize species turnover

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

sites <- dx$site
dx$site <- NULL
dx <- sapply(dx, function(x){1*(x>0)})
rownames(dx) <- sites

dz <- dz[order(dz$human),]
dx <- dx[dz$site,]

dm <- as.matrix(dx)
dm <- t(dm)

sporder <- apply(dm, 1, paste, collapse="")
sporder <- sort(sporder, decreasing=TRUE)
sporder

dm <- dm[names(sporder),]

sp <- data.frame(spp=rownames(dm), mn=NA, sd=NA,
                 n=NA)
hmi <- dz$human
names(hmi) <- rownames(dz)

for(i in 1:nrow(sp)){
    isites <- colnames(dm)[which(dm[i,] == 1)]
    sp$mn[i] <- mean(hmi[isites])
    sp$sd[i] <- sd(hmi[isites])
    sp$n[i] <- length(isites)
}
sp <- sp[order(sp$mn),]
sp
dm2 <- dm[sp$spp,]

y1 <- 1:nrow(dm2)
y2 <- rev(y1)

par(mar=c(5.1, 5.1, 1.1,1.1), bty="n")
plot(hmi, rep(1,23), type="n",
     xlim=c(0.2, 1), ylim=c(0.5, 12.5),
     yaxt="n")
for(i in 1:12){
    points(hmi, rep(y2[i], length(hmi)), pch=22,
           bg=ifelse(dm2[i,]==0, "white", "black"),
           cex=1.5)
}
axis(side=2, at=y1, labels=rownames(dm2)[y2], las=1)


par(mar=c(5.1, 5.1, 1.1,1.1), bty="n")
plot(hmi, rep(1,23), type="n",
     xlim=c(0.2, 1), ylim=c(0.5, 12.5),
     yaxt="n", ylab="")
for(i in 1:12){
    points(hmi, rep(y2[i], length(hmi)), pch=15,
           col=ifelse(dm2[i,]==0, "transparent", "black"),
           cex=1.5)
}
axis(side=2, at=y1, labels=rownames(dm2)[y2], las=1)

h <- hmi[colnames(dm)]
hmi.diff <- apply(dm, 1, function(x) {
    mean(h[x == 1]) - mean(h[x == 0])
})

sort(hmi.diff)

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


jpeg("rv-figure-04.jpg", width=10.25,
     height=5.3, units="in", res=500)
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
