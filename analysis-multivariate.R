library(vegan)

# species-based response matrix
in.name <- "dat_response_species_2026-04-27.csv"
dx <- read.csv(in.name)

# functional group response matrix
in.name <- "dat_response_traits_2026-04-27.csv"
kx <- read.csv(in.name)

# final explanatory variable matrices
in.name <- "dat_explanatory_0100_2026-04-27.csv"
ex <- read.csv(in.name)

# fix rownames and make sure dataframes are lined up
rownames(dx) <- dx$site
rownames(kx) <- kx$site
rownames(ex) <- ex0100$site

dx$site <- NULL
kx$site <- NULL
ex$site <- NULL

all(rownames(dx) == rownames(kx))
all(rownames(dx) == rownames(ex))

###### NMDS ordination on population density by site type ####

# species population densities only

# species abundance matrix dxm
dxm <- dx[,1:12]
set.seed(42)
n1 <- metaMDS(dxm, try=100, trymax=500)  
n1

# set up for plot
urbs <- sort(unique(ex$type))
cols <- c("olivedrab", "steelblue", "orchid3")
names(cols) <- urbs
plotdf <- data.frame(color=cols[ex$type])

pchs <- c(21, 22, 24)
names(pchs) <- names(cols)
plotdf$pch <- pchs[ex$type]


# get site scores
xy1 <- scores(n1, "sites", choices=c(1,2))
## get species scores
z1 <- scores(n1, choices=1:2, display="species")

# limits for plot given species scores
apply(rbind(z1, xy1), 2, range)

# 2 letter abbreviations for species names
spnames <- paste0(toupper(substr(rownames(z1), 1, 1)),
                  substr(rownames(z1), 3, 3))

######  envfit for explanatory variables
names(ex)[which(names(ex) == "pct_imp_perim")] <- "pct.imp"

use.x <- c(
    # IBT-related group
    "areaha", "shape", "island", "age", "pct.imp",
    # urban indicator group (matches Sydney ms)
    "imp", "tree", "popden",
    # general landcover group
    "open", "forest", "dev",
    # site-level group
    "elev.mean", "db.mn", "lux_day", 
    "tempc_day",
    # soceioeconomic group
    "income", "povrate"
)
xx <- ex[,use.x]

set.seed(12345)
e1 <- envfit(n1, xx, permutations=1e4)

as.matrix(e1$vectors$arrows)[which(e1$vectors$pvals < 0.1),]

ze1 <- as.matrix(e1$vectors$arrows)[which(e1$vectors$pvals < 0.05),]
ze1



# scalar to move vector labels out a bit
zfac <- 0.11
zcol <- "black"
zcex <- 1.2
pt.lwd <- 1.1
pt.cex <- 1.3
zecol <- "blue"

# add a constant to the absolute value of each z coord
zadj <- function(x, ad){
    res <- x
    res[res<0] <- res[res<0]-ad
    res[res>0] <- res[res>0]+ad
    return(res)
}

fig.name <- "figure_03.jpg"
jpeg(fig.name, width=7, height=6, units="in", res=600)
par(mfrow=c(1,1), bty="n", lend=1, las=1, cex.axis=1.3, cex.lab=1.3,
    mar=c(5.1, 5.1, 2.1, 1.1), xpd=NA, cex.main=1.5)
plot(xy1[,1], xy1[,2], pch=plotdf$pch, bg=plotdf$col,
     cex=pt.cex, lwd=pt.lwd,
     xlab="NMDS axis 1", ylab="NMDS axis 2",
     xlim=c(-2, 3), xaxt="n",
     ylim=c(-2, 3), yaxt="n")
axis(side=1, at=seq(-2, 3, by=1))
axis(side=2, at=seq(-2, 3, by=1))
ordiellipse(n1, ex0100$type, draw="polygon",
            col=cols, border=NA, alpha=50)
points(xy1[,1], xy1[,2], pch=plotdf$pch, bg=plotdf$col,
       cex=pt.cex, lwd=pt.lwd)
segments(0, 0, z1[,1], z1[,2], col=zcol, lwd=2)
## add constants to z1 to move text labels
text(zadj(z1[,1],zfac)+c(0.05, 0, 0, 0, 0.12, 0, 0.05, 0, 0, -0.1, 0, 0),
     zadj(z1[,2],zfac)+c(0.05, 0, 0, 0.1, 0, 0, -0.03, 0, 0, 0, 0, 0.15),
     spnames,
     cex=zcex, col=zcol, font=3)
segments(0,0, ze1[,1], ze1[,2], col=zecol, lwd=2)
text(
    zadj(ze1[,1], zfac)+c(0.0,  -0.2, 0.0, 0.0),
    zadj(ze1[,2], zfac)+c(-0.1, -0.1, -0.1, 0.0),
    c("Shape", "Age", "Open", "Elev."), col=zecol, cex=zcex)
legend("topright", legend=c("Rural", "Suburban", "Urban"),
       pch=pchs, cex=pt.cex, pt.bg=cols, pt.lwd=pt.lwd, bty="n")
dev.off()

# species significant, but low R2
# functional groups NS
set.seed(123)
adonis2(dxm~ex$type, permutations=9999)
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = dxm ~ ex$type, permutations = 9999)
#         Df SumOfSqs      R2      F Pr(>F)  
#Model     2   0.9913 0.17786 2.1634 0.0244 *
#Residual 20   4.5823 0.82214                
#Total    22   5.5737 1.00000  

# posthoc test for permanovas on ABUNDANCES
set.seed(321)
flag1 <- which(ex$type %in% c("rural", "suburban"))
flag2 <- which(ex$type %in% c("rural", "urban"))
flag3 <- which(ex$type %in% c("suburban", "urban"))
adonis2(dxm[flag1,]~ex$type[flag1], permutations=9999)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = dxm[flag1, ] ~ ex$type[flag1], permutations = 9999)
#          Df SumOfSqs      R2       F Pr(>F)  
# Model     1   0.7094 0.17053 2.8782 0.0208 *
# Residual 14   3.4505 0.82947                
# Total    15   4.1599 1.00000    

set.seed(567)
adonis2(dxm[flag2,]~ex$type[flag2], permutations=9999)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = dxm[flag2, ] ~ ex$type[flag2], permutations = 9999)
#          Df SumOfSqs     R2      F  Pr(>F)  
# Model     1   0.5072 0.1306 1.9528  0.099 .
# Residual 13   3.3762 0.8694                
# Total    14   3.8834 1.0000 

set.seed(123)
adonis2(dxm[flag3,]~ex$type[flag3], permutations=9999)
#adonis2(formula = dxm[flag3, ] ~ ex$type[flag3], permutations = 9999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     1  0.25623 0.09878 1.4248 0.2139
#Residual 13  2.33786 0.90122              
#Total    14  2.59409 1.00000  

## what about rural vs. nonrural?
ex$urb2 <- ifelse(ex$type=="rural", "rural", "nonrural")
set.seed(321)
adonis2(dxm~ex$urb2, permutations=9999)
#adonis2(formula = dxm ~ ex$urb2, permutations = 9999)
#         Df SumOfSqs      R2      F  Pr(>F)  
#Model     1   0.7351 0.13189 3.1905 0.0108 *
#Residual 21   4.8385 0.86811                 
#Total    22   5.5737 1.00000   


# beta dispersion test
# multivariate analogue of Levene's test
Dmat <- vegdist(dxm, method="bray")
bd <- betadisper(Dmat, ex$type)
anova(bd)
set.seed(42)
bd_perm <- permutest(bd, permutations=9999)
bd_perm

plot(bd)
boxplot(bd, xlab = "Site type", ylab = "Distance to group centroid")
#"Homogeneity of multivariate dispersion did not differ among site types (betadisper permutation test, p = …), supporting interpretation of PERMANOVA results as differences in group centroids."
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 9999
#
#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.12350 0.061751 2.0251   9999  0.158
#Residuals 20 0.60985 0.030492                     

##############################################################################
####  beta dispersion analysis  ##############################################
##############################################################################

grp <- ex$type

bd <- betadisper(Dmat, grp)

set.seed(42)
bd_perm <- permutest(bd, permutations=9999)
bd_perm
plot(bd)
boxplot(bd, xlab = "Site type", ylab = "Distance to group centroid")

##############################################################################
####  beta diversity analysis  ###############################################
##############################################################################

# convert to presence absence
dxm_pa <- decostand(dxm, method="pa")

library(betapart)

beta_pa <- beta.pair(dxm_pa, index.family="sorensen")
beta_pa
#$beta.sor = total beta diversity
#$beta.sim = turnover (replacement)
#$beta.sne = nestedness-resultant component

# Compare pairwise distances by site-type contrasts

# Convert distance matrix to long format
beta_mat <- as.matrix(beta_pa$beta.sim)
beta_long <- as.data.frame(as.table(beta_mat), stringsAsFactors=FALSE)
names(beta_long) <- c("site1", "site2", "beta_val")

# keep unique pairs
beta_long <- beta_long[which(beta_long$site1 < beta_long$site2),]

site_df <- data.frame(
    site = rownames(dxm),
    type = ex$type,
    stringsAsFactors = FALSE
)

beta_long$type1 <- site_df$type[match(beta_long$site1, site_df$site)]
beta_long$type2 <- site_df$type[match(beta_long$site2, site_df$site)]

beta_long$contrast <- ifelse(
    beta_long$type1 == beta_long$type2,
    paste(beta_long$type1, "within", sep = "_"),
    paste(
        pmin(beta_long$type1, beta_long$type2),
        pmax(beta_long$type1, beta_long$type2),
        sep = "_vs_"
    )#paste
)#ifelse

ag.turnover <- aggregate(
    beta_val ~ contrast,
    data = beta_long,
    FUN = function(x) {round(c(n = length(x),
                               median = median(x),
                               lo=quantile(x,0.25),up=quantile(x,0.75)),2)}
)#aggregate

### nestedness
beta_mat <- as.matrix(beta_pa$beta.sne)
beta_long <- as.data.frame(as.table(beta_mat), stringsAsFactors=FALSE)
names(beta_long) <- c("site1", "site2", "beta_val")

# keep unique pairs
beta_long <- beta_long[which(beta_long$site1 < beta_long$site2),]

site_df <- data.frame(
    site = rownames(dxm),
    type = ex0100$type,
    stringsAsFactors = FALSE
)

beta_long$type1 <- site_df$type[match(beta_long$site1, site_df$site)]
beta_long$type2 <- site_df$type[match(beta_long$site2, site_df$site)]

beta_long$contrast <- ifelse(
    beta_long$type1 == beta_long$type2,
    paste(beta_long$type1, "within", sep = "_"),
    paste(
        pmin(beta_long$type1, beta_long$type2),
        pmax(beta_long$type1, beta_long$type2),
        sep = "_vs_"
    )#paste
)#ifelse

ag.nested <- aggregate(
    beta_val ~ contrast,
    data = beta_long,
    FUN = function(x) {round(c(n = length(x),
                               median = median(x),
                               lo=quantile(x,0.25),up=quantile(x,0.75)),2)}
)#aggregate


### total beta diversity
beta_mat <- as.matrix(beta_pa$beta.sor)
beta_long <- as.data.frame(as.table(beta_mat), stringsAsFactors=FALSE)
names(beta_long) <- c("site1", "site2", "beta_val")

# keep unique pairs
beta_long <- beta_long[which(beta_long$site1 < beta_long$site2),]

site_df <- data.frame(
    site = rownames(dxm),
    type = ex0100$type,
    stringsAsFactors = FALSE
)

beta_long$type1 <- site_df$type[match(beta_long$site1, site_df$site)]
beta_long$type2 <- site_df$type[match(beta_long$site2, site_df$site)]

beta_long$contrast <- ifelse(
    beta_long$type1 == beta_long$type2,
    paste(beta_long$type1, "within", sep = "_"),
    paste(
        pmin(beta_long$type1, beta_long$type2),
        pmax(beta_long$type1, beta_long$type2),
        sep = "_vs_"
    )#paste
)#ifelse

ag.total <- aggregate(
    beta_val ~ contrast,
    data = beta_long,
    FUN = function(x) {round(c(n = length(x),
                               median = median(x),
                               lo=quantile(x,0.25),up=quantile(x,0.75)),2)}
)#aggregate

ag.total[c(1,2,4,3,5,6),]
#            contrast beta_val.n beta_val.median beta_val.lo.25% beta_val.up.75%
# 1 rural_vs_suburban      64.00            0.67            0.50            0.76
# 2    rural_vs_urban      56.00            0.63            0.42            0.83
# 4 suburban_vs_urban      56.00            0.50            0.33            0.60
# 3      rural_within      28.00            0.63            0.33            1.00
# 5   suburban_within      28.00            0.50            0.40            0.60
# 6      urban_within      21.00            0.50            0.33            0.60

ag.turnover[c(1,2,4,3,5,6),]
#            contrast beta_val.n beta_val.median beta_val.lo.25% beta_val.up.75%
# 1 rural_vs_suburban      64.00            0.50            0.00            0.50
# 2    rural_vs_urban      56.00            0.50            0.00            0.81
# 4 suburban_vs_urban      56.00            0.00            0.00            0.50
# 3      rural_within      28.00            0.38            0.00            1.00
# 5   suburban_within      28.00            0.00            0.00            0.50
# 6      urban_within      21.00            0.00            0.00            0.33

ag.nested[c(1,2,4,3,5,6),]
#            contrast beta_val.n beta_val.median beta_val.lo.25% beta_val.up.75%
# 1 rural_vs_suburban      64.00            0.17            0.00            0.33
# 2    rural_vs_urban      56.00            0.09            0.00            0.33
# 4 suburban_vs_urban      56.00            0.21            0.08            0.43
# 3      rural_within      28.00            0.19            0.00            0.33
# 5   suburban_within      28.00            0.30            0.10            0.45
# 6      urban_within      21.00            0.33            0.10            0.50

########################################

# beta diversity again, but with rural vs. nonrural

# Convert distance matrix to long format
beta_mat <- as.matrix(beta_pa$beta.sim)
beta_long <- as.data.frame(as.table(beta_mat), stringsAsFactors=FALSE)
names(beta_long) <- c("site1", "site2", "beta_val")

# keep unique pairs
beta_long <- beta_long[which(beta_long$site1 < beta_long$site2),]

site_df <- data.frame(
    site = rownames(dxm),
    type = ex$urb2,
    stringsAsFactors = FALSE
)

beta_long$type1 <- site_df$type[match(beta_long$site1, site_df$site)]
beta_long$type2 <- site_df$type[match(beta_long$site2, site_df$site)]

beta_long$contrast <- ifelse(
    beta_long$type1 == beta_long$type2,
    paste(beta_long$type1, "within", sep = "_"),
    paste(
        pmin(beta_long$type1, beta_long$type2),
        pmax(beta_long$type1, beta_long$type2),
        sep = "_vs_"
    )#paste
)#ifelse

ag.turnover <- aggregate(
    beta_val ~ contrast,
    data = beta_long,
    FUN = function(x) {round(c(n = length(x),
                               median = median(x),
                               lo=quantile(x,0.25),up=quantile(x,0.75)),2)}
)#aggregate

### nestedness
beta_mat <- as.matrix(beta_pa$beta.sne)
beta_long <- as.data.frame(as.table(beta_mat), stringsAsFactors=FALSE)
names(beta_long) <- c("site1", "site2", "beta_val")

# keep unique pairs
beta_long <- beta_long[which(beta_long$site1 < beta_long$site2),]

site_df <- data.frame(
    site = rownames(dxm),
    type = ex$urb2,
    stringsAsFactors = FALSE
)

beta_long$type1 <- site_df$type[match(beta_long$site1, site_df$site)]
beta_long$type2 <- site_df$type[match(beta_long$site2, site_df$site)]

beta_long$contrast <- ifelse(
    beta_long$type1 == beta_long$type2,
    paste(beta_long$type1, "within", sep = "_"),
    paste(
        pmin(beta_long$type1, beta_long$type2),
        pmax(beta_long$type1, beta_long$type2),
        sep = "_vs_"
    )#paste
)#ifelse

ag.nested <- aggregate(
    beta_val ~ contrast,
    data = beta_long,
    FUN = function(x) {round(c(n = length(x),
                               median = median(x),
                               lo=quantile(x,0.25),up=quantile(x,0.75)),2)}
)#aggregate


### total beta diversity
beta_mat <- as.matrix(beta_pa$beta.sor)
beta_long <- as.data.frame(as.table(beta_mat), stringsAsFactors=FALSE)
names(beta_long) <- c("site1", "site2", "beta_val")

# keep unique pairs
beta_long <- beta_long[which(beta_long$site1 < beta_long$site2),]

site_df <- data.frame(
    site = rownames(dxm),
    type = ex$urb2,
    stringsAsFactors = FALSE
)

beta_long$type1 <- site_df$type[match(beta_long$site1, site_df$site)]
beta_long$type2 <- site_df$type[match(beta_long$site2, site_df$site)]

beta_long$contrast <- ifelse(
    beta_long$type1 == beta_long$type2,
    paste(beta_long$type1, "within", sep = "_"),
    paste(
        pmin(beta_long$type1, beta_long$type2),
        pmax(beta_long$type1, beta_long$type2),
        sep = "_vs_"
    )#paste
)#ifelse

ag.total <- aggregate(
    beta_val ~ contrast,
    data = beta_long,
    FUN = function(x) {round(c(n = length(x),
                               median = median(x),
                               lo=quantile(x,0.25),up=quantile(x,0.75)),2)}
)#aggregate


ag.total
#           contrast beta_val.n beta_val.median beta_val.lo.25% beta_val.up.75%
#1 nonrural_vs_rural     120.00            0.67            0.50            0.78
#2   nonrural_within     105.00            0.50            0.33            0.60
#3      rural_within      28.00            0.63            0.33            1.00
ag.turnover
#           contrast beta_val.n beta_val.median beta_val.lo.25% beta_val.up.75%
#1 nonrural_vs_rural     120.00            0.50            0.00            0.69
#2   nonrural_within     105.00            0.00            0.00            0.50
#3      rural_within      28.00            0.38            0.00            1.00
ag.nested
#           contrast beta_val.n beta_val.median beta_val.lo.25% beta_val.up.75%
#1 nonrural_vs_rural     120.00            0.10            0.00            0.33
#2   nonrural_within     105.00            0.27            0.10            0.50
#3      rural_within      28.00            0.19            0.00            0.33

# end script!