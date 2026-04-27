# univariate analysis of urban rodent abundance 
# and presence/absence data
library(vegan)
library(dunn.test)

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

# species richness - Poisson GLM
d01 <- data.frame(y=dx$rich, dex)

mod01.01 <- glm(y~areaha,        data=d01, family=poisson)
mod01.02 <- glm(y~shape,         data=d01, family=poisson)
mod01.03 <- glm(y~age,           data=d01, family=poisson)
mod01.04 <- glm(y~island,        data=d01, family=poisson)
mod01.05 <- glm(y~pct_imp_perim, data=d01, family=poisson)
mod01.06 <- glm(y~tree,          data=d01, family=poisson)
mod01.07 <- glm(y~imp,           data=d01, family=poisson)
mod01.08 <- glm(y~forest,        data=d01, family=poisson)
mod01.09 <- glm(y~dev,           data=d01, family=poisson)
mod01.10 <- glm(y~open,          data=d01, family=poisson)
mod01.11 <- glm(y~popden,        data=d01, family=poisson)
mod01.12 <- glm(y~income,        data=d01, family=poisson)
mod01.13 <- glm(y~povrate,       data=d01, family=poisson)


# PELE pop density - linear regression
d02 <- data.frame(y=dx$pele, dex)

mod02.01 <- lm(y~areaha,        data=d02)
mod02.02 <- lm(y~shape,         data=d02)
mod02.03 <- lm(y~age,           data=d02)
mod02.04 <- lm(y~island,        data=d02)
mod02.05 <- lm(y~pct_imp_perim, data=d02)
mod02.06 <- lm(y~tree,          data=d02)
mod02.07 <- lm(y~imp,           data=d02)
mod02.08 <- lm(y~forest,        data=d02)
mod02.09 <- lm(y~dev,           data=d02)
mod02.10 <- lm(y~open,          data=d02)
mod02.11 <- lm(y~popden,        data=d02)
mod02.12 <- lm(y~income,        data=d02)
mod02.13 <- lm(y~povrate,       data=d02)


# ESG6 pop density - linear regression
d03 <- data.frame(y=kx$k6, dex)

mod03.01 <- lm(y~areaha,        data=d03)
mod03.02 <- lm(y~shape,         data=d03)
mod03.03 <- lm(y~age,           data=d03)
mod03.04 <- lm(y~island,        data=d03)
mod03.05 <- lm(y~pct_imp_perim, data=d03)
mod03.06 <- lm(y~tree,          data=d03)
mod03.07 <- lm(y~imp,           data=d03)
mod03.08 <- lm(y~forest,        data=d03)
mod03.09 <- lm(y~dev,           data=d03)
mod03.10 <- lm(y~open,          data=d03)
mod03.11 <- lm(y~popden,        data=d03)
mod03.12 <- lm(y~income,        data=d03)
mod03.13 <- lm(y~povrate,       data=d03)

# SIHI presence/absence - logistic regression
d04 <- data.frame(y=dx$sihi, dex)
d04$y <- ifelse(d04$y > 0, 1, 0)

mod04.01 <- glm(y~areaha,        data=d04, family=binomial)
mod04.02 <- glm(y~shape,         data=d04, family=binomial)
mod04.03 <- glm(y~age,           data=d04, family=binomial)
mod04.04 <- glm(y~island,        data=d04, family=binomial)
mod04.05 <- glm(y~pct_imp_perim, data=d04, family=binomial)
mod04.06 <- glm(y~tree,          data=d04, family=binomial)
mod04.07 <- glm(y~imp,           data=d04, family=binomial)
mod04.08 <- glm(y~forest,        data=d04, family=binomial)
mod04.09 <- glm(y~dev,           data=d04, family=binomial)
mod04.10 <- glm(y~open,          data=d04, family=binomial)
mod04.11 <- glm(y~popden,        data=d04, family=binomial)
mod04.12 <- glm(y~income,        data=d04, family=binomial)
mod04.13 <- glm(y~povrate,       data=d04, family=binomial)


# ESG2 pop density - linear regression
d05 <- data.frame(y=kx$k2, dex)
mod05.01 <- lm(y~areaha,        data=d05)
mod05.02 <- lm(y~shape,         data=d05)
mod05.03 <- lm(y~age,           data=d05)
mod05.04 <- lm(y~island,        data=d05)
mod05.05 <- lm(y~pct_imp_perim, data=d05)
mod05.06 <- lm(y~tree,          data=d05)
mod05.07 <- lm(y~imp,           data=d05)
mod05.08 <- lm(y~forest,        data=d05)
mod05.09 <- lm(y~dev,           data=d05)
mod05.10 <- lm(y~open,          data=d05)
mod05.11 <- lm(y~popden,        data=d05)
mod05.12 <- lm(y~income,        data=d05)
mod05.13 <- lm(y~povrate,       data=d05)


# TAST presence/absence - logistic regression
d06 <- data.frame(y=dx$tast, dex)
d06$y <- ifelse(d06$y > 0, 1, 0)

mod06.01 <- glm(y~areaha,        data=d06, family=binomial)
mod06.02 <- glm(y~shape,         data=d06, family=binomial)
mod06.03 <- glm(y~age,           data=d06, family=binomial)
mod06.04 <- glm(y~island,        data=d06, family=binomial)
mod06.05 <- glm(y~pct_imp_perim, data=d06, family=binomial)
mod06.06 <- glm(y~tree,          data=d06, family=binomial)
mod06.07 <- glm(y~imp,           data=d06, family=binomial)
mod06.08 <- glm(y~forest,        data=d06, family=binomial)
mod06.09 <- glm(y~dev,           data=d06, family=binomial)
mod06.10 <- glm(y~open,          data=d06, family=binomial)
mod06.11 <- glm(y~popden,        data=d06, family=binomial)
mod06.12 <- glm(y~income,        data=d06, family=binomial)
mod06.13 <- glm(y~povrate,       data=d06, family=binomial)
