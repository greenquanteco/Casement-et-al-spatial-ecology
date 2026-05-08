# Compare versions 3 and 4 of datasets:
# version 3 = dat_rodent_popden_2025-01-24
# version 4 = dat_rodent_popden_site_2025-02-25

library(vegan)
library(dunn.test)
library(pROC)

# final explanatory variable matrix @ 100 m
in.name <- "dat_explanatory_0100_2026-04-27.csv"
ex <- read.csv(in.name)

# final response variable matrix - by species
dx3 <- read.csv("dat_rodent_popden_site_2025-02-25.csv")
dx4 <- read.csv("dat_rodent_popden_site_2025-01-24.csv")

# final respons variable matrix - by species
tx <- read.csv("dat_response_traits_2026-04-27.csv")

rownames(dx3) <- dx3$site
rownames(dx4) <- dx4$site
rownames(ex) <- ex$site
rownames(tx) <- tx$site

dx3$site <- NULL
dx4$site <- NULL
ex$site <- NULL
tx$site <- NULL

dx3$rich <- specnumber(dx3)
dx4$rich <- specnumber(dx4)

all(rownames(ex) == rownames(dx3))
all(rownames(ex) == rownames(dx4))
all(rownames(ex) == rownames(tx))

# univariate analyses
## y vs. urbanization
kruskal.test(dx3$rich~ex$type)
kruskal.test(dx4$rich~ex$type)

## shows that dx3 was used to create manuscript!
kruskal.test(dx3$pele~ex$type)
kruskal.test(dx4$pele~ex$type)

## confirms!
kruskal.test(dx3$sihi~ex$type)
kruskal.test(dx4$sihi~ex$type)

kruskal.test(dx3$tast~ex$type)
kruskal.test(dx4$tast~ex$type)
                 
kruskal.test(tx$k6~ex$type)
kruskal.test(tx$k2~ex$type)

dunn.test(dx3$pele, ex$type)
dunn.test(dx3$tast, ex$type)

# conclusion: regenerate analysis with dx3 (matches univariate)