# Data for Casement et al. manuscript on urban small mammal spatial ecology

# Description

This repo contains the data supporting the manuscript by Bri Casement, Leslie Lopez, Katelyn Nestor, and Nicholas Green entitled, "Urbanization restructures small mammal communities without reducing species richness". Full details of how the data were collected are contained in the manuscript, a draft of which will be posted here once finalized and submitted to the journal. This manuscript is derived from the lead author's [masters thesis](https://digitalcommons.kennesaw.edu/masterstheses/31/).

# Datasets

There are 3 datasets used to run the analyses.

## dat_response_species_2026-04-27.csv

Contains species population densities at each site. Variables:

-   mipi: population density of Pine Vole (*Microtus pinetorum*).

-   mumu: population density of House Mouse (*Mus musculus*).

-   pego: population density of Cotton Mouse (*Peromyscus gossypinus*).

-   pele: population density of White-footed Mouse (*Peromyscus leucopus*).

-   rano: population density of Brown Rat (*Rattus norvegicus*).

-   rehu: population density of Eastern Harvest Mouse (*Reithrodontomys humilis*).

-   sihi: population density of Hispid Cotton Rat (*Sigmodon hispidus*).

-   tast: population density of Eastern Chipmunk (*Tamias striatus*).

-   blca: population density of Southern Short-Tailed Shrew (*Blarina Brevicauda*).

-   glvo: population density of Southern Flying Squirrel (*Glaucomys volans*).

-   ocnu: population density of Golden Mouse (*Ochrotomys nuttalli*).

-   scca: population density of Eastern Gray Squirrel (*Sciurus carolinensis*).

-   rich: small mammal species richness

-   simp: Simpson's diversity index of small mammals.

## dat_response_traits_2026-04-27.csv

Contains individual densities (n/ha) by ecological strategy group (ESGs).

-   k1: Individual density of ESG 1 (diurnal tree squirrels).

-   k2: Individual density of ESG 2 (specialist muroids).

-   k3: Individual density of ESG 3 (large, arboreal, long-lived squirrels and muroids).

-   k4: Individual density of ESG 4 (brown rats).

-   k5: Individual density of ESG 5 (shrews).

-   k6: Individual density of ESG 6 (generalist muroids).

-   rich: richness of ESGs.

-   simp: Simpson's diversity index of ESGs.

## dat_explanatory_2026-04-27.csv

Contains explanatory variables derived from on-site measurements and GIS.

-   areaha: area of site in ha

-   type: type of site (rural, suburban, urban)

-   elev.mean: mean elevation of site in m

-   tempc_day: mean daily temperature during trapping (deg. C)

-   lux_day: mean daily brightness reading during trapping (lux).

-   db.mn: mean daytime sound level (db) during trapping.

-   shape: shape complexity of site, calculated as $0.25P/\sqrt(A)$, where P is site perimeter and A is site area, given that P and $\sqrt(A)$ have the same units.

-   age: Years since site was last connected to nearby site with similar habitat.

-   island: distance to nearest patch of similar habitat, in m.

-   pct.imp: percentage of the site perimeter that is impervious surface.

-   tree: percentage of land surrounding site (within 100 m) covered by tree cover of any type (not used; slightly different than forest (below)).

-   imp: percentage of land in 100 m buffer surrounding site covered by impervious surface.

-   areaha100: area of the 100 m buffer surrounding the site, in ha.

-   popden: human population density, in people / ha, in the 100 m buffer surrounding each site.

-   forest: mean percentage forest of each 30 m pixel in the the 100 m buffer surrounding each site.

-   dev: mean percentage developed land cover of each 30 m pixel in the 100 m buffer surrounding each site.

-   wetland: mean percentage wetland land cover of each 30 m pixel in the 100 m buffer surrounding each site

-   open: mean percentage of open land cover of each 30 m pixel in the 100 m buffer surrounding each site

-   income: median household income within the 100 m buffer surrounding each site

-   povrate: percentage of households in the 100 m buffer surrounding each site with incomes below the federal poverty line.

# Scripts

There are 2 scripts that run the analyses. 1. "analysis-univariate.r": runs analyses of univariate responses such as species richness or *P. leucopus* population density. 2. "analysis-multivariate.r": runs multivariate analyses of community structure and beta diversity.

All analyses were performed using R version 4.5.1.

Script **analysis-univariate.r** used packages vegan 2.6-10, dunn.test version 1.3.6, and pROC version 1.18.5.

Script **analysis-multivariate.r** used packages vegan 2.6-10 and betapart version 1.6.1.
