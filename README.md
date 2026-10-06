# Data for Casement et al. manuscript on urban small mammal spatial ecology

## Description

This repo contains the data supporting the revised manuscript by Bri Casement, Leslie Lopez, Katelyn Nestor, and Nicholas Green entitled, "Urban small mammal communities are different, not depauperate: species-specific responses along an urban-to-rural gradien". This release corresponds to the revised version submitted to Journal of Mammalogy on 2026-10-06 as manuscript ID JMAMM-2026-091. This manuscript is derived from the lead author's [masters thesis](https://digitalcommons.kennesaw.edu/masterstheses/31/).

Please direct all correspondence to [ngreen62@kennesaw.edu](mailto:ngreen62@kennesaw.edu). Information about activities in the Green Quantitative Ecology lab at Kennesaw State University can be found on the [lab website](https://greenquanteco.github.io).

## Data files

There are 3 data files used to run the analyses. Inputs for the script are in the `data` folder.

### dat-mammal.csv


Contains species population densities at each site. Variables:

-   `mipi`: detection (1) or nondetection (0) of Pine Vole (*Microtus pinetorum*).
-   `mumu`: detection (1) or nondetection (0) of House Mouse (*Mus musculus*).
-   `pego`: detection (1) or nondetection (0) of Cotton Mouse (*Peromyscus gossypinus*).
-   `pele`: population density of White-footed Mouse (*Peromyscus leucopus*) in individuals/ha.
-   `rano`: detection (1) or nondetection (0) of Brown Rat (*Rattus norvegicus*).
-   `rehu`: detection (1) or nondetection (0) of Eastern Harvest Mouse (*Reithrodontomys humilis*).
-   `sihi`: detection (1) or nondetection (0) of Hispid Cotton Rat (*Sigmodon hispidus*).

-   `tast`: detection (1) or nondetection (0) of Eastern Chipmunk (*Tamias striatus*).
-   `blca`: detection (1) or nondetection (0) of Southern Short-Tailed Shrew (*Blarina Brevicauda*).
-   `glvo`: detection (1) or nondetection (0) of Southern Flying Squirrel (*Glaucomys volans*).
-   `ocnu`: detection (1) or nondetection (0) of Golden Mouse (*Ochrotomys nuttalli*).
-   `scca`: detection (1) or nondetection (0) of Eastern Gray Squirrel (*Sciurus carolinensis*).

### dat-explanatory.csv

Contains explanatory variables derived from on-site measurements and GIS.

-   `areaha`: area of site in ha
-   `shape`: shape complexity of site, calculated as $0.25P/\sqrt(A)$, where P is site perimeter and A is site area, given that P and $\sqrt(A)$ have the same units.
-   `age`: Years since site was last connected to nearby site with similar habitat.
-   `perim`: percentage of the site perimeter that is impervious surface.
-   `island`: distance to nearest patch of similar habitat, in m.
-   `imp`: percentage of land in 100 m buffer surrounding site covered by impervious surface.
-   `forest`: mean percentage forest of each 30 m pixel in the 100 m buffer surrounding each site.
-   `dev`: mean percentage developed land cover of each 30 m pixel in the 100 m buffer surrounding each site.
-   `tree`: percentage of land surrounding site (within 100 m) covered by tree cover of any type (not used; slightly different than forest (below)).
-   `open`: mean percentage of open land cover of each 30 m pixel in the 100 m buffer surrounding each site
-   `popden`: human population density, in people / ha, in the 100 m buffer surrounding each site.
-   `povrate`: percentage of households in the 100 m buffer surrounding each site with incomes below the federal poverty line.
-   `human`: Human Modification Index wtihin a 100 m buffer surrounding each site. Ranges from 0 (no human modification of landscape) to 1 (maximal modification).

### weights_15.7km.rds

R object containing row-standardized spatial weights based on a 15.7 km neighborhood for testing spatial autocorrelation in model residuals.

## Species silhouettes

Figure 3 utilizes 3 silhouettes downloaded from [PhyloPic](https://www.phylopic.org). These are stored in the `images` folder.

- `pele.png`: White-footed mouse (*Peromyscus leucopus*). Credit: Edwin Price ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/))
- `sihi.png`: Hispid cotton rat (*Sigmodon hispidus*). Public domain.
- `tamias.png`: Eastern chipmunk (*Tamias striatus*). Credit: Chloé Schmidt ([CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)). 

## Scripts

There is 1 script, `analysis-release.r`, that runs all analyses. The script is stored in the project root. Users should open this R script directly in RStudio so that the directory structure in the code is preserved.

## Software notes

All analyses were performed using R version 4.5.1. This analysis used packages `vegan` version 2.6-10, `betapart` version 1.6.1, `png` version 0.1-9, `spdep` version 1.4-1, and `MuMIn` version 1.48.11.

## Outputs

- Figure 2: Correlations between explanatory variable, grouped by hypotheses 1, 2, and 3.
- Figure 3: Small mammal population density and detections vs. the best supported spatial predictor variables.
- Figure 4: Detections of small mammal species at 23 sites along the Atlanta, Georgia, USA urban-to-rural gradient.
