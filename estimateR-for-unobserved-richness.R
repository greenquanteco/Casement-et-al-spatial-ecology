dx <- read.csv("dat-cap-final.csv")
spp <- c("blca", "glvo", "mipe", "mipi",
         "mumu", "ocnu", "pego", "pele",
         "rano", "rehu", "scca", "sihi", "tast")
dx$species[which(dx$species == "ocmv")] <- "ocnu"
dx$species[which(dx$species == "blbr")] <- "blca"

sites <- sort(unique(dx$site))

cap <- matrix(0, nrow=length(sites), ncol=length(spp))
colnames(cap) <- spp
rownames(cap) <- sites

cap <- as.matrix(ftable(species~site, data=dx))
cap <- cap[,spp]

library(vegan)
estimateR(cap)
