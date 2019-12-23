# Implementation of an imputation model using a linear model with one covariate
#
# The imputation model is a proper CAR model on several covariates

library(INLA)
library(MIINLA)

library(spdep)
library(rgdal)
library(sp)
library(INLA)
library(RColorBrewer)

#Load data
nc.sids <- readOGR(system.file("shapes/sids.shp", package="spData")[1])
proj4string(nc.sids) <- CRS("+proj=longlat +ellps=clrk66")

#Compute covariate and expected counts
nc.sids$NWPROP74 <- (nc.sids$NWBIR74 / nc.sids$BIR74)
# Convert to continuous
nc.sids$NWPROP74 <- log(nc.sids$NWPROP74 / (1- nc.sids$NWPROP74))
# Scale
nc.sids$NWPROP74 <- scale(nc.sids$NWPROP74)[, 1]
# Expected counts
nc.sids$EXP74 <- nc.sids$BIR74 * sum(nc.sids$SID74) / sum(nc.sids$BIR74)

#SMR74
nc.sids$SMR74 <- nc.sids$SID74 / nc.sids$EXP74

#Create indices of missing observations under MCAR and MNAR
set.seed(1)
idx.mis <- list(MCAR = sample(1:100, 50))

# MNAR
mis.prob <- log(0.5 / (1 - 0.5)) + 5 * nc.sids$NWPROP74
mis.prob <- 1 / (1 + exp(-mis.prob))
idx.mis$MNAR <- sample(1:100, 50, prob = mis.prob)


# Display SMR and covariate
pdf(file = "NCSIDS-SMR.pdf", width = 10, height = 5)
print(spplot(nc.sids, "SMR74",
  col.regions = colorRampPalette(brewer.pal (9, "Oranges"))(16)))
dev.off()

pdf(file = "NCSIDS-NWPROP74.pdf",, width = 10, height = 5)
print(spplot(nc.sids, "NWPROP74",
  col.regions = colorRampPalette(brewer.pal (9, "Blues"))(16)))
dev.off()


#Missing covariates
nc.sids$MCAR10 <- nc.sids$NWPROP74
nc.sids$MCAR30 <- nc.sids$NWPROP74
nc.sids$MCAR50 <- nc.sids$NWPROP74
nc.sids$MNAR10 <- nc.sids$NWPROP74
nc.sids$MNAR30 <- nc.sids$NWPROP74
nc.sids$MNAR50 <- nc.sids$NWPROP74

nc.sids$MCAR10[idx.mis$MCAR[1:10]] <- NA
nc.sids$MCAR30[idx.mis$MCAR[1:30]] <- NA
nc.sids$MCAR50[idx.mis$MCAR[1:50]] <- NA
nc.sids$MNAR10[idx.mis$MNAR[1:10]] <- NA
nc.sids$MNAR30[idx.mis$MNAR[1:30]] <- NA
nc.sids$MNAR50[idx.mis$MNAR[1:50]] <- NA


#Display missing observations in the covariates
pdf(file = "NCSIDS-NWPmiss.pdf", width = 10, height = 6)
print(spplot(nc.sids, c("MCAR10", "MNAR10", "MCAR30", "MNAR30", 
  "MCAR50", "MNAR50"),
  col.regions = colorRampPalette(brewer.pal (9, "Blues"))(16),
sp.layout = list(
        list("sp.polygons", nc.sids, first = TRUE, fill = "grey")
)))
dev.off()


# Standard model (fit t FULL DATASET)
m0 <- inla(SID74 ~ NWPROP74, family = "poisson", data = as.data.frame(nc.sids),
 E = EXP74)
summary(m0)


# Adjacency structure
adj <- poly2nb(nc.sids)
W <- as(nb2mat(adj, style = "B"), "Matrix")
W.scaled <- W / max(eigen(W)$values)

#Index for fitting models
nc.sids$idx <- 1:nrow(nc.sids)




save(file = "01_data.RData", list = ls())
