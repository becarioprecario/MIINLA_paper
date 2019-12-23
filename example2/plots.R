# PLots of imputed values

library(INLA)


load("03_MNAR50c.RData")

# Models in the data:
#
# r: Model with imputation under MCAR from 02_MI.R
# mcar.inla: MOdel with imputation and joint model under MCAR
# mnar.inla: Model with imputation and joint model under MNAR
#
# NOte that 'r' and 'mcar.inla' are the same model in practice.


# Set areas with missing values in BOTH MCAR and MNAR to compare
# We do this so that we have a higher variaety of counties
idx.plot <- sort(idx.mis$MCAR[(idx.mis$MCAR %in%  idx.mis$MNAR)])[1:9]

# Just take the first 9 areas with missing observations
#idx.plot <- sort(idx.mis$MNAR)[1:9]

pdf(file = "ncsids-imputed.pdf")
par(mfrow = c(3, 3))

for(i in 1:9) {
  j <- idx.plot[i]

  # Proportion of missing neighbors
  mis.neig <- round(mean(is.na(nc.sids$NWPROP74M[ W[j,]  == 1])), 2)

  plot(r$marginals.random$idxNA[[j]], type = "l", xlab = "", ylab = "",
    xlim = c(-3, 4), ylim = c(0, 1), lty = 2,
    main = paste0(nc.sids$NAME[j], " (", mis.neig, ")")) #paste0("county ", j))
  lines(mnar.inla$marginals.random$idxNA[[j]], type = "l", lty = 3)

  abline(v = nc.sids$NWPROP74[j])
}
dev.off()



# Try with other dataset

load("03_MNAR05c.RData")

idx.plot <- idx.mis$MNAR[1:5]

par(mfrow = c(3, 2))

for(i in 1:5) {
  j <- idx.plot[i]

  # Proportion of missing neighbors
  mis.neig <- round(mean(is.na(nc.sids$NWPROP74M[ W[j,]  == 1])), 2)

  plot(r$marginals.random$idxNA[[j]], type = "l", xlab = "", ylab = "",
    xlim = c(-3, 4), ylim = c(0, 1), lty = 2,
    main = paste0(nc.sids$NAME[j], " (", mis.neig, ")")) #paste0("county ", j))
  lines(mnar.inla$marginals.random$idxNA[[j]], type = "l", lty = 3)

  abline(v = nc.sids$NWPROP74[j])
}

