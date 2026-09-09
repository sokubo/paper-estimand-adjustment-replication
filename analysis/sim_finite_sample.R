# sim_finite_sample.R — Paper 4 finite-sample study.
# Family F of Theorem 2 (+ an instrument I): saturated nonparametric
# AIPW-ATT for four adjustment sets. Reports bias, SD, RMSE, EIF-SE
# calibration, coverage, and the bound ratio predicted by Lemma 1.
set.seed(2026); R <- 1000
aipw_att <- function(y, a, cell) {
  # saturated plug-ins: cell-wise e and mu0; EIF-based SE
  eh <- ave(a, cell); m0 <- ave(ifelse(a==0, y, NA), cell, FUN=function(x) mean(x, na.rm=TRUE))
  ok <- eh > 0 & eh < 1 & !is.na(m0)
  y <- y[ok]; a <- a[ok]; eh <- eh[ok]; m0 <- m0[ok]
  p <- mean(a); lam <- eh/(1-eh)
  phi0 <- (a*(y-m0) - lam*(1-a)*(y-m0))/p
  psi <- mean(phi0); phi <- phi0 - a*psi/p
  c(est = psi, se = sd(phi)/sqrt(length(y)))
}
h <- function(x) x^2/(1-x)
run <- function(n, e0, e1, delta = 1.5, gamma = 1, theta = 1, tag) {
  # Z ~ B(.5) confounder; W ~ B(.5) control-arm-only predictor; I ~ B(.5) instrument
  psi_true <- theta - delta*0.5
  sets <- list(Z="z", ZW=c("z","w"), ZI=c("z","i"), ZWI=c("z","w","i"))
  res <- array(NA, c(R, length(sets), 2), dimnames=list(NULL, names(sets), c("est","se")))
  for (r in 1:R) {
    z <- rbinom(n,1,.5); w <- rbinom(n,1,.5); i <- rbinom(n,1,.5)
    e <- plogis(qlogis(ifelse(z==1,e1,e0)) + 0.8*(i-.5))    # instrument shifts e within z
    a <- rbinom(n,1,e)
    y <- gamma*z + delta*(1-a)*w + theta*a + rnorm(n)
    d <- data.frame(z,w,i)
    for (s in names(sets)) {
      cell <- do.call(paste, d[, sets[[s]], drop=FALSE])
      res[r,s,] <- aipw_att(y, a, cell)
    }
  }
  out <- do.call(rbind, lapply(names(sets), function(s) {
    est <- res[,s,"est"]; se <- res[,s,"se"]
    data.frame(design=tag, n=n, set=s, bias=mean(est)-psi_true, sd=sd(est),
               rmse=sqrt(mean((est-psi_true)^2)), se_ratio=mean(se)/sd(est),
               cov95=mean(abs(est-psi_true) <= 1.96*se))
  }))
  out$var_rel_Z <- out$sd^2/out$sd[out$set=="Z"]^2
  out
}
t0 <- Sys.time()
tabs <- rbind(
  run( 500, .168, .401, tag="A: rare treatment (e<1/2)"),
  run(2000, .168, .401, tag="A: rare treatment (e<1/2)"),
  run( 500, .599, .832, tag="B: common treatment (e>1/2)"),
  run(2000, .599, .832, tag="B: common treatment (e>1/2)"))
cat("time:", round(difftime(Sys.time(), t0, units="mins"),1), "min\n")
print(tabs, digits=3, row.names=FALSE)
dir.create("output", showWarnings=FALSE)
write.csv(tabs, "output/finite_sample_table.csv", row.names=FALSE)
capture.output(print(tabs, digits=3, row.names=FALSE), file="output/finite_sample_summary.txt")
