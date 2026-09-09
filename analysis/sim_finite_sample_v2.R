# sim_finite_sample_v2.R -- Paper 4 finite-sample study, version 2.
# Same design and estimator as sim_finite_sample.R (family F of Theorem 2
# plus a binary instrument I; saturated nonparametric AIPW-ATT for four
# adjustment sets), with three additions requested in the revision:
#   (1) the exact ATT bound V_att(S) of each set, computed by enumerating
#       the eight (z, w, i) cells of the discrete design (no simulation),
#       and the exact bound ratio relative to {Z};
#   (2) Monte Carlo standard errors for bias, SD, the variance ratio
#       (paired bootstrap over replications) and coverage;
#   (3) a second seed and R = 2000 as a check on the R = 1000 run.
# Output: output/finite_sample_v2_table.csv, output/finite_sample_v2_summary.txt
args <- commandArgs(trailingOnly = TRUE)
R    <- if (length(args) >= 1) as.integer(args[1]) else 1000
SEED <- if (length(args) >= 2) as.integer(args[2]) else 2026
set.seed(SEED)

aipw_att <- function(y, a, cell) {
  eh <- ave(a, cell); m0 <- ave(ifelse(a == 0, y, NA), cell, FUN = function(x) mean(x, na.rm = TRUE))
  ok <- eh > 0 & eh < 1 & !is.na(m0)
  y <- y[ok]; a <- a[ok]; eh <- eh[ok]; m0 <- m0[ok]
  p <- mean(a); lam <- eh / (1 - eh)
  phi0 <- (a * (y - m0) - lam * (1 - a) * (y - m0)) / p
  psi <- mean(phi0); phi <- phi0 - a * psi / p
  c(est = psi, se = sd(phi) / sqrt(length(y)))
}

# exact ATT bound for adjustment set S in the discrete design ------------------
exact_bounds <- function(e0, e1, delta = 1.5, gamma = 1, theta = 1, shift = 0.8) {
  cells <- expand.grid(z = 0:1, w = 0:1, i = 0:1)
  cells$pr <- 1 / 8
  cells$e  <- plogis(qlogis(ifelse(cells$z == 1, e1, e0)) + shift * (cells$i - .5))
  cells$m1 <- gamma * cells$z + theta                 # E(Y | A = 1, cell)
  cells$m0 <- gamma * cells$z + delta * cells$w       # E(Y | A = 0, cell)
  cells$v1 <- 1; cells$v0 <- 1
  p   <- sum(cells$pr * cells$e)
  psi <- sum(cells$pr * cells$e * (cells$m1 - cells$m0)) / p
  sets <- list(Z = "z", ZW = c("z", "w"), ZI = c("z", "i"), ZWI = c("z", "w", "i"))
  V <- sapply(sets, function(S) {
    key <- do.call(paste, cells[, S, drop = FALSE])
    tot <- 0
    for (k in unique(key)) {
      cc <- cells[key == k, ]; wgt <- cc$pr / sum(cc$pr)
      eS <- sum(wgt * cc$e)
      w1 <- wgt * cc$e / eS; w0 <- wgt * (1 - cc$e) / (1 - eS)
      mu1 <- sum(w1 * cc$m1); mu0 <- sum(w0 * cc$m0)
      s1 <- sum(w1 * (cc$v1 + (cc$m1 - mu1)^2)); s0 <- sum(w0 * (cc$v0 + (cc$m0 - mu0)^2))
      tau <- mu1 - mu0
      tot <- tot + sum(cc$pr) * (eS * s1 + eS^2 / (1 - eS) * s0 + eS * (tau - psi)^2)
    }
    tot / p^2
  })
  list(V = V, ratio = V / V["Z"], p = p, psi = psi)
}

run <- function(n, e0, e1, delta = 1.5, gamma = 1, theta = 1, tag) {
  psi_true <- theta - delta * 0.5
  sets <- list(Z = "z", ZW = c("z", "w"), ZI = c("z", "i"), ZWI = c("z", "w", "i"))
  res <- array(NA, c(R, length(sets), 2), dimnames = list(NULL, names(sets), c("est", "se")))
  for (r in 1:R) {
    z <- rbinom(n, 1, .5); w <- rbinom(n, 1, .5); i <- rbinom(n, 1, .5)
    e <- plogis(qlogis(ifelse(z == 1, e1, e0)) + 0.8 * (i - .5))
    a <- rbinom(n, 1, e)
    y <- gamma * z + delta * (1 - a) * w + theta * a + rnorm(n)
    d <- data.frame(z, w, i)
    for (s in names(sets)) res[r, s, ] <- aipw_att(y, a, do.call(paste, d[, sets[[s]], drop = FALSE]))
  }
  ex <- exact_bounds(e0, e1, delta, gamma, theta)
  # paired bootstrap over replications for the variance-ratio MCSE
  # bootstrap indices from a separate stream so that the simulation draws match sim_finite_sample.R exactly
  saved <- .Random.seed; set.seed(SEED + 1000 * n + 7)
  B <- 500; idx <- matrix(sample.int(R, R * B, replace = TRUE), R, B)
  .Random.seed <<- saved
  out <- do.call(rbind, lapply(names(sets), function(s) {
    est <- res[, s, "est"]; se <- res[, s, "se"]; estZ <- res[, "Z", "est"]
    cov <- mean(abs(est - psi_true) <= 1.96 * se)
    vr  <- var(est) / var(estZ)
    vr_b <- apply(idx, 2, function(ii) var(est[ii]) / var(estZ[ii]))
    data.frame(design = tag, n = n, set = s, R = R, seed = SEED,
               bias = mean(est) - psi_true, bias_mcse = sd(est) / sqrt(R),
               sd = sd(est), sd_mcse = sd(est) / sqrt(2 * (R - 1)),
               rmse = sqrt(mean((est - psi_true)^2)),
               se_ratio = mean(se) / sd(est),
               cov95 = cov, cov_mcse = sqrt(cov * (1 - cov) / R),
               var_rel_Z = vr, var_rel_Z_mcse = if (s == "Z") 0 else sd(vr_b),
               n_var_over_V = n * var(est) / ex$V[[s]],
               V_exact = ex$V[[s]], bound_ratio_exact = ex$ratio[[s]])
  }))
  out
}
t0 <- Sys.time()
tabs <- rbind(
  run( 500, .168, .401, tag = "A: rare treatment (e<1/2)"),
  run(2000, .168, .401, tag = "A: rare treatment (e<1/2)"),
  run( 500, .599, .832, tag = "B: common treatment (e>1/2)"),
  run(2000, .599, .832, tag = "B: common treatment (e>1/2)"))
el <- round(difftime(Sys.time(), t0, units = "mins"), 1)
dir.create("output", showWarnings = FALSE)
suffix <- if (R == 1000 && SEED == 2026) "" else sprintf("_R%d_seed%d", R, SEED)
write.csv(tabs, sprintf("output/finite_sample_v2_table%s.csv", suffix), row.names = FALSE)
sink(sprintf("output/finite_sample_v2_summary%s.txt", suffix))
cat(sprintf("sim_finite_sample_v2.R  R = %d, seed = %d, %s min\n", R, SEED, el))
cat("Exact bounds (enumeration of the eight (z,w,i) cells; Hahn bound with unknown propensity):\n")
for (tg in unique(tabs$design)) {
  x <- tabs[tabs$design == tg & tabs$n == 500, ]
  cat(sprintf("  %s: V = %s; ratio to {Z} = %s\n", tg,
              paste(sprintf("%.4f", x$V_exact), collapse = ", "),
              paste(sprintf("%.3f", x$bound_ratio_exact), collapse = ", ")))
}
cat("\n")
print(format(tabs[, c("design", "n", "set", "bias", "bias_mcse", "sd", "se_ratio", "cov95", "cov_mcse",
                      "var_rel_Z", "var_rel_Z_mcse", "bound_ratio_exact", "n_var_over_V")], digits = 3),
      row.names = FALSE)
sink()
cat(readLines(sprintf("output/finite_sample_v2_summary%s.txt", suffix)), sep = "\n")
