# counterexample_search.R — Paper 4: random-structure scan (Appendix B; the
# 49-structure, 10-divergence result of Section 8) preceded by four Monte
# Carlo design checks T1-T4 that are not reported in the paper.
# Setup: binary pre-treatment covariates Z (SEM), binary A (logistic in
# parents), continuous Y (linear + optional effect modification).
# For each valid adjustment set S (adjustment criterion relative to the
# drawn DAG, enumerated with dagmv; all covariates declared as vertices)
# the semiparametric variance bounds
#   V_ATE(S) = E[ s1/e + s0/(1-e) + (tau - psiATE)^2 ]
#   V_ATT(S) = (1/p^2) E[ e*s1 + e^2/(1-e)*s0 + e*(tau - psiATT)^2 ]
#     (Hahn 1998 bound, propensity unknown; all quantities given S-cells)
# are computed EXACTLY by enumerating the 2^J covariate cells from the
# structural parameters (exact_bounds_for_S); Monte Carlo cell averages on
# a simulated population (bounds_for_S) are kept only for comparison.
# Exactly tied minimizers (relative tolerance 1e-9) are all listed.
# Seeds: 20260818 (design checks T1-T4), 99 (the 60-structure scan).
# Output: output/counterexample_results.txt. Runtime about 40 s.
# Superseded Monte Carlo-only version: counterexample_search_v1_mc.R.
suppressMessages(library(dagmv))
set.seed(20260818); N <- 2e6
dir.create("output", showWarnings = FALSE)

bounds_for_S <- function(Z, A, Y, S) {
  if (length(S) == 0) id <- rep(1L, length(A))
  else id <- as.integer(as.matrix(Z[, S, drop=FALSE]) %*% 2^(seq_along(S)-1)) + 1L
  n <- length(A)
  cells <- sort(unique(id)); K <- length(cells)
  ncc <- tabulate(id, nbins = max(id))[cells]
  n1c <- rowsum(A, id)[, 1]                       # named by sorted cells
  ec  <- pmin(pmax(n1c / ncc, 1e-4), 1 - 1e-4)
  idx1 <- id[A == 1]; idx0 <- id[A == 0]
  m1 <- tapply(Y[A == 1], idx1, mean); v1 <- tapply(Y[A == 1], idx1, var)
  m0 <- tapply(Y[A == 0], idx0, mean); v0 <- tapply(Y[A == 0], idx0, var)
  key <- as.character(cells)
  M1 <- m1[key]; M0 <- m0[key]; V1 <- v1[key]; V0 <- v0[key]
  # cells missing an arm entirely: fall back to overall arm moments
  M1[is.na(M1)] <- mean(Y[A == 1]); M0[is.na(M0)] <- mean(Y[A == 0])
  V1[is.na(V1)] <- 0; V0[is.na(V0)] <- 0
  tau <- M1 - M0; w <- ncc / n; p <- mean(A)
  psiATE <- sum(w * tau)
  psiATT <- sum(w * ec * tau) / sum(w * ec)
  psiATO <- sum(w * ec * (1 - ec) * tau) / sum(w * ec * (1 - ec))
  vate <- sum(w * (V1 / ec + V0 / (1 - ec) + (tau - psiATE)^2))
  vatt <- sum(w * (ec * V1 + ec^2 / (1 - ec) * V0 + ec * (tau - psiATT)^2)) / p^2
  c(V_ATE = vate, V_ATT = vatt, psiATE = psiATE, psiATT = psiATT,
    psiATO = psiATO)
}

lab <- function(S) if (length(S)) paste(S, collapse="+") else "(none)"

# exact bounds from the structural parameters: enumerate the 2^J covariate cells
# (no Monte Carlo); sigma_a^2(z) = 1 in both arms because the outcome error is N(0,1)
exact_bounds_for_S <- function(cells, S) {
  key <- if (length(S)) do.call(paste, cells[, S, drop = FALSE]) else rep("1", nrow(cells))
  p <- sum(cells$pr * cells$e)
  psiATE <- sum(cells$pr * (cells$m1 - cells$m0))
  psiATT <- sum(cells$pr * cells$e * (cells$m1 - cells$m0)) / p
  vate <- vatt <- 0
  for (k in unique(key)) {
    cc <- cells[key == k, ]; pg <- sum(cc$pr); wgt <- cc$pr / pg
    eS <- sum(wgt * cc$e)
    w1 <- wgt * cc$e / eS; w0 <- wgt * (1 - cc$e) / (1 - eS)
    mu1 <- sum(w1 * cc$m1); mu0 <- sum(w0 * cc$m0)
    s1 <- sum(w1 * (1 + cc$m1^2)) - mu1^2; s0 <- sum(w0 * (1 + cc$m0^2)) - mu0^2
    tau <- mu1 - mu0
    vate <- vate + pg * (s1 / eS + s0 / (1 - eS) + (tau - psiATE)^2)
    vatt <- vatt + pg * (eS * s1 + eS^2 / (1 - eS) * s0 + eS * (tau - psiATT)^2)
  }
  c(V_ATE = vate, V_ATT = vatt / p^2, p = p)
}


report <- function(title, Z, A, Y, sets) {
  cat("\n====", title, "====\n")
  out <- t(vapply(sets, function(S) bounds_for_S(Z, A, Y, S), numeric(5)))
  rownames(out) <- vapply(sets, lab, "")
  print(round(out, 4))
  invisible(out)
}

sink("output/counterexample_results.txt")

## ---- Design T1: control-arm-only predictor, rare treatment ------------
## Z1 confounder; P predicts Y(0) only (effect modification); e ~ 0.2.
## Claim to test: V_ATT({Z1,P}) > V_ATT({Z1}) while V_ATE improves.
Z1 <- rbinom(N, 1, .5); P <- rbinom(N, 1, .5)
A  <- rbinom(N, 1, plogis(-1.6 + 1.2*Z1))            # e in ~(.17,.40)
Y  <- 1 + 1.0*Z1 + (1-A)*1.5*P + 1.0*A + rnorm(N, sd=1)
t1 <- report("T1: P predicts Y(0) only, treatment rare (e<1/2 a.s.)",
             data.frame(Z1=Z1, P=P), A, Y,
             list("Z1", c("Z1","P")))

## ---- Design T1b: same but treatment common (e>1/2 a.s.) ---------------
A  <- rbinom(N, 1, plogis( 1.6 - 1.2*Z1))            # e in ~(.60,.83)
Y  <- 1 + 1.0*Z1 + (1-A)*1.5*P + 1.0*A + rnorm(N, sd=1)
t1b <- report("T1b: P predicts Y(0) only, treatment common (e>1/2 a.s.)",
              data.frame(Z1=Z1, P=P), A, Y,
              list("Z1", c("Z1","P")))

## ---- Design T2: treated-arm-only predictor ----------------------------
## Claim: V_ATT identical with/without P; V_ATE improves with P.
A  <- rbinom(N, 1, plogis(-0.4 + 1.2*Z1))
Y  <- 1 + 1.0*Z1 + A*1.5*P + 1.0*A + rnorm(N, sd=1)
t2 <- report("T2: P predicts Y(1) only",
             data.frame(Z1=Z1, P=P), A, Y,
             list("Z1", c("Z1","P")))

## ---- Design T3: instrument deletion (both estimands should improve) ---
W  <- rbinom(N, 1, .5)
A  <- rbinom(N, 1, plogis(-0.8 + 1.0*Z1 + 1.4*W))
Y  <- 1 + 1.2*Z1 + 1.0*A + rnorm(N, sd=1)
t3 <- report("T3: W pure instrument (deletion should help both)",
             data.frame(Z1=Z1, W=W), A, Y,
             list("Z1", c("Z1","W")))

## ---- Design T4: O_ATE vs ATT-argmin divergence with both arms mixed ---
## P0 predicts control arm strongly + treated arm weakly; rare treatment.
P0 <- rbinom(N, 1, .5)
A  <- rbinom(N, 1, plogis(-1.7 + 1.1*Z1))
Y  <- 1 + 1.0*Z1 + 1.6*P0 - A*1.35*P0 + 1.0*A + rnorm(N, sd=1)
t4 <- report("T4: P0 in Pa(Y) both arms (O_ATE includes it); rare treatment",
             data.frame(Z1=Z1, P0=P0), A, Y,
             list("Z1", c("Z1","P0")))

## ---- ATO estimand drift across valid sets -----------------------------
cat("\n==== ATO estimand drift (psiATO by adjustment set) ====\n")
cat(sprintf("T1  : psiATO(Z1) = %.4f  vs psiATO(Z1,P) = %.4f\n",
            t1[1,"psiATO"],  t1[2,"psiATO"]))
cat(sprintf("T4  : psiATO(Z1) = %.4f  vs psiATO(Z1,P0) = %.4f\n",
            t4[1,"psiATO"],  t4[2,"psiATO"]))

## ---- Random-graph scan -------------------------------------------------
## 60 random structures, J=4 binary covariates; count divergences of
## argmin_S V_ATT vs argmin_S V_ATE over valid sets (dagmv validity).
cat("\n==== Random scan: argmin divergence over valid adjustment sets ====\n")
set.seed(99); NR <- 4e5; n_div <- 0; n_ok <- 0; examples <- list()
n_skip_par <- 0; n_skip_valid <- 0; gaps <- c(); prev <- c(); mc_dev <- c(); n_tie <- 0; tie_spread <- 0; ties <- list()   # all argmin divergences, prevalences and skip reasons are recorded
for (g in 1:60) {
  J <- 4
  # random upstream DAG among Z (lower-triangular), A-parents, Y-parents
  Bz <- matrix(rbinom(J*J, 1, .3), J, J); Bz[upper.tri(Bz, TRUE)] <- 0
  paA <- rbinom(J, 1, .6); paY <- rbinom(J, 1, .6)
  if (sum(paY) == 0 || sum(paA) == 0) { n_skip_par <- n_skip_par + 1; next }
  mod <- rbinom(J, 1, .5) * paY                     # effect modifiers
  cf  <- runif(J, .6, 1.2); cfA <- runif(J, .8, 1.4)*paA
  cfY <- runif(J, .6, 1.4)*paY; cfM <- runif(J, -1.6, 1.2)*mod  # arm-asymmetric modification allowed
  Z <- matrix(0L, NR, J)
  for (j in 1:J) {
    lin <- if (j > 1) Z[, 1:(j-1), drop=FALSE] %*% (Bz[j, 1:(j-1)] * cf[1:(j-1)]) else 0
    Z[, j] <- rbinom(NR, 1, plogis(-.3 + lin))
  }
  A <- rbinom(NR, 1, plogis(-1.9 + Z %*% cfA))  # rarer treatment
  Y <- as.vector(1 + Z %*% cfY + A * (1 + Z %*% cfM)) + rnorm(NR)
  # DAG string for dagmv validity
  zn <- paste0("Z", 1:J)
  es <- c()
  for (j in 1:J) for (i in seq_len(j-1)) if (Bz[j,i]) es <- c(es, paste(zn[i], "->", zn[j]))
  es <- c(es, paste(zn[paA==1], "-> A"), paste(zn[paY==1 | mod==1], "-> Y"), "A -> Y")
  # declare every simulated covariate as a vertex (an isolated covariate is a valid addition to any
  # valid set; without the declaration such sets were silently rejected in the 2026-09-08 run)
  gdag <- dag_parse(paste("dag {", paste(c(zn, "A", "Y", es), collapse=" ; "), "}"))
  subs <- lapply(0:(2^J-1), function(i) zn[bitwAnd(i, 2^(0:(J-1))) > 0])
  valid <- Filter(function(S) adjustment_valid(gdag, S, "A", "Y"), subs)
  if (length(valid) < 2) { n_skip_valid <- n_skip_valid + 1; next }
  Zdf <- as.data.frame(Z); names(Zdf) <- zn
  bb_mc <- t(vapply(valid, function(S) bounds_for_S(Zdf, A, Y, S), numeric(5)))   # Monte Carlo cell averages (kept for comparison)
  # exact cell table for this structure
  cells <- expand.grid(rep(list(0:1), J)); names(cells) <- zn
  pr <- rep(1, nrow(cells))
  for (j in 1:J) {
    lin <- if (j > 1) as.matrix(cells[, 1:(j-1), drop=FALSE]) %*% (Bz[j, 1:(j-1)] * cf[1:(j-1)]) else 0
    pj <- plogis(-.3 + lin); pr <- pr * ifelse(cells[, j] == 1, pj, 1 - pj)
  }
  cells$pr <- as.vector(pr)
  cells$e  <- as.vector(plogis(-1.9 + as.matrix(cells[, zn]) %*% cfA))
  cells$m0 <- as.vector(1 + as.matrix(cells[, zn]) %*% cfY)
  cells$m1 <- as.vector(cells$m0 + 1 + as.matrix(cells[, zn]) %*% cfM)
  bb <- t(vapply(valid, function(S) exact_bounds_for_S(cells, S), numeric(3)))
  prev <- c(prev, bb[1, "p"])
  aATE <- which.min(bb[,"V_ATE"]); aATT <- which.min(bb[,"V_ATT"])
  # exact ties: every valid set whose bound is within 1e-9 (relative) of the minimum
  tieATT <- which(bb[,"V_ATT"] <= bb[aATT,"V_ATT"] * (1 + 1e-9))
  tieATE <- which(bb[,"V_ATE"] <= bb[aATE,"V_ATE"] * (1 + 1e-9))
  if (length(tieATT) > 1 || length(tieATE) > 1) {
    n_tie <- n_tie + 1
    ties[[n_tie]] <- list(graph = paste(es, collapse = " ; "),
                          ate_min = paste(vapply(valid[tieATE], lab, ""), collapse = " = "),
                          att_min = paste(vapply(valid[tieATT], lab, ""), collapse = " = "))
  }
  # do tied minimizers agree on the other bound? (if so, the penalty below does not depend on which tied set is named)
  tie_spread <- max(tie_spread, diff(range(bb[tieATE, "V_ATT"])) / bb[aATT, "V_ATT"],
                    diff(range(bb[tieATT, "V_ATE"])) / bb[aATE, "V_ATE"])
  n_ok <- n_ok + 1
  # exact ATT penalty of the ATE-optimal set (0 when the minimizers coincide or tie)
  gapATT <- (bb[aATE,"V_ATT"] - bb[aATT,"V_ATT"]) / bb[aATT,"V_ATT"]
  gapMC  <- (bb_mc[aATE,"V_ATT"] - bb_mc[aATT,"V_ATT"]) / bb_mc[aATT,"V_ATT"]
  mc_dev <- c(mc_dev, gapMC - gapATT)
  if (gapATT > 1e-9) {
    gaps <- c(gaps, gapATT)
    if (gapATT > 0.01) {
      n_div <- n_div + 1
      if (length(examples) < 3)
        examples[[length(examples)+1]] <- list(
          graph = paste(es, collapse=" ; "),
          ate_min = paste(vapply(valid[tieATE], lab, ""), collapse = " = "),
          att_min = paste(vapply(valid[tieATT], lab, ""), collapse = " = "),
          att_penalty_pct = round(100*gapATT, 2))
    }
  }
}
cat("bounds: exact enumeration of the 2^J covariate cells (Monte Carlo cell averages on 4e5 draws kept only for comparison)\n")
cat(sprintf("skipped: no treatment or outcome parent = %d; fewer than two valid sets = %d\n", n_skip_par, n_skip_valid))
cat("exact ATT penalty (pct) of the ATE-optimal set in every structure whose argmins differ (penalty > 0):",
    paste(round(100*sort(gaps), 2), collapse = ", "), "\n")
cat(sprintf("median exact penalty among structures with penalty > 1%%: %.2f%%\n", 100*median(gaps[gaps > 0.01])))
cat(sprintf("structures evaluated: %d; argmin(ATE) != argmin(ATT) with >1%% ATT penalty: %d (%.0f%%)\n",
            n_ok, n_div, 100*n_div/max(n_ok,1)))
cat(sprintf("treatment prevalence P(A=1) across evaluated structures: %.3f to %.3f (%d above 0.5)\n",
            min(prev), max(prev), sum(prev > 0.5)))
cat(sprintf("largest |Monte Carlo - exact| penalty difference: %.3f percentage points\n", 100*max(abs(mc_dev))))
cat(sprintf("structures with exactly tied minimizers (relative tolerance 1e-9; every tied set is listed at the end of this file): %d\n", n_tie))
cat(sprintf("largest relative spread of the other bound among tied minimizers: %.2e (0 means the penalty is the same whichever tied set is named)\n", tie_spread))
for (ex in examples) {
  cat("\n-- example --\n graph:", ex$graph,
      "\n  ATE-optimal:", ex$ate_min, " | ATT-optimal:", ex$att_min,
      " | ATT penalty of using ATE-optimal:", ex$att_penalty_pct, "%\n")
}
cat("\n==== All structures with exactly tied minimizers (relative tolerance 1e-9) ====\n")
for (k in seq_along(ties)) {
  cat(sprintf("[%d] graph: %s\n     ATE-minimizers: %s\n     ATT-minimizers: %s\n",
              k, ties[[k]]$graph, ties[[k]]$ate_min, ties[[k]]$att_min))
}
sink()
cat(readLines("output/counterexample_results.txt"), sep="\n")
