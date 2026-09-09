# counterexample_search.R — Paper 4 theory design: numerical verification
# of estimand-dependence of efficient adjustment sets.
# Setup: binary pre-treatment covariates Z (SEM), binary A (logistic in
# parents), continuous Y (linear + optional effect modification).
# For each valid adjustment set S we compute the semiparametric variance
# bounds by exact cell averages on a large simulated population:
#   V_ATE(S) = E[ s1/e + s0/(1-e) + (tau - psiATE)^2 ]
#   V_ATT(S) = (1/p^2) E[ e*s1 + e^2/(1-e)*s0 + e*(tau - psiATT)^2 ]
#     (Hahn 1998 bound, propensity unknown; all quantities given S-cells)
# plus the ATO estimand psi_h(S) = E[e(1-e)tau]/E[e(1-e)] to measure
# estimand drift across adjustment sets.
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
n_skip_par <- 0; n_skip_valid <- 0; gaps <- c()   # all argmin divergences and skip reasons are recorded
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
  gdag <- dag_parse(paste("dag {", paste(es, collapse=" ; "), "}"))
  subs <- lapply(0:(2^J-1), function(i) zn[bitwAnd(i, 2^(0:(J-1))) > 0])
  valid <- Filter(function(S) adjustment_valid(gdag, S, "A", "Y"), subs)
  if (length(valid) < 2) { n_skip_valid <- n_skip_valid + 1; next }
  Zdf <- as.data.frame(Z); names(Zdf) <- zn
  bb <- t(vapply(valid, function(S) bounds_for_S(Zdf, A, Y, S), numeric(5)))
  aATE <- which.min(bb[,"V_ATE"]); aATT <- which.min(bb[,"V_ATT"])
  n_ok <- n_ok + 1
  if (aATE != aATT) {
    # require the gap to exceed MC noise (relative 1%)
    gapATT <- (bb[aATE,"V_ATT"] - bb[aATT,"V_ATT"]) / bb[aATT,"V_ATT"]
    gaps <- c(gaps, gapATT)
    if (gapATT > 0.01) {
      n_div <- n_div + 1
      if (length(examples) < 3)
        examples[[length(examples)+1]] <- list(
          graph = paste(es, collapse=" ; "),
          ate_min = lab(valid[[aATE]]), att_min = lab(valid[[aATT]]),
          att_penalty_pct = round(100*gapATT, 1))
    }
  }
}
cat(sprintf("skipped: no treatment or outcome parent = %d; fewer than two valid sets = %d\n", n_skip_par, n_skip_valid))
cat("ATT penalty (pct) of the ATE-optimal set in every structure whose argmins differ:",
    paste(round(100*sort(gaps), 1), collapse = ", "), "\n")
cat(sprintf("structures evaluated: %d; argmin(ATE) != argmin(ATT) with >1%% ATT penalty: %d (%.0f%%)\n",
            n_ok, n_div, 100*n_div/max(n_ok,1)))
for (ex in examples) {
  cat("\n-- example --\n graph:", ex$graph,
      "\n  ATE-optimal:", ex$ate_min, " | ATT-optimal:", ex$att_min,
      " | ATT penalty of using ATE-optimal:", ex$att_penalty_pct, "%\n")
}
sink()
cat(readLines("output/counterexample_results.txt"), sep="\n")
