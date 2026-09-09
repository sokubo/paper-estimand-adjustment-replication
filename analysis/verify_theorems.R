# verify_theorems.R — Paper 4: machine-precision verification of the
# algebra behind Theorems 1-2 BEFORE the proofs are written up.
# Engine: exact computation of the S-indexed ATT bound on discrete
# covariate designs given as probability tables (no Monte Carlo except
# Part D). Covariates (Z, W) binary; cells c = (z, w) with pmf pi(c),
# propensity e(c), arm means mu0(c), mu1(c), arm noise v0(c), v1(c).
options(digits = 12)
h <- function(x) x^2 / (1 - x)

bound_exact <- function(tab, S) {
  # tab: data.frame(z, w, pi, e, mu0, mu1, v0, v1); S in {"z","w",c("z","w"),character(0)}
  p   <- sum(tab$pi * tab$e)
  psi <- sum(tab$pi * tab$e * (tab$mu1 - tab$mu0)) / p
  if (length(S) == 0) key <- rep(1, nrow(tab)) else
    key <- do.call(paste, tab[, S, drop = FALSE])
  out <- 0
  for (s in unique(key)) {
    i <- key == s
    pis <- sum(tab$pi[i])
    w_c <- tab$pi[i] / pis                       # P(c | S=s)
    eS  <- sum(w_c * tab$e[i])                   # P(A=1 | S=s)
    w1  <- w_c * tab$e[i] / eS                   # P(c | S=s, A=1)
    w0  <- w_c * (1 - tab$e[i]) / (1 - eS)       # P(c | S=s, A=0)
    m1  <- sum(w1 * tab$mu1[i]); m0 <- sum(w0 * tab$mu0[i])
    s1  <- sum(w1 * (tab$v1[i] + tab$mu1[i]^2)) - m1^2
    s0  <- sum(w0 * (tab$v0[i] + tab$mu0[i]^2)) - m0^2
    out <- out + pis * (eS * s1 + h(eS) * s0 + eS * (m1 - m0 - psi)^2)
  }
  out / p^2
}

cat("== Part A: Theorem 2 family — closed form vs exact engine ==\n")
# Family: Z~Bern(r), W~Bern(q) indep; e = e(Z); Y = g*Z + d*(1-A)*W + t*A + N(0,1)
fam_tab <- function(e0, e1, r, q, g, d, t) {
  tab <- expand.grid(z = 0:1, w = 0:1)
  tab$pi <- ifelse(tab$z == 1, r, 1 - r) * ifelse(tab$w == 1, q, 1 - q)
  tab$e  <- ifelse(tab$z == 1, e1, e0)
  tab$mu0 <- g * tab$z + d * tab$w
  tab$mu1 <- g * tab$z + t
  tab$v0 <- 1; tab$v1 <- 1
  tab
}
closed_forms <- function(e0, e1, r, q, g, d, t) {
  p  <- (1 - r) * e0 + r * e1
  Eh <- (1 - r) * h(e0) + r * h(e1)
  vw <- d^2 * q * (1 - q)
  c(V_S1 = (p + Eh * (1 + vw)) / p^2,          # S1 = {Z}
    V_S2 = (p + Eh + p * vw) / p^2,            # S2 = {Z, W}
    diff_sign = sign((p - Eh)))                # V_S2 - V_S1 sign = sign(E[e(1-2e)/(1-e)])
}
for (cfg in list(c(.168, .401, .5, .5, 1, 1.5, 1),          # main example of Section 5 (e<1/2)
                 c(.60, .80, .5, .5, 1, 1.5, 1),            # T1b-type (e>1/2)
                 c(.25, .45, .3, .6, .7, -1.2, .5))) {      # asymmetric extra
  tab <- fam_tab(cfg[1], cfg[2], cfg[3], cfg[4], cfg[5], cfg[6], cfg[7])
  cf  <- closed_forms(cfg[1], cfg[2], cfg[3], cfg[4], cfg[5], cfg[6], cfg[7])
  ex1 <- bound_exact(tab, "z"); ex2 <- bound_exact(tab, c("z", "w"))
  cat(sprintf(" e=(%.3f,%.3f): closed (%.8f, %.8f) exact (%.8f, %.8f)  max|err|=%.2e  argminATT=%s\n",
      cfg[1], cfg[2], cf[1], cf[2], ex1, ex2,
      max(abs(cf[1] - ex1), abs(cf[2] - ex2)),
      ifelse(ex1 < ex2, "{Z}", "{Z,W}")))
}

cat("\n== Part B: Theorem 1(b) supplementation increment identity ==\n")
# Random designs with W ⊥ A | Z (e depends on z only); arbitrary mu tables.
set.seed(1); errB <- 0
for (rep in 1:300) {
  tab <- expand.grid(z = 0:1, w = 0:1)
  pz <- runif(1, .2, .8); pw <- matrix(runif(2, .2, .8), 2)  # P(w=1|z)
  tab$pi <- ifelse(tab$z == 1, pz, 1 - pz) *
            ifelse(tab$w == 1, pw[tab$z + 1], 1 - pw[tab$z + 1])
  ez <- runif(2, .05, .95)
  tab$e <- ez[tab$z + 1]
  tab$mu0 <- rnorm(4, 0, 2); tab$mu1 <- rnorm(4, 0, 2)
  tab$v0 <- runif(4, .5, 2); tab$v1 <- runif(4, .5, 2)
  lhs <- (sum(tab$pi * tab$e))^2 *
         (bound_exact(tab, c("z", "w")) - bound_exact(tab, "z"))
  # RHS: E{ e(S) [ (1-2e)/(1-e) Var(mu0|S) - 2 Cov(mu1, mu0|S) ] }
  rhs <- 0
  for (z in 0:1) {
    i <- tab$z == z; wz <- tab$pi[i] / sum(tab$pi[i])   # P(w|z)
    m0 <- sum(wz * tab$mu0[i]); m1 <- sum(wz * tab$mu1[i])
    V0 <- sum(wz * tab$mu0[i]^2) - m0^2
    CV <- sum(wz * tab$mu1[i] * tab$mu0[i]) - m1 * m0
    e  <- ez[z + 1]
    rhs <- rhs + sum(tab$pi[i]) * e * ((1 - 2*e)/(1 - e) * V0 - 2 * CV)
  }
  errB <- max(errB, abs(lhs - rhs))
}
cat(sprintf(" 300 random designs: max |LHS - RHS| = %.3e\n", errB))

cat("\n== Part C: Theorem 1(a) deletion increment identity ==\n")
# W ⊥ Y(a) | Z (mu, v free of w) but e depends on (z, w).
set.seed(2); errC <- 0; negC <- 0
for (rep in 1:300) {
  tab <- expand.grid(z = 0:1, w = 0:1)
  pz <- runif(1, .2, .8); pw <- runif(2, .2, .8)
  tab$pi <- ifelse(tab$z == 1, pz, 1 - pz) *
            ifelse(tab$w == 1, pw[tab$z + 1], 1 - pw[tab$z + 1])
  tab$e  <- runif(4, .05, .95)
  mu0z <- rnorm(2, 0, 2); mu1z <- rnorm(2, 0, 2); v0z <- runif(2, .5, 2); v1z <- runif(2, .5, 2)
  tab$mu0 <- mu0z[tab$z + 1]; tab$mu1 <- mu1z[tab$z + 1]
  tab$v0  <- v0z[tab$z + 1];  tab$v1  <- v1z[tab$z + 1]
  Delta <- (sum(tab$pi * tab$e))^2 *
           (bound_exact(tab, c("z", "w")) - bound_exact(tab, "z"))
  # RHS: E_S[ sigma0^2(S) * ( E[h(e_SW)|S] - h(e_S) ) ]
  rhs <- 0
  for (z in 0:1) {
    i <- tab$z == z; wz <- tab$pi[i] / sum(tab$pi[i])
    eS <- sum(wz * tab$e[i])
    # sigma0^2(S=z): W|z,A=0 tilts by (1-e)!
    w0 <- wz * (1 - tab$e[i]) / (1 - eS)
    s0 <- sum(w0 * (tab$v0[i] + tab$mu0[i]^2)) - (sum(w0 * tab$mu0[i]))^2
    rhs <- rhs + sum(tab$pi[i]) * (sum(wz * h(tab$e[i])) - h(eS)) * s0
  }
  errC <- max(errC, abs(Delta - rhs)); if (Delta < -1e-12) negC <- negC + 1
}
cat(sprintf(" 300 random designs: max |LHS - RHS| = %.3e ; Delta<0 count = %d (expect 0)\n",
            errC, negC))

cat("\n== Part D: finite-sample AIPW-ATT variance vs bound (T1 config) ==\n")
set.seed(3); R <- 400; n <- 20000
tab <- fam_tab(.168, .401, .5, .5, 1, 1.5, 1)
for (S in list("z", c("z", "w"))) {
  est <- numeric(R)
  for (r in 1:R) {
    z <- rbinom(n, 1, .5); w <- rbinom(n, 1, .5)
    e <- ifelse(z == 1, .401, .168); a <- rbinom(n, 1, e)
    y <- z + 1.5 * (1 - a) * w + a + rnorm(n)
    cell <- if (length(S) == 1) z else 2 * z + w
    eh <- ave(a, cell); m0 <- ave(y * 0, cell)   # placeholders
    m0 <- ave(ifelse(a == 0, y, NA), cell, FUN = function(x) mean(x, na.rm = TRUE))
    ph <- mean(a); lam <- eh / (1 - eh)
    est[r] <- mean(a * (y - m0) - lam * (1 - a) * (y - m0)) / ph
  }
  Vb <- bound_exact(tab, S)
  cat(sprintf(" S={%s}: n*Var_MC = %.3f vs bound = %.3f (ratio %.3f); mean est = %.4f (psi = %.4f)\n",
      paste(S, collapse = ","), n * var(est), Vb, n * var(est) / Vb,
      mean(est), sum(tab$pi * tab$e * (tab$mu1 - tab$mu0)) / sum(tab$pi * tab$e)))
}
