# verify_thm3.R — machine-precision check of Theorem 3 (paper 4):
# under no effect modification, V_att(Z) - V_att(O) equals the two-term
# expression (5) of the paper and is >= 0 for every graphically valid Z.
# Design: Z1 confounder (Z1->A, Z1->Y), W outcome-only parent (W->Y),
# I instrument (I->A). O = {Z1, W}. Valid sets: {Z1},{Z1,W},{Z1,I},{Z1,W,I}.
options(digits = 12); h <- function(x) x^2/(1-x)

bound_exact <- function(tab, S) {
  p <- sum(tab$pi * tab$e); psi <- sum(tab$pi*tab$e*(tab$mu1-tab$mu0))/p
  key <- if (length(S)) do.call(paste, tab[, S, drop=FALSE]) else rep(1, nrow(tab))
  out <- 0
  for (s in unique(key)) {
    i <- key == s; pis <- sum(tab$pi[i]); wc <- tab$pi[i]/pis
    eS <- sum(wc*tab$e[i]); w1 <- wc*tab$e[i]/eS; w0 <- wc*(1-tab$e[i])/(1-eS)
    m1 <- sum(w1*tab$mu1[i]); m0 <- sum(w0*tab$mu0[i])
    s1 <- sum(w1*(tab$v1[i]+tab$mu1[i]^2)) - m1^2
    s0 <- sum(w0*(tab$v0[i]+tab$mu0[i]^2)) - m0^2
    out <- out + pis*(eS*s1 + h(eS)*s0 + eS*(m1-m0-psi)^2)
  }
  out/p^2
}
cond_stats <- function(tab, S) {   # per-cell-of-S: pi_s, e_S, and helpers
  key <- if (length(S)) do.call(paste, tab[, S, drop=FALSE]) else rep(1, nrow(tab))
  list(key = key)
}

set.seed(7); errE <- 0; negE <- 0; nO_best <- 0; N <- 300
for (rep in 1:N) {
  tab <- expand.grid(z1 = 0:1, w = 0:1, i = 0:1)
  pz <- runif(1,.2,.8); pw <- runif(1,.2,.8); pi_ <- runif(1,.2,.8)
  tab$pi <- ifelse(tab$z1==1,pz,1-pz)*ifelse(tab$w==1,pw,1-pw)*ifelse(tab$i==1,pi_,1-pi_)
  ez <- matrix(runif(4,.05,.95),2)                # e(z1, i): confounder + instrument
  tab$e <- ez[cbind(tab$z1+1, tab$i+1)]
  base <- rnorm(4,0,2)                            # mu0(z1,w): both parents
  tab$mu0 <- base[2*tab$z1 + tab$w + 1]
  tau <- rnorm(1)                                 # NO effect modification
  tab$mu1 <- tab$mu0 + tau
  tab$v0 <- runif(8,.5,2)[2*tab$z1+tab$w+1]; tab$v1 <- runif(8,.5,2)[2*tab$z1+tab$w+1]
  O <- c("z1","w")
  VO <- bound_exact(tab, O)
  for (Z in list("z1", c("z1","i"), c("z1","w","i"))) {
    VZ <- bound_exact(tab, Z)
    # RHS of (5): term1 over cells of Z; term2 over cells of O
    p <- sum(tab$pi*tab$e)
    keyZ <- do.call(paste, tab[, Z, drop=FALSE]); keyU <- do.call(paste, tab[, union(Z,O), drop=FALSE])
    keyO <- do.call(paste, tab[, O, drop=FALSE])
    # mu0(Z∪O) per union cell (W ⊥ A | Z holds: e free of w) -> Var(mu0(Z∪O)|Z) under P(w|z)
    t1 <- 0
    for (s in unique(keyZ)) {
      iz <- keyZ == s; wz <- tab$pi[iz]/sum(tab$pi[iz]); eZ <- sum(wz*tab$e[iz])
      # mu0(Z∪O) on union cells inside s: since union adds w only, mu0 = tab$mu0 (given z1,w; free of i)
      # conditional on Z-cell s, mu0(Z∪O) takes value mu0(z1,w) with prob P(w|s)
      ku <- keyU[iz]; mu_u <- tapply(tab$mu0[iz]*wz, ku, sum)/tapply(wz, ku, sum)
      pw_u <- tapply(wz, ku, sum)
      V0 <- sum(pw_u*mu_u^2) - (sum(pw_u*mu_u))^2
      t1 <- t1 + sum(tab$pi[iz]) * eZ/(1-eZ) * V0
    }
    t2 <- 0
    for (s in unique(keyO)) {
      io <- keyO == s; wo <- tab$pi[io]/sum(tab$pi[io]); eO <- sum(wo*tab$e[io])
      ku <- keyU[io]; e_u <- tapply(tab$e[io]*wo, ku, sum)/tapply(wo, ku, sum); p_u <- tapply(wo, ku, sum)
      w0 <- wo*(1-tab$e[io])/(1-eO)
      s0 <- sum(w0*(tab$v0[io]+tab$mu0[io]^2)) - (sum(w0*tab$mu0[io]))^2
      t2 <- t2 + sum(tab$pi[io]) * s0 * (sum(p_u*h(e_u)) - h(eO))
    }
    rhs <- (t1 + t2)/p^2
    errE <- max(errE, abs((VZ - VO) - rhs)); if (VZ - VO < -1e-12) negE <- negE + 1
  }
  Vall <- c(bound_exact(tab,"z1"), VO, bound_exact(tab,c("z1","i")), bound_exact(tab,c("z1","w","i")))
  if (which.min(Vall) == 2) nO_best <- nO_best + 1
}
cat(sprintf("Theorem 3 check on %d random no-modification designs:\n  max |(V(Z)-V(O)) - RHS(5)| = %.3e ;  V(Z) < V(O) count = %d (expect 0) ;  O is argmin in %d/%d\n",
            N, errE, negE, nO_best, N))
