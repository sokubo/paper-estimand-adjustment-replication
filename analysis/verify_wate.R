# verify_wate.R — Paper 4 §WATE: verification of (i) Lemma 2 (WATE bound
# with propensity-dependent weight h(e)) against finite-sample Monte Carlo,
# (ii) Theorem 4 (WATE supplementation identity) exactly, (iii) Corollary 6
# ATO thresholds (1/3, 2/3), (iv) Proposition 1 (ATO estimand drift under
# instrument adjustment) exactly.
options(digits = 10)
H <- list(ate = list(h = function(e) 1 + 0*e,      hp = function(e) 0*e),
          att = list(h = function(e) e,            hp = function(e) 1 + 0*e),
          ato = list(h = function(e) e*(1-e),      hp = function(e) 1 - 2*e))

# exact cell engine for the WATE bound V_h(S) on a probability table
bound_h <- function(tab, S, w) {
  key <- if (length(S)) do.call(paste, tab[, S, drop=FALSE]) else rep(1, nrow(tab))
  cells <- unique(key); K <- length(cells)
  pis <- eS <- m1 <- m0 <- s1 <- s0 <- numeric(K)
  for (k in seq_len(K)) {
    i <- key == cells[k]; pis[k] <- sum(tab$pi[i]); wc <- tab$pi[i]/pis[k]
    eS[k] <- sum(wc*tab$e[i]); w1 <- wc*tab$e[i]/eS[k]; w0 <- wc*(1-tab$e[i])/(1-eS[k])
    m1[k] <- sum(w1*tab$mu1[i]); m0[k] <- sum(w0*tab$mu0[i])
    s1[k] <- sum(w1*(tab$v1[i]+tab$mu1[i]^2)) - m1[k]^2
    s0[k] <- sum(w0*(tab$v0[i]+tab$mu0[i]^2)) - m0[k]^2
  }
  h <- w$h(eS); hp <- w$hp(eS); tau <- m1 - m0
  D <- sum(pis*h); psi <- sum(pis*h*tau)/D
  V <- sum(pis*(h^2*s1/eS + h^2*s0/(1-eS) + (tau-psi)^2*(h^2 + hp^2*eS*(1-eS))))/D^2
  c(V = V, psi = psi)
}

cat("== (i) Lemma 2: ATO bound vs finite-sample Monte Carlo (saturated plug-in) ==\n")
set.seed(11); n <- 40000; R <- 400
tab <- expand.grid(z=0:1, w=0:1); tab$pi <- c(.3,.3,.2,.2)*1; tab$pi <- tab$pi/sum(tab$pi)
tab$e <- c(.2,.35,.55,.7); tab$mu0 <- c(0,1,.5,2); tab$mu1 <- c(1,1.5,2.5,2); tab$v0 <- 1; tab$v1 <- 1.5
for (nm in c("att","ato")) {
  w <- H[[nm]]; est <- numeric(R)
  for (r in 1:R) {
    cell <- sample(1:4, n, TRUE, tab$pi); e <- tab$e[cell]; a <- rbinom(n,1,e)
    y <- ifelse(a==1, tab$mu1[cell] + rnorm(n, 0, sqrt(tab$v1)), tab$mu0[cell] + rnorm(n, 0, 1))
    eh <- tapply(a, cell, mean); m1 <- tapply(y[a==1], cell[a==1], mean); m0 <- tapply(y[a==0], cell[a==0], mean)
    pc <- tabulate(cell, 4)/n; hh <- w$h(eh)
    est[r] <- sum(pc*hh*(m1-m0))/sum(pc*hh)
  }
  b <- bound_h(tab, c("z","w"), w)
  cat(sprintf("  %s: n*Var_MC = %.4f  bound = %.4f  ratio = %.3f  (mean est %.4f vs psi %.4f)\n",
              nm, n*var(est), b["V"], n*var(est)/b["V"], mean(est), b["psi"]))
}

cat("\n== (ii) Theorem 4: WATE supplementation identity (W ⊥ A | Z), 300 random designs ==\n")
set.seed(12); err <- c(att=0, ato=0, ate=0)
for (rep in 1:300) {
  tab <- expand.grid(z=0:1, w=0:1); pz <- runif(1,.2,.8); pw <- runif(2,.2,.8)
  tab$pi <- ifelse(tab$z==1,pz,1-pz)*ifelse(tab$w==1,pw[tab$z+1],1-pw[tab$z+1])
  ez <- runif(2,.05,.95); tab$e <- ez[tab$z+1]
  tab$mu0 <- rnorm(4,0,2); tab$mu1 <- rnorm(4,0,2); tab$v0 <- runif(4,.5,2); tab$v1 <- runif(4,.5,2)
  for (nm in names(H)) {
    w <- H[[nm]]
    lhs <- bound_h(tab, c("z","w"), w)["V"] - bound_h(tab, "z", w)["V"]
    D <- sum(tab$pi*w$h(tab$e)); rhs <- 0
    for (z in 0:1) {
      i <- tab$z==z; wz <- tab$pi[i]/sum(tab$pi[i]); e <- ez[z+1]; h <- w$h(e); hp <- w$hp(e)
      m0 <- sum(wz*tab$mu0[i]); m1 <- sum(wz*tab$mu1[i])
      V0 <- sum(wz*tab$mu0[i]^2)-m0^2; V1 <- sum(wz*tab$mu1[i]^2)-m1^2; C <- sum(wz*tab$mu1[i]*tab$mu0[i])-m1*m0
      rhs <- rhs + sum(tab$pi[i])*( V1*(hp^2*e*(1-e) - h^2*(1-e)/e) + V0*(hp^2*e*(1-e) - h^2*e/(1-e)) - 2*C*(h^2 + hp^2*e*(1-e)) )
    }
    err[nm] <- max(err[nm], abs(lhs - rhs/D^2))
  }
}
print(signif(err, 3))

cat("\n== (iii) Corollary 6: ATO thresholds — sign of increment vs e for arm-specific predictors ==\n")
for (e0 in c(.2, .3, .34, .5, .66, .7, .8)) {
  tab <- expand.grid(z=0, w=0:1); tab$pi <- c(.5,.5); tab$e <- e0
  # control-arm-only predictor
  tab$mu0 <- c(0,1); tab$mu1 <- c(1,1); tab$v0 <- tab$v1 <- 1
  dc <- bound_h(tab, "w", H$ato)["V"] - bound_h(tab, character(0), H$ato)["V"]
  # treated-arm-only predictor
  tab$mu0 <- c(0,0); tab$mu1 <- c(1,2)
  dt <- bound_h(tab, "w", H$ato)["V"] - bound_h(tab, character(0), H$ato)["V"]
  cat(sprintf("  e=%.2f  control-only: %+8.4f (theory: hurts iff e<1/3)   treated-only: %+8.4f (theory: hurts iff e>2/3)\n", e0, dc, dt))
}

cat("\n== (iv) Proposition 1: estimand drift under instrument adjustment ==\n")
# W pure instrument: e depends on (z,w); mu, v free of w. tau varies with z.
# (assign by explicit (z,w) indexing — expand.grid varies z fastest)
tab <- expand.grid(z=0:1, w=0:1); tab$pi <- rep(.25,4)
tab$e   <- with(tab, ifelse(z==0, ifelse(w==0,.10,.50), ifelse(w==0,.30,.70)))  # e_S: .30 (z=0), .50 (z=1)
tab$mu0 <- ifelse(tab$z==0, 0, 1); tab$mu1 <- ifelse(tab$z==0, 1, 3); tab$v0 <- tab$v1 <- 1   # tau: 1, 2
for (nm in names(H)) {
  a <- bound_h(tab, "z", H[[nm]]); b <- bound_h(tab, c("z","w"), H[[nm]])
  cat(sprintf("  %s: psi(S={Z}) = %.5f  psi(S={Z,W}) = %.5f  drift = %+.5f\n", nm, a["psi"], b["psi"], b["psi"]-a["psi"]))
}
# closed-form drift for ATO: E[{E[h(e_SW)|S]-h(e_S)}(tau_S-psi_S)]/E[h(e_SW)]
hS <- function(e) e*(1-e)
EhSW <- c(mean(hS(c(.10,.50))), mean(hS(c(.30,.70)))); hSz <- hS(c(.30,.50)); tau <- c(1,2)
psiS <- sum(.5*hSz*tau)/sum(.5*hSz); D2 <- sum(.5*EhSW)
cat(sprintf("  closed-form ATO drift = %+.5f\n", sum(.5*(EhSW-hSz)*(tau-psiS))/D2))
