# apply_lalonde_att.R — Paper 4 empirical illustration: the ATT of NSW
# training (185 treated) with 15,992 CPS controls — the canonical ATT
# setting, with treatment prevalence ~1.1% (deep in the e < 1/2 regime).
# For each leave-one-covariate-out adjustment set we report the AIPW-ATT
# estimate, its EIF standard error, and arm-specific relevance diagnostics,
# to confront the calculus of Theorem 1 with data.
load("data/nsw_dw.rda"); load("data/cps_controls.rda")
tr <- nsw_dw[nsw_dw$treat == 1, ]; co <- cps_controls
d <- rbind(tr, co); d$u74 <- as.integer(d$re74 == 0); d$u75 <- as.integer(d$re75 == 0)
d$re74k <- d$re74/1000; d$re75k <- d$re75/1000; d$re78k <- d$re78/1000
X_full <- c("age","education","black","hispanic","married","nodegree","re74k","re75k","u74","u75")
n <- nrow(d); A <- d$treat; Y <- d$re78k; p <- mean(A)
cat(sprintf("n = %d, treated = %d, prevalence = %.4f\n", n, sum(A), p))

aipw_att <- function(S, trim = 0.9) {
  f_ps <- as.formula(paste("treat ~", paste(S, collapse = "+")))
  e <- fitted(glm(f_ps, data = d, family = binomial()))
  e <- pmin(e, trim)
  f_om <- as.formula(paste("re78k ~", paste(S, collapse = "+")))
  m0 <- predict(lm(f_om, data = d[A == 0, ]), newdata = d)
  lam <- e/(1-e)
  phi0 <- (A*(Y - m0) - lam*(1-A)*(Y - m0))/p
  psi <- mean(phi0); phi <- phi0 - A*psi/p
  c(est = psi, se = sd(phi)/sqrt(n), max_e = max(e[A==0]))
}
# arm-specific relevance: |t| of each covariate in the treated-arm and
# control-arm outcome regressions and in the propensity logit (full set)
t_arm <- function(sub) { m <- lm(as.formula(paste("re78k ~", paste(X_full, collapse="+"))), data = sub)
  abs(summary(m)$coefficients[X_full, "t value"]) }
t1 <- t_arm(d[A==1,]); t0 <- t_arm(d[A==0,])
tps <- abs(summary(glm(as.formula(paste("treat ~", paste(X_full, collapse="+"))), data=d, family=binomial()))$coefficients[X_full, "z value"])

full <- aipw_att(X_full)
rows <- data.frame(set = "full (all 10)", dropped = "", est = full["est"], se = full["se"],
                   t_treated = NA, t_control = NA, z_ps = NA)
for (v in X_full) {
  r <- aipw_att(setdiff(X_full, v))
  rows <- rbind(rows, data.frame(set = paste("drop", v), dropped = v, est = r["est"], se = r["se"],
                                 t_treated = t1[v], t_control = t0[v], z_ps = tps[v]))
}
rows$se_change_pct <- 100*(rows$se/full["se"] - 1)
rows$est_change <- rows$est - full["est"]
rownames(rows) <- NULL
print(rows, digits = 3)
cat(sprintf("\nExperimental benchmark (Dehejia-Wahba experimental subsample): 1.794 thousand 1982 dollars (1978 earnings), unpooled SE 0.671\n"))
dir.create("output", showWarnings = FALSE)
write.csv(rows, "output/lalonde_att_table.csv", row.names = FALSE)
capture.output(print(rows, digits=3), file = "output/lalonde_att_summary.txt")
