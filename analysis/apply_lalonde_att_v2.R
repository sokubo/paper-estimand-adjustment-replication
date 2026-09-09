# apply_lalonde_att_v2.R -- Paper 4 empirical illustration, version 2.
# NSW treated (185) with CPS-1 controls (15,992): AIPW-ATT under the full
# ten-covariate set and under each single-covariate deletion, with the
# diagnostics used in Section 9 of the revision:
#   * fitted propensities of the treated units (median, mean, quartiles) and
#     the largest control propensity (whether the 0.9 truncation binds);
#   * the control-prediction sensitivity measure of Table 3: the mean over
#     treated units of {m0_full(X) - m0_reduced(X)}^2, the squared change in
#     the fitted control-outcome regression when the covariate is dropped,
#     as a percentage of the treated-arm term var(Y - m0_full | A = 1). This
#     is a regression-difference sensitivity measure, NOT an estimate of the
#     conditional-variance sign functional of Theorem 1 / Remark 4 (which
#     would require the conditional variance of mu_0 given the reduced set);
#   * the experimental benchmark: the Dehejia-Wahba experimental subsample
#     difference in means with its UNPOOLED difference-in-means standard
#     error (the published SE 633 uses the pooled convention).
# Validity of the deletions is not established by any of these numbers; the
# table describes the estimate and its standard error under each set.
load("data/nsw_dw.rda"); load("data/cps_controls.rda")
tr <- nsw_dw[nsw_dw$treat == 1, ]; co <- cps_controls
d <- rbind(tr, co); d$u74 <- as.integer(d$re74 == 0); d$u75 <- as.integer(d$re75 == 0)
d$re74k <- d$re74 / 1000; d$re75k <- d$re75 / 1000; d$re78k <- d$re78 / 1000
X_full <- c("age", "education", "black", "hispanic", "married", "nodegree", "re74k", "re75k", "u74", "u75")
n <- nrow(d); A <- d$treat; Y <- d$re78k; p <- mean(A)
cat(sprintf("n = %d, treated = %d, mixing proportion = %.4f\n", n, sum(A), p))

fit_parts <- function(S, trim = 0.9) {
  f_ps <- as.formula(paste("treat ~", paste(S, collapse = "+")))
  e_raw <- fitted(glm(f_ps, data = d, family = binomial()))
  e <- pmin(e_raw, trim)
  f_om <- as.formula(paste("re78k ~", paste(S, collapse = "+")))
  m0 <- predict(lm(f_om, data = d[A == 0, ]), newdata = d)
  lam <- e / (1 - e)
  phi0 <- (A * (Y - m0) - lam * (1 - A) * (Y - m0)) / p
  psi <- mean(phi0); phi <- phi0 - A * psi / p
  list(est = psi, se = sd(phi) / sqrt(n), e_raw = e_raw, m0 = m0,
       max_e_control = max(e_raw[A == 0]), n_trunc_control = sum(e_raw[A == 0] > trim))
}
full <- fit_parts(X_full)
e_tr <- full$e_raw[A == 1]
cat(sprintf("treated propensities (full set): median %.3f, mean %.3f, IQR %.3f-%.3f, share > 0.5 = %.2f, share > 0.1 = %.2f\n",
            median(e_tr), mean(e_tr), quantile(e_tr, .25), quantile(e_tr, .75), mean(e_tr > .5), mean(e_tr > .1)))
cat(sprintf("largest control propensity %.4f; controls above the 0.9 truncation: %d (truncation does not bind)\n",
            full$max_e_control, full$n_trunc_control))
treated_term <- var((Y - full$m0)[A == 1])
cat(sprintf("treated-arm term var(Y - m0_full | A = 1) = %.2f ($1000^2); residual SD %.2f ($1000)\n", treated_term, sqrt(treated_term)))

# arm-specific relevance in the full-set regressions (descriptive)
t_arm <- function(sub) { m <- lm(as.formula(paste("re78k ~", paste(X_full, collapse = "+"))), data = sub)
  abs(summary(m)$coefficients[X_full, "t value"]) }
t1 <- t_arm(d[A == 1, ]); t0 <- t_arm(d[A == 0, ])
tps <- abs(summary(glm(as.formula(paste("treat ~", paste(X_full, collapse = "+"))), data = d, family = binomial()))$coefficients[X_full, "z value"])

rows <- data.frame(set = "full (all 10)", dropped = "", est = full$est, se = full$se,
                   se_change_pct = 0, est_change = 0, ctrl_pred_shift_treated = NA, ctrl_pred_sensitivity_pct = NA,
                   t_treated = NA, t_control = NA, z_ps = NA)
for (v in X_full) {
  r <- fit_parts(setdiff(X_full, v))
  inc <- mean(((full$m0 - r$m0)[A == 1])^2)     # squared shift of the fitted control regression on the treated (regression-difference measure)
  rows <- rbind(rows, data.frame(set = paste("drop", v), dropped = v, est = r$est, se = r$se,
                                 se_change_pct = 100 * (r$se / full$se - 1), est_change = r$est - full$est,
                                 ctrl_pred_shift_treated = inc, ctrl_pred_sensitivity_pct = 100 * inc / treated_term,
                                 t_treated = t1[v], t_control = t0[v], z_ps = tps[v]))
}
rownames(rows) <- NULL
print(rows, digits = 3)

# experimental benchmark from the NSW sample itself
bm <- nsw_dw
y1 <- bm$re78[bm$treat == 1] / 1000; y0 <- bm$re78[bm$treat == 0] / 1000
bench <- c(est = mean(y1) - mean(y0), se = sqrt(var(y1) / length(y1) + var(y0) / length(y0)), n1 = length(y1), n0 = length(y0))
cat(sprintf("\nExperimental benchmark (Dehejia-Wahba experimental subsample, %d treated vs %d controls): %.3f (unpooled SE %.3f) thousand 1982 dollars (1978 earnings)\n",
            bench["n1"], bench["n0"], bench["est"], bench["se"]))
dir.create("output", showWarnings = FALSE)
write.csv(rows, "output/lalonde_att_v2_table.csv", row.names = FALSE)
sink("output/lalonde_att_v2_summary.txt")
cat(sprintf("n = %d, treated = %d, mixing proportion = %.4f\n", n, sum(A), p))
cat(sprintf("treated propensities (full set): median %.3f, mean %.3f, IQR %.3f-%.3f, share > 0.5 = %.2f, share > 0.1 = %.2f\n",
            median(e_tr), mean(e_tr), quantile(e_tr, .25), quantile(e_tr, .75), mean(e_tr > .5), mean(e_tr > .1)))
cat(sprintf("largest control propensity %.4f; controls above the 0.9 truncation: %d\n", full$max_e_control, full$n_trunc_control))
cat(sprintf("treated-arm term var(Y - m0_full | A = 1) = %.2f; residual SD %.2f\n", treated_term, sqrt(treated_term)))
print(rows, digits = 3)
cat(sprintf("\nExperimental benchmark (Dehejia-Wahba experimental subsample, %d treated vs %d controls): %.3f (unpooled SE %.3f) thousand 1982 dollars (1978 earnings)\n", bench["n1"], bench["n0"], bench["est"], bench["se"]))
sink()
