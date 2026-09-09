# fig_signmap.R — Figure 1: sign of the outcome-side increment (Theorem 4)
# for arm-specific predictors, by estimand and propensity e. A sign map.
suppressMessages(library(ggplot2))
e <- seq(.005, .995, by = .005)
sgn <- function(v) ifelse(abs(v) < 1e-12, "neutral", ifelse(v > 0, "harms", "helps"))
rowsf <- function(lab, h, hp) rbind(
  data.frame(estimand = lab, arm = "treated-arm-only predictor", e = e,
             s = sgn(hp(e)^2*e*(1-e) - h(e)^2*(1-e)/e)),
  data.frame(estimand = lab, arm = "control-arm-only predictor", e = e,
             s = sgn(hp(e)^2*e*(1-e) - h(e)^2*e/(1-e))))
d <- rbind(rowsf("ATE   h(e) = 1",      function(e) 1+0*e,   function(e) 0*e),
           rowsf("ATT   h(e) = e",      function(e) e,       function(e) 1+0*e),
           rowsf("ATO   h(e) = e(1-e)", function(e) e*(1-e), function(e) 1-2*e))
d$estimand <- factor(d$estimand, levels = rev(unique(d$estimand)))
d$arm <- factor(d$arm, levels = c("treated-arm-only predictor", "control-arm-only predictor"))
d$s <- factor(d$s, levels = c("helps", "neutral", "harms"))
thr <- data.frame(arm = factor(c("control-arm-only predictor","control-arm-only predictor","treated-arm-only predictor"),
                               levels = levels(d$arm)),
                  x = c(.5, 1/3, 2/3), lab = c("1/2 (ATT)", "1/3 (ATO)", "2/3 (ATO)"))
fig <- ggplot(d, aes(x = e, y = estimand, fill = s)) +
  geom_tile(height = .62, width = .006) +
  geom_vline(data = thr, aes(xintercept = x), linetype = "dotted", color = "#0b0b0b", linewidth = .45) +
  geom_text(data = thr, aes(x = x, y = c(3.62, 3.62, 3.62) - c(0, 0, 0), label = lab, hjust = c(1.05, 1.05, -0.05)), inherit.aes = FALSE, size = 2.9, color = "#52514e") +
  facet_wrap(~arm, nrow = 1) +
  scale_fill_manual(values = c(helps = "#e9edf3", neutral = "#c9c7c0", harms = "#3d3c39"), name = "adding the predictor") +
  scale_x_continuous(breaks = c(0, 1/3, .5, 2/3, 1), labels = c("0", "1/3", "1/2", "2/3", "1"), expand = c(.01, 0)) +
  coord_cartesian(ylim = c(.5, 3.9)) +
  labs(x = expression(propensity ~ e[S]), y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "#c9c7c0", linewidth = .4),
        strip.text = element_text(face = "bold", color = "#52514e"),
        plot.background = element_rect(fill = "white", color = NA))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig1_signmap.png", fig, width = 8, height = 2.9, dpi = 220, bg = "white")
cat("fig1 written\n")
