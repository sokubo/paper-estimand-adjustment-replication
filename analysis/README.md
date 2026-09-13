# Replication code for "Optimal Covariate Adjustment beyond the Average Treatment Effect"

Archive version 1.1 (September 2026), matching arXiv version 2 of the paper,
at https://github.com/sokubo/paper-estimand-adjustment-replication.

Changes from version 1.0 (arXiv v1): `counterexample_search.R` now prints,
next to the divergence count (structures whose ATE- and ATT-optimal sets
differ with an exact ATT penalty above 1 %), the number of structures with
distinct minimizers and *any* positive penalty, and accepts a scan seed
(`SCAN_SEED=<n> Rscript counterexample_search.R`); outputs for seeds 123 and
2026 are added (`output/counterexample_results_seed123.txt`,
`output/counterexample_results_seed2026.txt`). At the paper's seed 99 the two
counts coincide (10 = 10) and every shipped number is unchanged. The
threshold was applied in version 1.0 but not stated in the paper; we thank
Jürgen Degenfellner (ZHAW) for reproducing the archive independently and
pointing this out.

Every table entry, the figure, and every number in the paper attributed to a
script are produced by one of the scripts below (closed-form constants worked
by hand in the text — for example the 3/533 of Proposition 1 or the 5 versus
8.5 of the Section 6 counterexample — are not). All scripts are plain R (base R plus
`dagmv` for the random-structure scan and `ggplot2` for Figure 1); none needs
more than a minute on a laptop. Run them from this directory (`analysis/`)
unless stated otherwise; outputs are written to `output/` (already included,
so that every script-produced number can be checked without rerunning anything).

## Table / figure / number → script → output

| Where in the paper | Script (entry point) | Output file(s) | Seed |
|---|---|---|---|
| Section 5, Lemma 3 instance ($V_{att}(\{Z\}) = 6.433$ vs $7.360$); Appendix B parts A–D (closed forms vs exact engine, identities (3) and (5), deletion inequality, saturated AIPW at $n = 20{,}000$) | `verify_theorems.R` | `output/verify_theorems.txt` (stdout, redirected) | 1, 2, 3 |
| Appendix B, Theorem 3 check (300 no-modification designs) | `verify_thm3.R` | `output/verify_thm3.txt` | 7 |
| Appendix B, Lemma 2 / identity (6) / Corollary 6 thresholds / Proposition 1 drift (7) | `verify_wate.R` | `output/verify_wate.txt` | 11, 12 |
| Section 8 last paragraph and Appendix B: random-structure scan (49 structures, 10 divergences = 10 structures with any positive penalty, penalties 2.92–25.86 %, median 14.73 %, prevalence 0.190–0.559, 21 structures with exact ties) | `counterexample_search.R` (needs `dagmv`) | `output/counterexample_results.txt` | 20260818 (design checks T1–T4, not reported), 99 (the scan) |
| Appendix B, supplementary scan seeds (seed 123: 53 structures, 5 with positive penalty of which 4 above 1 %, penalties 0.35, 16.47, 18.99, 24.27, 24.57 %; seed 2026: 56 structures, 6 with positive penalty of which 5 above 1 %, penalties 0.01, 1.49, 4.94, 10.69, 18.21, 26.87 %) | `SCAN_SEED=123 Rscript counterexample_search.R`, `SCAN_SEED=2026 Rscript counterexample_search.R` | `output/counterexample_results_seed123.txt`, `output/counterexample_results_seed2026.txt` | 123, 2026 |
| Table 2 (finite-sample simulation, designs A and B, $n = 500, 2000$, 1,000 replications) | `sim_finite_sample_v2.R` (`Rscript sim_finite_sample_v2.R [R] [seed]`) | `output/finite_sample_v2_table.csv`, `output/finite_sample_v2_summary.txt` | 2026 (bootstrap streams 2026 + 1000 n + 7) |
| Table 3 and the Section 9 diagnostics (LaLonde ATT, 185 NSW trainees vs 15,992 CPS-1 controls) | `apply_lalonde_att_v2.R` | `output/lalonde_att_v2_table.csv` (columns: `est`, `se` in thousands of 1982 dollars; `se_change_pct`, `est_change`; `ctrl_pred_sensitivity_pct` = Table 3 "Control-prediction sensitivity"; `t_treated`, `t_control`, `z_ps` = absolute t / Wald z statistics), `output/lalonde_att_v2_summary.txt` (propensity diagnostics, the 0.64 maximum control propensity, the benchmark 1.794 with unpooled SE 0.671) | none (deterministic) |
| Figure 1 (sign map) | `fig_signmap.R` — run from the directory above (`Rscript analysis/fig_signmap.R`), it writes `figures/fig1_signmap.png` | `../figures/fig1_signmap.png` | none |

The `verify_*.R` scripts print to standard output; the shipped `.txt` files
are those printouts (`Rscript verify_theorems.R > output/verify_theorems.txt`).

## Data (Section 9 only)

`data/nsw_dw.rda` (445 rows: 185 treated, 260 experimental controls) and
`data/cps_controls.rda` (15,992 rows) are Rajeev Dehejia's distribution files
`nsw_dw.dta` and `cps_controls.dta` (the Dehejia–Wahba 1999 experimental
subsample and LaLonde's CPS-1 comparison group; earnings in 1982 dollars),
read into R and saved unchanged. SHA-256:

```
6b4637bdadac7e4c7d95e9d2a1a8a3afcdf0de8b15cb6bb47731d6efcb25d49e  data/nsw_dw.rda
0c41b180b79c01eb77f201a24b5366a205b72f0b7bdabc822fbee505fd2a5a23  data/cps_controls.rda
```

`data/get_lalonde_data.R` verifies the two files (row counts, treated counts,
the experimental benchmark 1,794.34) and prints their checksums; if they are
missing it rebuilds them from the CRAN package `causaldata` (datasets
`nsw_mixtape`, `cps_mixtape`) or from Dehejia's site. The `causaldata` copies
were checked to be value-identical to the shipped files in every variable, and
`apply_lalonde_att_v2.R` run on files rebuilt from them reproduces
`output/lalonde_att_v2_table.csv` exactly. No other data enter the paper; the
simulations generate their own.

Sample filters and variable construction (all in `apply_lalonde_att_v2.R`):
the analysis sample is the 185 treated rows of `nsw_dw` stacked on the 15,992
rows of `cps_controls` ($n = 16{,}177$); the ten covariates are age,
education, black, hispanic, married, nodegree, re74/1000, re75/1000,
u74 = 1{re74 = 0}, u75 = 1{re75 = 0}; the outcome is re78/1000; the propensity
is a logistic regression on the ten covariates entered linearly, truncated at
0.9 (the truncation never binds); the control-outcome regression is linear on
the same covariates fitted to the 15,992 controls.

## Software versions used for the shipped outputs

R 4.3.3 (2024-02-29); `dagmv` 0.1.2 (`remotes::install_github("sokubo/dagmv")`,
then check `packageVersion("dagmv")`; the package is used only to enumerate
the sets satisfying the adjustment criterion in `counterexample_search.R`);
`ggplot2` 3.4.4 (Figure 1 only). Everything else is base R. Rerunning the
seeded scripts under R 4.3.3 reproduces the shipped output files exactly
(checked from a clean copy of this directory). An independent rerun of the
archive under R 4.6.0 (Jürgen Degenfellner, September 2026) reproduced the
reported results; other R versions may still change the last floating-point
digits of the machine-precision checks, but the tables are unaffected at
their reported rounding.

## Conventions worth knowing

* `sim_finite_sample_v2.R`: a covariate cell with no treated or no control
  unit is dropped from that replication's estimator (the replication is kept).
  With seed 2026 this happens in 2 of 1,000 replications of design A at
  $n = 500$ for the set $\{Z, W, I\}$ (control-only cells of 62 and 59
  observations); the point estimate is unchanged and the standard error
  differs from retaining zero contributions by about 0.014 %.
* `counterexample_search.R`: all four covariates are declared as vertices of
  each drawn DAG; valid sets are those satisfying the adjustment criterion
  relative to the drawn graph; bounds are exact (enumeration of the $2^4$
  cells); the 21 structures with exactly tied minimizers (relative tolerance
  $10^{-9}$) are listed in full at the end of the output file, and the
  printed "largest relative spread" line confirms that sets tied for one
  minimum share the other bound (relative spread $2.2 \times 10^{-15}$), so
  the penalty does not depend on which tied set is named. A structure is
  *counted* as a divergence only when the exact ATT penalty of the
  ATE-optimal set exceeds 1 % (`gapATT > 0.01`); every positive penalty
  (`gapATT > 1e-9`, the tie tolerance) is nevertheless listed, and the line
  "structures with distinct minimizers and any positive ATT penalty" gives
  the unthresholded count. At seed 99 both counts are 10; at seeds 123 and
  2026 one structure each has a penalty below 1 % (0.35 % and 0.01 %).

## Legacy scripts (kept for the record, not used by the paper)

`sim_finite_sample.R` and `apply_lalonde_att.R` are the first versions of the
simulation and the application (outputs `output/finite_sample_*`,
`output/lalonde_att_summary.txt`, `output/lalonde_att_table.csv`);
`counterexample_search_v1_mc.R` is the Monte Carlo-only version of the scan
with an edge-only graph declaration (see Appendix B). The paper's tables and
numbers come from the `_v2` scripts and from `counterexample_search.R`.
