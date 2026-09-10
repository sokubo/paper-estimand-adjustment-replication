# get_lalonde_data.R — public acquisition / verification step for the two
# data files used in Section 9 (apply_lalonde_att_v2.R):
#   nsw_dw.rda        Dehejia–Wahba experimental subsample of the NSW data
#                     (445 rows: 185 treated, 260 controls; 11 variables)
#   cps_controls.rda  CPS-1 comparison group (15,992 rows; 11 variables)
# Both are Rajeev Dehejia's distribution files nsw_dw.dta and cps_controls.dta
# read into R and saved unchanged (variables data_id, treat, age, education,
# black, hispanic, married, nodegree, re74, re75, re78; earnings in 1982 USD).
# Both are stored UNCOMPRESSED (save(..., compress = FALSE)): arXiv's upload
# processing tries to decompress any file whose bytes look like a compressed
# archive, which corrupts a bzip2- or gzip-compressed .rda in the ancillary
# bundle. load() reads either form, so this affects only the stored bytes.
# SHA-256 of the files shipped with the archive:
#   nsw_dw.rda        f5c210eb2a1caddda736a6decff9a2b0e50a5e4cc58dd42ef8ac3a3058a80757
#   cps_controls.rda  5c770cf1e9c0b780dc58e0de693157839e27f05e0fe38b9837ace613c443a83c
#
# If the .rda files are present this script only verifies them. If they are
# absent it rebuilds them from a public copy, in this order:
#   (1) the CRAN package `causaldata` (Huntington-Klein), datasets nsw_mixtape
#       and cps_mixtape, whose nine numeric columns and data_id were checked
#       to be value-identical to the shipped files (causaldata 0.1.5 Stata
#       files; columns educ/hisp/marr are renamed education/hispanic/married);
#   (2) Dehejia's site, https://users.nber.org/~rdehejia/data/{nsw_dw,cps_controls}.dta,
#       read with haven::read_dta.
# Every path ends with the same checks: row counts, treated counts, and the
# experimental benchmark 1794.34 (difference in mean 1978 earnings).
# Run from analysis/: Rscript data/get_lalonde_data.R
here <- if (dir.exists("data")) "data" else "."
f_nsw <- file.path(here, "nsw_dw.rda"); f_cps <- file.path(here, "cps_controls.rda")
vars <- c("data_id", "treat", "age", "education", "black", "hispanic", "married",
          "nodegree", "re74", "re75", "re78")

check <- function(nsw_dw, cps_controls) {
  stopifnot(identical(names(nsw_dw), vars), identical(names(cps_controls), vars))
  stopifnot(nrow(nsw_dw) == 445, sum(nsw_dw$treat == 1) == 185, sum(nsw_dw$treat == 0) == 260)
  stopifnot(nrow(cps_controls) == 15992, all(cps_controls$treat == 0))
  bench <- mean(nsw_dw$re78[nsw_dw$treat == 1]) - mean(nsw_dw$re78[nsw_dw$treat == 0])
  stopifnot(abs(bench - 1794.34) < 0.01)
  cat(sprintf("checks passed: nsw_dw 185/260, cps_controls 15,992, benchmark %.2f\n", bench))
}

if (file.exists(f_nsw) && file.exists(f_cps)) {
  load(f_nsw); load(f_cps); check(nsw_dw, cps_controls)
  if (requireNamespace("digest", quietly = TRUE)) {
    cat("sha256 nsw_dw.rda      ", digest::digest(f_nsw, algo = "sha256", file = TRUE), "\n")
    cat("sha256 cps_controls.rda", digest::digest(f_cps, algo = "sha256", file = TRUE), "\n")
  } else cat("(install the digest package, or run `shasum -a 256 data/*.rda`, to print the checksums)\n")
  quit(save = "no")
}

as_df <- function(x) { x <- as.data.frame(x); names(x) <- sub("^educ$", "education", names(x))
  names(x) <- sub("^hisp$", "hispanic", names(x)); names(x) <- sub("^marr$", "married", names(x))
  x[vars] }
if (requireNamespace("causaldata", quietly = TRUE)) {
  cat("rebuilding from CRAN package causaldata\n")
  nsw_dw <- as_df(causaldata::nsw_mixtape); cps_controls <- as_df(causaldata::cps_mixtape)
} else if (requireNamespace("haven", quietly = TRUE)) {
  cat("rebuilding from https://users.nber.org/~rdehejia/data/\n")
  base <- "https://users.nber.org/~rdehejia/data/"
  nsw_dw <- as_df(haven::read_dta(paste0(base, "nsw_dw.dta")))
  cps_controls <- as_df(haven::read_dta(paste0(base, "cps_controls.dta")))
} else stop("install either the causaldata or the haven package to rebuild the data files")
for (v in vars[-1]) { nsw_dw[[v]] <- as.numeric(nsw_dw[[v]]); cps_controls[[v]] <- as.numeric(cps_controls[[v]]) }
nsw_dw$data_id <- as.character(nsw_dw$data_id); cps_controls$data_id <- as.character(cps_controls$data_id)
check(nsw_dw, cps_controls)
save(nsw_dw, file = f_nsw, compress = FALSE); save(cps_controls, file = f_cps, compress = FALSE)
cat("wrote", f_nsw, "and", f_cps, "(rebuilt copies; checksums may differ from the shipped files, values are identical)\n")
