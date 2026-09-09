# Replication archive: Optimal Covariate Adjustment beyond the Average Treatment Effect: Treated-Population and Overlap-Weighted Estimands

Code and (public) data to reproduce the tables and figures of the paper.
Preprint: arXiv:XXXX.XXXXX (to be filled at posting). Author: Shoki Okubo (Toyo University).

## Contents
- analysis
- analysis/data
- analysis/output
- figures

## How to run
1. Install R (>= 4.1) and, where required, `remotes::install_github("sokubo/dagmv")`.
2. Run the scripts in each folder in the order given by their headers; outputs are written to `output/`.
   Where the folder contains its own `analysis/README.md`, that file maps every table, figure and
   script-produced number to a script, an output file and a seed, and records data provenance,
   checksums and software versions.
3. The application uses the public NSW and CPS data included in `analysis/data/`; the simulations generate their own data. See `analysis/README.md` for provenance and checksums.

## Citation
Okubo, S. (2026). Optimal Covariate Adjustment beyond the Average Treatment Effect: Treated-Population and Overlap-Weighted Estimands. Working paper. arXiv:XXXX.XXXXX.
