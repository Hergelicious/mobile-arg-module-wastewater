# Supplementary Code

Code for: *Mobile antibiotic resistance genes form a reproducible, integron-associated module in wastewater*

## Input

Additional inputs for Section 2.11: from Zenodo record 14652833, `panres_counts.csv` and `kingdom_motus_pad_agg.csv`; the NCBI run table for PRJNA509305 (`SraRunTable.csv`); and national antibiotic consumption (WHO GLASS, via Our World in Data: `antibiotic-consumption-rate.csv`).


From Zhu et al. (2025), *Nat. Commun.* 16, 4006: the Source Data workbook (`41467_2025_59019_MOESM8_ESM.xlsx`), Supplementary Data 1 (`MOESM3`, sequencing output) and Supplementary Data 3 (`MOESM5`, MGE abundances). From Martiny et al. (2025), *Nat. Commun.* 16, 10278: Supplementary Data 3 and 4 (`41467_2025_66070_MOESM6_ESM.xlsx`, `MOESM7`; network nodes and links). None is redistributed here; download them from the publishers and check their licences. The scripts expect it at `/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx`; edit that path at the top of each script if needed. No other data are used.

## Environment

Python 3.12 with numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib, openpyxl. Node.js with the `docx` package (manuscript build only). Figures use the Liberation Sans font (metrically identical to Arial).

## Run order

Each step reads only files written by earlier steps. Random seeds are fixed, so all statistics are reproducible; the only non-deterministic element is the horizontal jitter of points in some figures, which is cosmetic.

| Step | Script | Writes | Used for |
|---|---|---|---|
| 1 | `analysis.py` | `results.json`, `gene_table.csv`, `sample_layers.csv`, `city_table.csv` | Sample counts, load statistics, ICC, socio-economic models, MAG host breadth, mobility indices |
| 2 | `make_gene_table2.py` | `gene_table2.csv` | Adds gene abundance, variability and detection columns |
| 3 | `controls.py` | `controls.json` | MAG representation and class × mechanism mixing (Section 3.7, Note S1) |
| 4 | `canonical.py` | `canonical.json` | All primary statistics: D (9,999 permutations), controls, bootstrap CI, jackknife, family exclusion, sum-of-squares shares, layer variability, REML |
| 5 | `arrays.py` | `gene_table_final.csv`, `Dnull.npy`, `rhonull.npy` | Per-gene results (Table S1) and permutation nulls for figures |
| 6 | `rf_fold.py` | adds `transfer_R2` to `canonical.json` | Random forests with imputation fitted inside folds (Fig. 1c) |
| 7 | `decode_final.py` | `decode.json`, `rand_acc.npy` | Country decoding (Fig. 4b) |
| 8 | `module.py` | `module.json`, `module_membership.csv` | Label-free module discovery; design fixed in the file header before running (Section 3.4, Fig. 3) |
| 9 | `effort.py` | `effort.json` | Sampling-effort adjustment (Section 3.6) |
| 9b | `freeze_spec.py` | `frozen_spec.json` (SHA-256 recorded) | Pre-specified tests, frozen with a recorded hash before the test script was written (Section 2.10) |
| 9c | `run_frozen.py` | `frozen_results.json` | Runs the frozen tests; refuses to run if the specification's hash has changed |
| 9d | `intI1_pergene.py` | `intI1_rho.csv` | Per-gene intI1 correlations for Table S1 (checked against the frozen run) |
| 9f | `freeze_spec2.py` | `frozen_spec2.json` (SHA-256 recorded) | Second frozen specification: external ResFinder classification test (Section 2.10, test C) |
| 9g | `run_frozen2.py` | `frozen_results2.json` | Runs test C; refuses to run if the specification's hash has changed |
| 9h | `build_martiny.py` | `martiny_matrix.pkl` | Sample-by-gene matrix from Martiny et al. `panres_counts.csv` (Zenodo 14652833) |
| 9i | `freeze_spec3.py` | `frozen_spec3.json` (SHA-256 recorded) | Third frozen specification: external validation, temporal stability, drivers, cross-compartment, batch (Section 2.11) |
| 9j | `run_frozen3.py`, then `finish_frozen3.py` | `frozen_results3.json` | Runs tests D and E, then variance components and test G; the CLR deviation is documented in the script header |
| 9k | `make_bacterial_reference.py` | `bacterial_reference.csv` | Bacterial mOTU fragments per sewage sample, for the pre-specified ALR normalisation |
| 9l | `run_frozen3_alr.py` | `frozen_results3_alr.json`, `F_country_table.csv` | Tests D, E, F and G with the pre-specified ALR normalisation (primary; CLR results are the sensitivity analysis) |
| 9m | `f_specificity.py` | `f_specificity.json` | Exploratory, after test F: consumption association for other acquired and latent genes |
| 9n | `test_H.py` | `test_H.json` | Test H (run-level adjustment) plus the collection-date analyses added once dates were found |
| 9e | `integron_depth.py` | `integron_depth.json`, `intI1_corr.npy` | Exploratory, added after the frozen tests: partial intI1 and qacEΔ1 correlations and the non-annotated module comparison |
| 10 | `figs2.py` | Figs 1, 2, 4, 5, 6 and S1 (PNG, PDF, TIFF in `figout/`) | |
| 11 | `fig_module.py` | Fig. 3 | |
| 12 | `ga.py` | graphical abstract | |
| 13 | `build_tables.py`, then `make_xlsx.py` | `Supplementary_Tables_S1_S2.xlsx` | Tables S1 and S2 |
| 14 | `node build15.js`, `node build_supp.js`, `node hl.js`, `node cover.js` | manuscript, supplementary information, highlights, cover letter | Every number in the text is read from the JSON files above, not typed by hand |

`build15.js` is the concatenation `head.js + body15.js + tail.js`; `build_supp.js` is `head.js + suppbody.js + tail_supp.js`.

## Notes

- `freeze_spec.py` and `freeze_spec2.py` regenerate the frozen specifications; re-running them should produce byte-identical files (same SHA-256). The shipped `frozen_spec.json` and `frozen_spec2.json` are the records of what was specified.
- `figs2.py` also adds the 95% range of the random 45-gene layer SDs to `canonical.json` (used for Fig. 1b).

## Independent check

`verify.py` re-derives the headline statistics with separately written code and can be run after step 5 to confirm them.

## Where each result comes from

| Manuscript item | Source |
|---|---|
| D, excess R², all D variants, bootstrap CI, jackknife, family exclusion | `canonical.json` |
| Sum-of-squares shares, REML components, layer SDs, pairwise correlations | `canonical.json` |
| Label-free module, enrichment, stability, sensitivity | `module.json` |
| Country decoding | `decode.json` |
| Sampling-effort adjustment | `effort.json` |
| Load statistics, ICC, socio-economic results, host breadth, coupling vs MAG agreement | `results.json` |
| Pre-specified intI1 and external network tests; sequencing-depth adjustment | `frozen_results.json` |
| Exploratory integron sensitivity analyses | `integron_depth.json` |
| Pre-specified external classification test (test C) | `frozen_results2.json` |
| External validation in untreated sewage (Section 3.9), primary | `frozen_results3_alr.json` |
| External validation, CLR sensitivity | `frozen_results3.json` |
| Drivers (test F) and campaign/sequencing (test H) | `frozen_results3_alr.json`, `f_specificity.json`, `test_H.json` |
| Figure 2d (per-gene intI1 correlations, pre-specified test B) | `intI1_rho.csv` |
| MAG representation and cell mixing | `controls.json` |
