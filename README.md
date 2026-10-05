# Analysis code: wastewater resistome module

This archive contains the analysis code and the result files it produces, together with the
released gene panel and the scoring tool. Figure-drawing and document-preparation scripts are
not included; everything here either computes a reported result or is a result.

## Input data

None of the input data are redistributed here. All are published and must be downloaded from
their original sources, whose licences apply.

| Source | Files used |
|---|---|
| Zhu et al. (2025) *Nature Communications* 16, 4006 | Source Data workbook (`41467_2025_59019_MOESM8_ESM.xlsx`); Supplementary Data 1 (per-sample sequencing output); Supplementary Data 3 (mobile genetic element abundances) |
| Martiny et al. (2025) *Nature Communications* 16, 10278 | Supplementary Data 1 (sample metadata); Supplementary Data 2 (gene annotation); Supplementary Data 3 and 4 (co-occurrence network nodes and links) |
| Zenodo record 14652833 | `panres_counts.csv` (gene-level counts) and `kingdom_motus_pad_agg.csv` (bacterial marker-gene counts) |
| NCBI BioProject PRJNA509305 | run table (`SraRunTable.csv`) |
| WHO GLASS, distributed by Our World in Data | `antibiotic-consumption-rate.csv` |
| WHO GLASS 2022 report, compiled by Leclerc (GitHub `qleclerc/GLASS2022`) | `compiled_WHO_GLASS_2022.csv` |
| World Bank World Development Indicators | GDP per capita, PPP (`NY.GDP.PCAP.PP.KD`) |

Each script reads its inputs from a directory set at the top of the file. Edit that path once
per script, or place the downloaded files in the expected directory.

## Environment

Python 3.12 with numpy, pandas, scipy, scikit-learn, statsmodels and openpyxl. Random seeds are
fixed throughout, so every reported statistic is reproducible.

## Frozen specifications

Five analysis specifications were written out in full and hashed (SHA-256) before the
corresponding test scripts were run. Each `freeze_spec*.py` writes a specification file; the
matching `run_frozen*.py` recomputes the hash and refuses to run if the specification has
changed. The hashes are recorded inside the specification files.

Note on interpretation: hashing fixes an analysis plan after the discovery data have been
inspected. It is specification freezing, not prospective registration of the underlying
hypothesis.

## Run order

Each step reads only files written by earlier steps.

| Step | Script | Writes | Computes |
|---|---|---|---|
| 1 | `analysis.py` | `results.json`, `gene_table.csv`, `sample_layers.csv`, `city_table.csv` | Sample counts, load statistics, reliability, socio-economic models, host breadth, mobility indices |
| 2 | `make_gene_table2.py` | `gene_table2.csv` | Gene abundance, variability and detection columns |
| 3 | `controls.py` | `controls.json` | Genome-bin representation and class-by-mechanism mixing |
| 4 | `canonical.py` | `canonical.json` | Primary statistics: the between-country contrast with 9,999 permutations, controls, bootstrap interval, jackknife, family exclusion, variance shares, layer variability, REML |
| 5 | `arrays.py` | `gene_table_final.csv`, `Dnull.npy`, `rhonull.npy` | Per-gene results and permutation nulls |
| 6 | `rf_fold.py` | adds `transfer_R2` to `canonical.json` | Random forests with imputation fitted inside folds |
| 7 | `decode_final.py` | `decode.json`, `rand_acc.npy` | Country decoding from held-out cities |
| 8 | `module.py` | `module.json`, `module_membership.csv` | Label-free module discovery; the design is fixed in the file header before running |
| 9 | `effort.py` | `effort.json` | Sampling-effort adjustment |
| 10 | `freeze_spec.py`, `run_frozen.py` | `frozen_spec.json`, `frozen_results.json` | Network co-occurrence and integron-integrase tests |
| 11 | `intI1_pergene.py` | `intI1_rho.csv` | Per-gene integron-integrase correlations |
| 12 | `freeze_spec2.py`, `run_frozen2.py` | `frozen_spec2.json`, `frozen_results2.json` | External acquired-gene classification test |
| 13 | `build_martiny.py` | `martiny_matrix.pkl` | Sample-by-gene matrix for the external cohort |
| 14 | `make_bacterial_reference.py` | `bacterial_reference.csv` | Bacterial marker-gene totals used for normalisation |
| 15 | `freeze_spec3.py`, `run_frozen3.py`, `finish_frozen3.py`, `run_frozen3_alr.py` | `frozen_spec3.json`, `frozen_results3.json`, `frozen_results3_alr.json` | External validation, temporal stability, exposure association, cross-compartment comparison and batch tests |
| 16 | `f_specificity.py` | `f_specificity.json` | Whether the exposure association is specific to the module |
| 17 | `test_H.py` | `test_H.json` | Run-level and collection-date adjustments |
| 18 | `integron_depth.py` | `integron_depth.json` | Integron-integrase abundance against sequencing depth |
| 19 | `freeze_spec4.py`, `run_frozen4.py`, `clinical_added.py` | `frozen_spec4.json`, `frozen_results4.json`, `clinical_added.json` | Clinical resistance tests and income-adjusted models |
| 20 | `freeze_spec5.py`, `run_frozen5.py`, `class_placebo.py` | `frozen_spec5.json`, `frozen_results5.json`, `class_placebo.json` | Within-country panel and drug-class-matched tests, and the class-by-phenotype grid |
| 21 | `crossfit.py` | `crossfit.json` | Cross-fitted module discovery and country-balanced discovery |
| 22 | `benchmark.py` | `benchmark.json` | Benchmark of five candidate metrics |
| 23 | `bootstrap_compare.py` | `bootstrap_compare.json` | Bootstrapped differences between metrics and the interval for the exclusion analysis |
| 24 | `ml_benchmark.py city`, `ml_benchmark.py country` | `ml_benchmark_city.json`, `ml_benchmark_country.json` | Prediction of each metric from plant metadata under grouped cross-validation |
| 25 | `specificity.py` | `specificity.json` | Comparison against size-, abundance-, prevalence- and class-matched gene sets |
| 26 | `robustness_metric.py` | `robustness_metric.json` | Sensitivity to sequencing depth, normalisation and annotation scheme |
| 27 | `minimal_panel.py` | `minimal_panel.json` | Reduced panels and their transfer to the external cohort |
| 28 | `mdd.py` | `mdd.json` | Smallest detectable difference between plants |
| 29 | `design.py aliased`, `design.py crossed`, `design2.py grid`, `design2.py season` | `design_aliased.json`, `design_crossed.json`, `design2_grid.json`, `design2_season.json` | Survey design simulations, including seasonal sensitivity |
| 30 | `make_panel.py` | `resistome_module_panel.csv` | The released gene panel with cross-database identifiers |

Five intermediate files written by steps 1–8 (`results.json`, `canonical.json`, `module.json`,
`gene_table_final.csv`, `module_membership.csv`) are not included in this archive and are
regenerated by running those steps.

## Scoring tool

`module_score.py` computes the module score, the module share and the mobile share from any ARG
abundance table.

```
python module_score.py --table my_abundances.csv --out scores.csv
```

The table must have samples as rows and gene names as columns. Gene names are matched to the
panel through the normalised names and cross-database identifiers in
`resistome_module_panel.csv`. Genes absent from the input are omitted rather than imputed, so
the score is a mean over matched genes; the tool reports how many matched and warns below ten,
where agreement with the full-panel score becomes unreliable.

## Result files

Every value reported in the manuscript is traceable to one of the JSON files listed above. The
specification files carry the recorded hashes.
