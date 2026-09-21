Mobile antibiotic resistance genes form a reproducible, integron-associated module in wastewater

Supplementary computational code and frozen results for the manuscript:

“Mobile antibiotic resistance genes form a reproducible, integron-associated module in wastewater”**

This repository contains the computational analyses, frozen specifications, verification procedures, external validation code, and numerical outputs underlying the manuscript.

The repository is intended to provide a transparent record of the analyses and to allow the reported computational results to be independently checked.

## Repository scope

This repository contains:

* core analysis and data-processing code;
* ARG-module construction and analysis;
* random-forest cross-validation;
* country decoding;
* sampling-effort analysis;
* integron-associated analyses;
* control analyses;
* pre-specified frozen tests;
* external validation analyses;
* independent verification code;
* scripts for constructing numerical tables;
* frozen specifications and frozen numerical results.
---

Repository structure

```text
mobile-arg-module-wastewater/
│
├── external_validation/
│   ├── build_martiny.py
│   ├── f_specificity.py
│   ├── make_bacterial_reference.py
│   └── test_H.py
│
├── results/
│   ├── frozen_results.json
│   ├── frozen_results2.json
│   ├── frozen_results3.json
│   └── frozen_results3_alr.json
│
├── src/
│   ├── analysis.py
│   ├── arrays.py
│   ├── canonical.py
│   ├── controls.py
│   ├── decode_final.py
│   ├── effort.py
│   ├── integron_depth.py
│   ├── module.py
│   ├── rf_fold.py
│   │
│   ├── frozen tests/
│   │   ├── finish_frozen3.py
│   │   ├── freeze_spec.py
│   │   ├── freeze_spec2.py
│   │   ├── freeze_spec3.py
│   │   ├── run_frozen.py
│   │   ├── run_frozen2.py
│   │   ├── run_frozen3.py
│   │   └── run_frozen3_alr.py
│   │
│   └── tables/
│       ├── build_tables.py
│       ├── make_gene_table2.py
│       └── make_xlsx.py
│
├── results/
│   ├── frozen_results.json
│   ├── frozen_results2.json
│   ├── frozen_results3.json
│   └── frozen_results3_alr.json
│
└── verify.py
```

---

## Input data

The analyses use external datasets that are not redistributed in this repository.

Wastewater ARG and metagenomic data

Additional inputs used for the analyses include:

* `panres_counts.csv`;
* `kingdom_motus_pad_agg.csv`;
* the NCBI run table for PRJNA509305 (`SraRunTable.csv`);
* national antibiotic-consumption data (`antibiotic-consumption-rate.csv`).

The relevant source data should be obtained from their original repositories or publishers, subject to their respective access conditions and licences.

Zhu et al. (2025)

The analyses use data associated with:

**Zhu et al. (2025), Nature Communications 16, 4006.**

Relevant materials include:

* Source Data workbook: `41467_2025_59019_MOESM8_ESM.xlsx`;
* Supplementary Data 1 (`MOESM3`);
* Supplementary Data 3 (`MOESM5`).

These files are not redistributed here.

Martiny et al. (2025)

External validation uses data associated with:

Martiny et al. (2025), Nature Communications 16, 10278.**

Relevant materials include:

* Supplementary Data 3 (`41467_2025_66070_MOESM6_ESM.xlsx`);
* Supplementary Data 4 (`MOESM7`).

These files are not redistributed here.

The scripts that require external inputs should be configured with the local paths to the downloaded source files.

---

## Core analysis

The principal computational implementation is contained in `src/`.

| Script              | Function                                                                                           |
| ------------------- | -------------------------------------------------------------------------------------------------- |
| `analysis.py`       | Main sample, abundance, load, ICC, socio-economic, host-breadth and mobility analyses              |
| `arrays.py`         | Construction of final gene-level arrays and permutation null distributions                         |
| `canonical.py`      | Primary statistical analyses, including D statistics, controls, bootstrap and jackknife procedures |
| `controls.py`       | Control analyses for MAG representation and class × mechanism mixing                               |
| `decode_final.py`   | Country-decoding analysis                                                                          |
| `effort.py`         | Sampling-effort adjustment                                                                         |
| `integron_depth.py` | Integron-associated sensitivity analyses                                                           |
| `module.py`         | Label-free module discovery and module membership analysis                                         |
| `rf_fold.py`        | Random-forest cross-validation with imputation performed within folds                              |

---

Frozen specifications and tests

Three frozen specifications are included:

```text
src/frozen tests/freeze_spec.py
src/frozen tests/freeze_spec2.py
src/frozen tests/freeze_spec3.py
```

with their corresponding frozen specification records:

```text
frozen_spec.json
frozen_spec2.json
frozen_spec3.json
```

The frozen specifications preserve the pre-specified computational tests used in the analysis.

The corresponding execution scripts are:

```text
run_frozen.py
run_frozen2.py
run_frozen3.py
run_frozen3_alr.py
```

`finish_frozen3.py` completes the third frozen analysis where required.

The frozen specifications should be treated as fixed records of the analyses that were specified before evaluation of the corresponding results.

---

Frozen numerical results

The repository contains four frozen result files:

```text
results/
├── frozen_results.json
├── frozen_results2.json
├── frozen_results3.json
└── frozen_results3_alr.json
```

These files preserve the numerical outputs of the frozen analyses.

Broadly:

| File                       | Analysis                                                                       |
| -------------------------- | ------------------------------------------------------------------------------ |
| `frozen_results.json`      | Pre-specified integron-associated and external-network tests                   |
| `frozen_results2.json`     | Pre-specified external classification test                                     |
| `frozen_results3.json`     | External validation and associated sensitivity analyses using the CLR workflow |
| `frozen_results3_alr.json` | External validation using the pre-specified ALR normalisation                  |

The frozen JSON files are retained as reference outputs for computational verification and reproducibility.

---

External validation

The `external_validation/` directory contains analyses supporting the external validation component.

`build_martiny.py`

Constructs the sample-by-gene matrix required for the external validation using the Martiny et al. data.

### `make_bacterial_reference.py`

Constructs the bacterial reference used for the pre-specified ALR normalisation.

### `f_specificity.py`

Performs the specificity analysis associated with the external validation workflow.

### `test_H.py`

Runs the run-level adjustment analysis and associated collection-date analyses.

---

Numerical tables

The `src/tables/` directory contains scripts for constructing numerical/tabular outputs:

```text
build_tables.py
make_gene_table2.py
make_xlsx.py
```

These scripts are included to document the generation of the numerical tables associated with the computational analyses.

They are not figure-generation scripts.

---

Independent verification

The top-level script:

```text
verify.py
```

provides an independently written verification of the principal reported statistics.

It is intended as an additional check on the computational implementation and does not replace the primary analysis scripts.

---

Reproducibility

The computational workflow uses fixed random seeds where stochastic procedures are involved.

Reproducibility requires:

1. the Python environment specified below;
2. the external input datasets listed above;
3. the source code in this repository;
4. the frozen specification files;
5. the frozen result files.

Because external datasets are not redistributed here, users must obtain those datasets from their original sources and configure the corresponding local file paths in the scripts where necessary.

The frozen JSON files provide fixed reference outputs against which the implementation can be checked.

---

## Python environment

The analyses were developed and tested using **Python 3.12**.

The principal Python dependencies are:

```text
numpy
pandas
scipy
scikit-learn
statsmodels
openpyxl
```

A `requirements.txt` file should be used to record the tested package environment.

If package versions are important for exact reproduction, the versions used for the final analysis should be pinned in `requirements.txt`.

---
Relationship to the manuscript

The repository contains the computational work underlying the analyses reported in the manuscript, including:

* ARG abundance and variability analyses;
* module discovery;
* integron-associated analyses;
* controls and sensitivity analyses;
* random-forest transfer analysis;
* country decoding;
* sampling-effort adjustment;
* frozen pre-specified tests;
* external validation;
* numerical table construction.

---

Data and code availability

Code developed specifically for the computational analyses is provided in this repository.

External datasets are not redistributed where they are already publicly available from their original repositories or publishers. Users should obtain the source data directly from those providers and comply with the applicable data-access conditions and licences.

The frozen specification and result files are included to preserve the computational record of the reported analyses.

---

Citation

If you use this code or the associated computational workflow, please cite the accompanying manuscript:

> Hassan et al. *Mobile antibiotic resistance genes form a reproducible, integron-associated module in wastewater.*

Please use the final published citation once available.

---

License

Unless otherwise specified by the repository owner, the code in this repository is provided under the **MIT License**.

External datasets referenced by the code are **not covered by this software licence** and remain subject to their original source-specific licences and terms of use.
