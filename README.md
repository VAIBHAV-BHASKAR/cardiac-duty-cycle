# Cardiac Duty Cycle — Python Reproduction

**Independent Python reproduction** of all analyses and figures from:

> Froese, T. et al. (in prep.). The healthy heart converges on a thermodynamic optimum (1/*e*) that predicts survival. *Nature Aging*.

This repository verifies every statistical result and manuscript figure using Python (NumPy, SciPy, pandas, Matplotlib), independent of the original MATLAB implementation.

[![Python](https://img.shields.io/badge/Python-3.10+-blue)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Zenodo Data](https://img.shields.io/badge/Data-Zenodo-blue)](https://doi.org/10.5281/zenodo.19246123)

---

## How to reproduce

1. Clone this repository:

   ```bash
   git clone https://github.com/VAIBHAV-BHASKAR/cardiac-duty-cycle.git
   cd cardiac-duty-cycle
   ```

2. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

3. Run the reproduction:

   ```bash
   python reproduce_paper.py
   ```

   On first run, the script automatically downloads the large CSV files (~460 MB) from Zenodo. Three small CSVs (LUDB, QTDB, PTB) are included in the repository.

All results appear in `results/` and figures in `results/figures/` (PNG + PDF).

---

## What this reproduces

### Statistical analyses (7)

| # | Analysis | Key result |
|---|----------|------------|
| 1 | Gold standard (LUDB + QTDB) | Healthy CDC mode near 1/*e*; pathological shifted upward |
| 2 | Large-scale (Fantasia, AA, PTB, PTB-XL) | CDC stable with age in healthy controls |
| 3 | Hierarchical model (all databases pooled) | Age x Group interaction: healthy flat, pathological rising |
| 4 | CDC vs heart rate | CDC predicts resting HR = ~64 bpm at 1/*e* |
| 5 | CODE-15% standalone | Normal vs pathological CDC distributions diverge |
| 6 | Age-stratified mortality | T3/T1 mortality fold ~2x at every age decade |
| 7 | Survival curves | Kaplan-Meier separation by CDC-deviation tertile |

### Figures (9)

| Figure | Description |
|--------|-------------|
| Fig 1 | CDC aging trajectories (healthy vs pathological) |
| Fig 2 | Age-stratified mortality by CDC-deviation tertile |
| SI Fig 1 | CDC vs heart rate + hidden mortality risk |
| SI Fig 2 | Systole & RR interval across the lifespan |
| SI Fig 4 | Gold-standard distributions (LUDB + QTDB) |
| SI Fig 5 | Large-scale distributions (4 databases) |
| SI Fig 6 | CODE-15% distributions + age trends |
| SI Fig 7 | Kaplan-Meier survival curves (overall, F, M) |
| SI Fig 8 | Sex-stratified CDC aging |

---

## Folder structure

```
cardiac-duty-cycle/
├── reproduce_paper.py          ← one-command reproduction
├── code/
│   ├── config.py               ← paths, constants, colour palette
│   ├── data_loader.py          ← quality filters, aggregation, bootstrap mode
│   ├── analyses/               ← 7 analysis modules
│   │   ├── gold_standard.py
│   │   ├── large_scale.py
│   │   ├── hierarchical.py
│   │   ├── cdc_vs_hr.py
│   │   ├── code15.py
│   │   ├── mortality.py
│   │   └── survival.py
│   └── figures/                ← 9 figure modules
│       ├── fig1_aging.py
│       ├── fig2_mortality.py
│       ├── si_fig1_hr.py
│       ├── si_fig2_systole.py
│       ├── si_fig4_gold.py
│       ├── si_fig5_large.py
│       ├── si_fig6_code15.py
│       ├── si_fig7_km.py
│       └── si_fig8_sex.py
├── data/
│   └── preprocessed/           ← beat-level CSVs (small ones committed, large from Zenodo)
├── results/                    ← generated outputs (gitignored)
├── requirements.txt
└── LICENSE
```

---

## Data

All analyses use the preprocessed beat-level CSVs from:

> Froese, T. (2026). CDC Analysis — Preprocessed Beat Data (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19246123

Seven databases, ~250,000 subjects, unified 14-column format (one row per beat). Quality filters are applied at analysis time, not at export — matching the MATLAB pipeline exactly.

---

## Relationship to the MATLAB implementation

This is a **fully independent reimplementation** in Python. It was developed separately from the original MATLAB toolbox ([tom-froese/cdc-ecg-analysis](https://github.com/tom-froese/cdc-ecg-analysis)) to provide cross-language verification of all results.

The same preprocessed data, the same quality filters, and the same statistical tests — different language, different researcher.

---

## Citation

If you use this code, please cite both the paper and the data:

> Froese, T. et al. (in prep.). The healthy heart converges on a thermodynamic optimum (1/*e*) that predicts survival. *Nature Aging*.

> Froese, T. (2026). CDC Analysis — Preprocessed Beat Data (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19246123

---

## License

[MIT](LICENSE) — Vaibhav Bhaskar, 2026

---

## Author

Vaibhav Bhaskar — KIIT University
