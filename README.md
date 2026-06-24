# ASICS NMR Pipeline — Absolute Metabolite Quantification from CHO Spent Media

An automated R pipeline for **absolute metabolite quantification (mM)** from 800 MHz
Bruker NMR spectra of CHO cell culture spent media. It combines automated profiling
(via [ASICS](https://bioconductor.org/packages/ASICS/)) with a TSP-anchored manual
qHNMR panel, and produces QC reports and tidy outputs ready for downstream
visualisation (e.g. GraphPad).

---

## What it does

The pipeline runs end to end in a single script and produces two complementary
quantification layers:

1. **Manual qHNMR panel .** 30 metabolites quantified from
   resolved reporter peaks, anchored to the TSP internal standard (0.01% w/v).
   Library-independent and confidence-tiered (high / moderate / upper-bound).
   **This is the layer to report.**
2. **ASICS cross-calibrated table (broad screening).** The full ASICS library
   profile, cross-calibrated to mM against a set of manual-panel anchors.
   Wide coverage but absolute values carry scaling error — **use for relative,
   pathway-level interpretation, not for standalone absolute numbers.**


## Pipeline stages

| Section | Purpose |
|--------:|---------|
| 0–2 | Packages, configuration, helper functions |
| 3   | Load the built-in ASICS pure library (curated, ~193 metabolites) |
| 4   | Import Bruker spectra + **up-front dilution/metadata validation** + preprocess |
| 5   | ASICS quantification (full library searched) |
| 6   | Save results + reconstruction plots |
| 7   | TSP-anchored manual qHNMR panel → tube + dilution-corrected real mM |
| 8   | Cross-calibrated full ASICS table to mM (Route B) |
| 9   | QC / reproducibility (window plots, run manifest, `sessionInfo`) |
| 10  | Focused CHO target extraction |

## Requirements

- **R** ≥ 4.2 (developed on R 4.4)
- **ASICS** ≥ 2.28.0 — install from Bioconductor
- **CRAN packages:** `ggplot2`, `dplyr`, `tidyr`, `pheatmap`
- Bruker NMR spectra that have been **Fourier-transformed and processed** in TopSpin
  (each spectrum folder must contain a real-part `1r` file — raw FIDs alone will not work)

The script bootstraps its own dependencies on first run (`BiocManager::install("ASICS")`
and `install.packages()` for the CRAN packages), so a clean R install should work
without manual setup.

## Installation

```r
# from R / RStudio
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ASICS")
install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap"))
```

Then clone or download this repository and open the pipeline script.

## Configuration

Open `ASICS_NMR_pipeline_v25.R` and edit the **`SETTINGS YOU MUST CHECK`** block near
the top. The four paths/values you almost always need to change:

| Setting | What it is |
|---------|------------|
| `SPECTRA_DIR`   | Folder containing your processed Bruker spectra |
| `OUT_DIR`       | Where results are written (created if missing) |
| `DILUTION_FILE` | Path to your dilution CSV (see below) |
| `DEFAULT_DF`    | Fallback dilution factor for any unlisted sample |

An early `.check_inputs()`-style validator runs **before** the long quantification
step and reports common path/file problems in plain English, so you find mistakes
in seconds rather than after a multi-hour run.

## Input files

### Dilution file (required)

A CSV with **exactly** these two column names (lowercase):

```csv
sample,dilution_factor
10,1.1
13,5.5
```

- `sample` must be the **Bruker folder numbers** (10, 11, 12 …) — these have to match
  the spectrum folder names, or no samples will be matched.
- `dilution_factor` is the prep-specific factor (e.g. `1.1` for 500:50, `5.5` for 1:4.5).
- Any sample not listed inherits `DEFAULT_DF`.

A ready-to-edit template is provided: asics_dilution.csv

### Metadata file (optional)

Set `META_FILE` to a CSV with a `sample` column plus any of `condition`, `time_h`,
`biorep`. Column-name normalisation is forgiving (it auto-corrects a few common typos),
and missing optional columns are filled with sensible defaults.


## Running

Source the whole script (RStudio: *Source*, or `Rscript ASICS_NMR_pipeline_v25.R`).
Long runs write a live `progress_log.txt` into `OUT_DIR` so you can monitor progress.

## Outputs

All written to `OUT_DIR`. The key files:

| File | Contents |
|------|----------|
| `manual_panel_real_mM.csv` | **Report from this** — defensible absolute concentrations |
| `manual_panel_tube_mM.csv` | Manual panel before dilution correction |
| `absolute_ASICS_crosscalibrated_real_mM.csv` | Broad screening table (relative use) |
| `quantifications_all_samples.csv` | Raw ASICS quantification matrix |
| `CHO_target_metabolites_relative_intensity.csv` | Focused CHO target panel |
| `qc_panel_windows.pdf` | **Verify integration windows here** (colour-coded by tier) |
| `tsp_qc.csv` | TSP integral per sample + CV (internal-standard QC) |
| `reconstruction_sample*.pdf` | ASICS fit vs. original spectrum |
| `run_manifest.txt`, `sessionInfo.txt` | Full parameter + environment record |

## Default ASICS parameters

Tuned for clean 800 MHz baselines: `NOISE_THRES = 0.02`, `ADD_NOISE = 0.05`,
`MULT_NOISE = 0.172`, `MAX_SHIFT = 0.05`, `CLEAN_THRES = 10`, water exclusion
4.7–5.0 ppm, PQN normalisation, alignment on.

## Troubleshooting

Most first-run failures fall into one of these:

- **"DILUTION_FILE must have columns: sample, dilution_factor"** — the header is wrong.
  The names must match exactly and be lowercase. Watch for `Sample`, trailing spaces, or
  Excel renaming the column.
- **0 samples matched in the dilution file** — `sample` values aren't the Bruker folder
  numbers. They must equal the spectrum folder names.
- **Import errors / empty spectra** — the spectra haven't been processed. Each folder
  needs a real-part `1r` file; Fourier-transform and process in TopSpin first.
- **Windows path issues** — use forward slashes in paths (`C:/...`), not backslashes.
- **Slow or flaky I/O** — running directly off OneDrive can stall on file reads.
  Point `SPECTRA_DIR`/`OUT_DIR` at a local (synced-then-paused or non-cloud) folder
  for the duration of a run.

## Citation

If you use this pipeline, please cite the ASICS package alongside this repository:

> Tardivel P. *et al.* ASICS: an R package for a whole analysis workflow of 1D ¹H NMR
> spectra. *Bioinformatics* (2017).

A `CITATION.cff` can be added if you'd like GitHub to surface a "Cite this repository"
button.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).

ASICS and its dependencies are installed from Bioconductor/CRAN under their own
licenses and are not redistributed here.

## Acknowledgements

Developed for NMR metabolomics of CHO cell spent media at the University of Manchester.
Built on [ASICS](https://bioconductor.org/packages/ASICS/).
