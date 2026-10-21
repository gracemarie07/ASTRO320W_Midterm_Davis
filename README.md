# ASTRO 320W Midterm Observing Report - Grace Davis

Start with `report/main.pdf`. Editable LaTeX is in `report/main.tex`.
This folder contains the final paper and accompanying analysis code. The Artemis group collected the observations; Grace Davis performed the analysis and wrote the report. The final paper retains both the 6.86% Ogma scatter and the supplementary 1.33% unsaturated-star scatter, with the principal results boxed.

The LaTeX source and PDF in `report/` include the latest wording, individual analysis attribution, references, and boxed values. Upload the contents of this folder at the root of the GitHub repository. This package has not been pushed automatically.

## What changed

- Reconciled M29 numbering against the current AstroImageJ export. The report's original ten-star measurements were correct; the previous repository CSV used a different selection. The corrected export is `data/M29.csv`; the superseded selection is retained as `data/M29_superseded_selection.csv` and is not used in the report.
- Rebuilt the finding chart, color diagrams, formal covariance, and numerical notebook outputs from consistent star identifiers. `data/m29_reconciled.csv` records all twelve stars and inclusion status.
- Added an 8-by-8 detector comparison, formal dark-signal uncertainty estimates, raw focus and exposure-time peak checks, and the explicit joint paired-flat gain regression.
- Clarified the distinction between nonlinear stellar cores and the ADC ceiling, nominal dark signal and intrinsic dark current, and measured stellar width and atmospheric seeing.
- Added a status table for all ten binary targets. Four remain verified calibration pairs. HR6446 remains a candidate and is excluded; unconfirmed separations are not invented.
- Shortened and edited the report, grouped related figures, replaced the printed code listing with an analysis-file guide, and expanded the AI-assistance disclosure.

## Files and inputs

| Exercise | Main notebooks | Inputs |
|---|---|---|
| A | `A320_bf_stdev_final.ipynb`, `A320_bf_stdev_ROI.ipynb`, `A_dark_signal.ipynb` | Raw bias frames; recorded dark measurements |
| B | `B_gain_workbook.ipynb`, `B_photon_transfer_curves.ipynb` | Paired-flat workbook and CSV; original sky measurements; digitized dome PTC |
| C | `exercise_C_omicron_aquilae.ipynb`, `exercise_C_ogma.ipynb`, `C_unsaturated_field_star.ipynb` | Raw/calibrated stellar images and AIJ exports |
| D | `D_M29_colors.ipynb`, `D_NGC7129.ipynb` | Current `M29.csv`; raw nebula images and dome flats |
| E | `E_field_identifications.ipynb`, `E_supplementary_centroids.ipynb`, `E_plate_scale.ipynb`, `E_observed_fields.ipynb` | Double-star observations, bundled survey cutouts, verified centroid CSV |
| Revision | `Report_revision.ipynb`; `scripts/revise_analysis.py`, `scripts/report_panels.py`, `scripts/flat_profiles.py` | Bundled measurements plus optional external raw frames; panel script needs focus Image 6 |

The two A320 bias algorithms are instructor-supplied. Other notebook filenames reflect the original workflow, and some contain historical or supplementary analysis. `Results_summary.ipynb` now prints the current reconciled report results. The old excerpt `report/report_calculations_listing.txt` is a historical partial listing, not the revised appendix or an independent reproduction of every result.

## Installation

Use Python 3.10 or newer in a virtual environment:

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

For notebooks, launch Jupyter from the `notebooks/` directory so their relative project root resolves correctly:

```sh
cd notebooks
python3 -m jupyterlab
```

The scripts locate the project root relative to their own file, so they can be run from any directory.

## External raw-data configuration

Raw observing collections and the M29 reference FITS stacks are not bundled in this repository. Saved plots and CSV results are included, so the PDF and tabulated-result notebooks can be read without those images. Set paths before running image-dependent notebooks or scripts:

```sh
export ASTRO320W_DATA='/absolute/path/to/original_observing_collection'
export ASTRO320W_BIAS='/absolute/path/to/Bias_Frames'
export ASTRO320W_ORGANIZED='/absolute/path/to/organized/04_Data'
export ASTRO320W_E_IMAGES='/absolute/path/to/Exercise E/EXERCISE E'
```

`ASTRO320W_DATA` expects the original folder structure, including:

- `EXERCISE C-2/CCD Image 1.fit` through focus Images 2, 3, 4, 6, 7, 8, 9;
- `EXERCISE C-2/Ogma-0001.fit` through `Ogma-0015.fit`, and their existing calibration products for the supplementary notebook;
- `EXERCISE C-2/ogma_long_exposures/CCD Image 21.fit` through 26, with Images 15--20 under `ogma_long_exposures/iforgot/`;
- the original Exercise D dome/science folders and `Exercise E/EXERCISE E/` used by the legacy pipelines.

`ASTRO320W_BIAS` expects exactly twenty `biasframe-*.fit` images. `ASTRO320W_ORGANIZED` expects:

- `A/Dark_Time/` and `A/Dark_Temperature/`;
- `D/AIJ_Reference_Stacks/M29_B_stack_bias_corrected.fit`, and matching V and R files.

If not set, paths default to `raw/`, `raw/Bias_Frames/`, and `organized_data/` within this package. Not all legacy notebooks use every environment variable. Errors about missing files indicate which raw collection is needed. No raw frames are downloaded automatically.

## Reproducing the revisions

```sh
python3 scripts/revise_analysis.py
python3 scripts/report_panels.py
python3 scripts/flat_profiles.py
```

The first script always recomputes the reconciled M29 table, joint gain fit, and verified plate-scale result from bundled measurements. It additionally regenerates the finding chart, noise maps, formal dark errors, and raw stellar checks when their required external files exist. The flat-profile script requires the identified sky/dome files under `ASTRO320W_ORGANIZED/B/` and records their filenames and profile ranges in `results/flat_profile_provenance.json`. The second script composes report panels and requires the raw focus Image 6, plus the derived CSVs already included in this package.

`data/dark_*_with_errors.csv` contain conditional independent-pixel errors, not total physical dark-current errors. `data/m29_reconciled.csv` preserves formal AIJ errors; they are not claimed as validated total errors on averaged stacks. `data/ogma_raw_maxima.csv` distinguishes whole-image maxima from stellar-region peaks.

## Execution status and limitations

See `REVISION_VALIDATION.md` for what was actually recomputed. `notebooks/EXECUTION_historical.csv` is the original saved execution log and must not be interpreted as evidence that every notebook was rerun during this revision. Historical outputs are retained in other notebooks for reference.

Six binary targets still lack verified calibration values. The -20 C and -25 C files are not valid temperature-series exposures. The current stacks lack a verified variance model. These limitations are stated in the report; no missing observations or calibration records were fabricated.

## Building the report

A TeX Live installation with the `newtx` fonts and the packages in `report/main.tex` is required:

```sh
cd report
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The source enforces the manual's 25-page maximum. The final PDF was visually checked after compilation. ChatGPT assistance is disclosed in the report appendix and in the revision notes above.
