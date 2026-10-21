# Revision validation

## Recomputed

- `scripts/revise_analysis.py`: current M29 measurements and covariance; joint twenty-pair gain regression; four-binary plate scale and empirical SEM; all 64 bias-map blocks from twenty raw bias frames; formal regional dark uncertainties; raw focus and Ogma peak checks.
- `scripts/report_panels.py`: compact dark, PTC, focus, repeatability, exposure-time, and color-diagram panels; final identification figures.
- `scripts/flat_profiles.py`: identified sky/dome flat profiles, with normalized central 11-pixel-wide median strips, and matching representative images. These supersede the untraceable original profile plots in the revised report.
- Eight notebooks' code cells were executed sequentially in Python: `A_dark_signal`, `B_gain_workbook`, `B_photon_transfer_curves`, `C_unsaturated_field_star`, `D_M29_colors`, `E_plate_scale`, `E_supplementary_centroids`, and `Results_summary`. Only the inline display magic was skipped. Standard-output results and execution counts were saved. This was direct Python execution, not a Jupyter-kernel run. All eight passed. Details are in `results/revision_execution.json`.
- The main reported results were reproduced: gain 1.2901515 with slope-propagated error 0.0062917 electrons/DN; plate scale 1.2513698 with between-binary SEM 0.0115021 arcseconds/pixel; unsaturated-star scatter 1.33125 percent.
- The current M29 export reproduces the original report's Tables 10--11 (renumbered in the revision). The older repository selection was preserved but is no longer the notebook input.

## PDF checks

The report was compiled with pdfLaTeX/latexmk. Cross-references were resolved, the page count checked against the 25-page limit, and rendered pages inspected for overlaps, clipping, figure legibility, and appendix layout. The final edition is supplied separately as a PDF and in `report/main.pdf`.

## What was not established

- The original instructor pixel-loop bias notebooks, full NGC 7129 pipeline, and all historical image-alignment or survey-identification notebooks were not rerun in full. The raw bias maps independently use the available twenty-frame stack. Previously saved outputs in other notebooks remain historical.
- The conditional dark uncertainties exclude temperature mismatch, stray light, spatial correlations, and gain systematics. They are not total errors on physical dark current.
- M29's formal AIJ errors are preserved and qualified; a validated variance model for the averaged/calibrated stacks is unavailable. Dividing all quoted errors by sqrt(15) would not account for shared calibration or resampling correlations.
- Four binaries remain verified. HR6446 is a candidate; five others lack secure usable separations. The status table describes the observations without treating unknown separations as measured calibration results.
- Files named -20C and -25C have incompatible exposures and no recorded CCD temperature; they do not complete the temperature series.
- This is a local revision package, not a published GitHub update. Read the edited prose and assistance disclosure before submission to confirm it accurately describes your work.
