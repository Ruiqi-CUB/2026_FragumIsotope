# 2026_FragumIsotope

Public scripts for the Fraginae isotope manuscript.

## Contents

- `scripts/20260413_fragum_unedo_isotope_R.r`: fits the isotope models for *Fragum unedo*, runs basic diagnostics, and generates the size and height plots.
- `scripts/20260417_metadata_figure2.R`: builds the cardiid-focused isotope comparison plot used for Figure 2.
- `scripts/20260420_metadata_figure1.R`: builds the broad metadata isotope plot used for Figure 1, plus the labeled supplementary version.

## Expected inputs

The scripts expect the source data files to be available locally with the same filenames used in the code. If you keep the default paths, place them in a `data/` folder next to `scripts/`, or edit the file paths in the scripts.

## Outputs

Figures are written to a local `figures/` folder if it exists or can be created.
