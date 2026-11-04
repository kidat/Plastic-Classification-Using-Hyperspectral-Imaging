# Plastic Classification from SWIR Hyperspectral Images

This project classifies unknown plastic samples from SWIR hyperspectral data
using reference transmittance spectra measured in the Chemistry Lab.

## Method

1. Apply white/dark reference correction.
2. Estimate single-pass transmittance from the white-background measurement.
3. Extract spectra from manually selected unknown-plastic ROIs.
4. Compare normalized ROI spectra with known Chemistry Lab references.
5. Create an optimized false-color SWIR image using bands with high ROI
   separation.

## Results

| Selected unknown-plastic ROIs | ROI transmittance spectra |
|---|---|
| ![Selected ROIs](Result/minimal_roi_selection.png) | ![ROI spectra](Result/minimal_roi_spectra.png) |
| The scanner image with the unknown-plastic ROIs marked for analysis. | Single-pass transmittance spectra extracted from the selected ROIs. |

| Optimized SWIR visualization | Unknown-to-reference classification |
|---|---|
| ![Optimized SWIR image](Result/minimal_rgb_optimized.png) | ![Reference comparison](Result/minimal_roi_reference_comparison.png) |
| False-color SWIR composite using bands selected for strong ROI variation. | Each unknown ROI compared with its closest known Chemistry Lab plastic reference. |

## Run

Install the required packages:

```bash
pip install matplotlib numpy pandas spectral
```

Run:

```bash
python Plastic_classificationl.py
```

The script writes result images to `Result/`.

## Data

The script expects the SWIR scan and reference CSV files outside this folder:

- `Group4/Scanner/ir plastics_scan_0552/capture/`
- `Transmittance/`

Update the paths in the script if the data are stored elsewhere. Raw scanner
data are not included in this repository.
