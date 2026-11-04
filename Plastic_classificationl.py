"""Read SWIR plastic spectra and apply white-dark reflectance correction."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from spectral import open_image


# ============================================================================
# LOAD HYPERSPECTRAL DATA
# ============================================================================
base_path = (
    Path(__file__).resolve().parents[3]
    / "Group4"
    / "Scanner"
    / "ir plastics_scan_0552"
    / "capture"
)
scan_name = "ir plastics_scan_0552"
output_path = Path(__file__).resolve().parent / "Result"
scanner_image_path = base_path.parent / f"{scan_name}.png"
output_path.mkdir(parents=True, exist_ok=True)

sample_path = base_path / f"{scan_name}.hdr"
white_path = base_path / f"WHITEREF_{scan_name}.hdr"
dark_path = base_path / f"DARKREF_{scan_name}.hdr"

sample_img = open_image(sample_path)
wavelength = np.asarray(sample_img.bands.centers)

sample_cube = sample_img.load().astype(np.float32)
white_cube = open_image(white_path).load().astype(np.float32)
dark_cube = open_image(dark_path).load().astype(np.float32)

lines, samples, bands = sample_cube.shape
print(f"Data dimensions - Lines: {lines}, Samples: {samples}, Bands: {bands}")

if white_cube.shape[1:] != (samples, bands) or dark_cube.shape[1:] != (samples, bands):
    raise ValueError(
        "White and dark reference cubes must have the same samples and bands "
        "as the sample cube."
    )


# ============================================================================
# REFLECTANCE CALIBRATION
# ============================================================================
# Average the references; NumPy broadcasts them across all image lines.
white_ref = white_cube.mean(axis=0)
dark_ref = dark_cube.mean(axis=0)

denominator = white_ref - dark_ref
with np.errstate(divide="ignore", invalid="ignore"):
    reflectance = np.divide(
        sample_cube - dark_ref,
        denominator,
        out=np.full_like(sample_cube, np.nan),
        where=np.abs(denominator) > 1e-6,
    )

# The white-background measurement contains two passes through the plastic,
# so the corrected ratio is approximately T^2. Recover single-pass
# transmittance with the square root.
transmittance = np.sqrt(np.maximum(reflectance, 0.0))
transmittance2 = 2*reflectance

# ============================================================================
# ROI SPECTRA
# ============================================================================
# ROI coordinates use half-open slices: (x1, x2, y1, y2).
rois = {
    "P-01": (33, 41, 45, 53),
    "P-02": (205, 213, 73, 81),
    "P-05_1": (50, 58, 170, 178),
    "P-05_2": (206, 214, 105, 113),
    "P-06_1": (93, 101, 227, 235),
    "P-06_2": (190, 198, 195, 203),
    "P-07_1": (204, 212, 35, 43),
    "P-07_2": (145, 153, 32, 40),
    "P-07_3": (199, 207, 149, 157),
    "P-bg": (161, 169, 175, 183),
}

spectra = {}
for name, (x1, x2, y1, y2) in rois.items():
    if not (0 <= x1 < x2 <= samples and 0 <= y1 < y2 <= lines):
        raise ValueError(f"ROI {name} is outside the image bounds.")
    spectra[name] = np.nanmedian(transmittance[y1:y2, x1:x2, :], axis=(0, 1))

# Compare each ROI with the known plastic reference spectra.
reference_path = Path(__file__).resolve().parent.parent / "Transmittance"
reference_files = {
    "HD-PE": "HD-PE-transmittance.Sample.Raw.csv",
    "PC": "PC-transmittance.Sample.Raw.csv",
    "POM": "POM-transmittance.Sample.Raw.csv",
    "PP(H)": "PP(H)-transmittance.Sample.Raw.csv",
    "PS": "PS-transmittance.Sample.Raw.csv",
}

def normalize(spectrum):
    spectrum = np.asarray(spectrum, dtype=float)
    if not np.all(np.isfinite(spectrum)):
        raise ValueError("Spectrum contains non-finite values.")
    low, high = np.min(spectrum), np.max(spectrum)
    return (spectrum - low) / (high - low) if high > low else np.zeros_like(spectrum)

reference_spectra = {}
for plastic, filename in reference_files.items():
    reference = pd.read_csv(reference_path / filename).sort_values("nm")
    reference_wavelength = reference["nm"].to_numpy(dtype=float)
    reference_transmittance = reference["T"].to_numpy(dtype=float) / 100
    reference_spectra[plastic] = (reference_wavelength, reference_transmittance)

scanner_min = float(wavelength.min())
scanner_max = float(wavelength.max())
ref_min = max(rn.min() for rn, _ in reference_spectra.values())
ref_max = min(rn.max() for rn, _ in reference_spectra.values())
w_lo = max(scanner_min, ref_min)
w_hi = min(scanner_max, ref_max)

if w_hi <= w_lo:
    raise ValueError(
        f"No wavelength overlap: scanner {scanner_min:.0f}-{scanner_max:.0f} nm, "
        f"references {ref_min:.0f}-{ref_max:.0f} nm."
    )

shared_bands = (wavelength >= w_lo) & (wavelength <= w_hi)
comparison_wavelength = wavelength[shared_bands]
print(
    f"Comparison range: {w_lo:.0f}-{w_hi:.0f} nm "
    f"({comparison_wavelength.size} bands)"
)

classification_lines = ["ROI classification"]
classifications = {}
for name, spectrum in spectra.items():
    if name == "P-bg":
        continue
    roi_spectrum = normalize(spectrum[shared_bands])
    distances = {
        plastic: np.linalg.norm(
            roi_spectrum
            - normalize(np.interp(comparison_wavelength, ref_nm, ref_spectrum))
        )
        for plastic, (ref_nm, ref_spectrum) in reference_spectra.items()
    }
    plastic = min(distances, key=distances.get)
    classifications[name] = plastic
    classification_lines.append(f"{name}: {plastic} ({distances[plastic]:.4f})")

# ============================================================================
# SAVE ROI IMAGE AND SPECTRA
# ============================================================================
# Select bands that separate the plastic ROIs most strongly for the RGB image.
rgb_spectra = {
    name: np.nanmedian(transmittance2[y1:y2, x1:x2, :], axis=(0, 1))
    for name, (x1, x2, y1, y2) in rois.items()
    if name != "P-bg"
}
plastic_spectra = np.array(list(rgb_spectra.values()))
band_score = np.nanstd(plastic_spectra, axis=0)
ranked_bands = np.argsort(band_score)[::-1]
rgb_bands = []
for band in ranked_bands:
    if all(abs(wavelength[band] - wavelength[selected]) >= 100 for selected in rgb_bands):
        rgb_bands.append(band)
    if len(rgb_bands) == 3:
        break
rgb_bands.sort()
rgb_wavelengths = wavelength[rgb_bands]

rgb = transmittance2[:, :, rgb_bands].copy()
for channel in range(3):
    low, high = np.nanpercentile(rgb[:, :, channel], [1, 99])
    if high > low:
        rgb[:, :, channel] = np.clip((rgb[:, :, channel] - low) / (high - low), 0, 1)

figure, axis = plt.subplots(figsize=(9, 6))
axis.imshow(rgb)
axis.axis("off")
figure.tight_layout()
figure.savefig(output_path / "minimal_rgb_optimized.png", dpi=200)
plt.close(figure)

# Show the same ROIs on the scanner-provided image.
roi_image = plt.imread(scanner_image_path)
figure, axis = plt.subplots(figsize=(9, 6))
axis.imshow(roi_image)
for index, (name, (x1, x2, y1, y2)) in enumerate(rois.items()):
    color = plt.cm.tab10(index % 10)
    axis.add_patch(
        plt.Rectangle(
            (x1, y1),
            x2 - x1,
            y2 - y1,
            fill=False,
            edgecolor=color,
            linewidth=1.5,
        )
    )
    axis.text(x1, y1 - 2, name, color=color, fontsize=8)
axis.axis("off")
figure.tight_layout()
figure.savefig(output_path / "minimal_roi_selection.png", dpi=200)
plt.close(figure)

figure, axis = plt.subplots(figsize=(9, 5))
for name, spectrum in spectra.items():
    axis.plot(wavelength, spectrum, label=name)

axis.set_xlabel("Wavelength (nm)")
axis.set_ylabel("Single-pass transmittance")
axis.grid(alpha=0.3)
axis.legend(fontsize=8, ncol=2)
figure.tight_layout()
figure.savefig(output_path / "minimal_roi_spectra.png", dpi=200)
plt.close(figure)

figure, axes = plt.subplots(3, 3, figsize=(14, 11), sharex=True, sharey=True)
for axis, (name, plastic) in zip(axes.flat, classifications.items()):
    roi_spectrum = normalize(spectra[name][shared_bands])
    ref_nm, ref_spectrum = reference_spectra[plastic]
    reference_spectrum = normalize(
        np.interp(comparison_wavelength, ref_nm, ref_spectrum)
    )
    axis.plot(comparison_wavelength, roi_spectrum, label=f"{name} measured")
    axis.plot(
        comparison_wavelength,
        reference_spectrum,
        "--",
        label=f"{plastic} reference",
    )
    axis.set_title(f"{name} -> {plastic}")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8)

for axis in axes[-1, :]:
    axis.set_xlabel("Wavelength (nm)")
for axis in axes[:, 0]:
    axis.set_ylabel("Normalized transmittance")
figure.suptitle("Classification of Unknown Plastics Using Chemistry Lab References")
figure.tight_layout()
figure.savefig(output_path / "minimal_roi_reference_comparison.png", dpi=200)
plt.close(figure)

figure, axis = plt.subplots(figsize=(12, 7))
for plastic, (ref_nm, ref_spectrum) in reference_spectra.items():
    axis.plot(
        ref_nm,
        normalize(ref_spectrum),
        linewidth=2,
        label=f"{plastic} reference",
    )

axis.set_xlabel("Wavelength (nm)")
axis.set_ylabel("Normalized transmittance")
axis.set_title("Known Chemistry Lab Plastic Reference Spectra")
axis.grid(alpha=0.3)
axis.legend(fontsize=9)
figure.tight_layout()
figure.savefig(output_path / "minimal_reference_spectra.png", dpi=200)
plt.close(figure)

print(f"Wavelength range: {wavelength[0]:.0f}-{wavelength[-1]:.0f} nm")
print("\n".join(classification_lines))
print(f"Optimized RGB bands: {rgb_wavelengths.round(1)} nm")
print(f"Saved: {output_path / 'minimal_rgb_optimized.png'}")
print(f"Saved: {output_path / 'minimal_roi_selection.png'}")
print(f"Saved: {output_path / 'minimal_roi_spectra.png'}")
print(f"Saved: {output_path / 'minimal_roi_reference_comparison.png'}")
print(f"Saved: {output_path / 'minimal_reference_spectra.png'}")
