"""
nanoDiff_BatchProcess.py
========================
Batch processing pipeline for nanodiffraction analysis.

Automates the workflow from nanoDiff_BasicAnaly.ipynb:
  - Iterates over a configurable list of scans
  - For each scan and each ROI, computes:
      * Total Intensity
      * X-Center of Mass (weighted centroid)
      * Y-Center of Mass (weighted centroid)
      * (Optional) 2D Gaussian fit parameters
  - Saves outputs as 600 DPI PNGs, SVGs, .npy arrays, and JSON metadata

Usage:
    python nanoDiff_BatchProcess.py
"""


# --- --- Imports --- ---


import os
import sys
import time
import json
from pathlib import Path
from datetime import datetime

import numpy as np
import h5py
import hdf5plugin
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for batch processing
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from matplotlib_scalebar.scalebar import ScaleBar

from Nanodiffraction_s26MDALoad import load_scan

# Try to import cmcrameri for scientific colormaps (managua_r, berlin)
try:
    import cmcrameri  # registers colormaps with matplotlib
except ImportError:
    pass

# Optional: scipy for Gaussian fitting
try:
    from scipy.optimize import curve_fit
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


# --- --- Configuration --- ---


# Data root — base path to the experiment data directory
DATA_ROOT = "/Users/smith.scott/Documents/Postdoc/Data/S26_20260211"

# Output root — base directory for all saved results
OUTPUT_ROOT = "/Users/smith.scott/Documents/Postdoc/Data/S26_20260211/analysis"

# Scan list — supports ranges and individual numbers:
#   e.g. list(range(150, 161)) + [175, 180]

SCAN_LIST = list(range(85, 89)) + [92,98,100] + list(range(103,111)) + list(range(112,126)) + [127,137,140,141,142,143,144] #Non bost scans
# SCAN_LIST = list(range(78, 85)) + [89,90,91,94,97] #Bost Scans
# SCAN_LIST = [89] # Bost Scan
# SCAN_LIST = [126] Bost Scan

# ROI definitions — hardcoded, embedded in JSON output for traceability.
# Each ROI is a dict with name and detector-pixel-space coordinates.
# Multiple ROIs per scan are supported.
ROI_LIST = [
    {"name": "ROI_LeftPk",  "xmin": 110, "xmax": 180, "ymin": 660, "ymax": 755},
    {"name": "ROI_MidPk",   "xmin": 199, "xmax": 264, "ymin": 673, "ymax": 756},
    {"name": "ROI_RightPk", "xmin": 271, "xmax": 331, "ymin": 673, "ymax": 756},
    {"name": "ROI_xPtRing",  "xmin": 421, "xmax": 518, "ymin": 568, "ymax": 986},
]
    # Add more ROIs here as needed:
    # {"name": "ROI_Secondary", "xmin": 600, "xmax": 700, "ymin": 200, "ymax": 300},

# Projective Distortion Correction — applies 1/sin(theta) stretch to x-axis.
# Set False for scans that don't need it.
ENABLE_PROJECTIVE_CORRECTION = True

# Boustrophedon Scan Correction — reverses every other row to account
# for snake (boustrophedon) raster scan geometry.
ENABLE_BOUSTROPHEDON_CORRECTION = False

# Optional 2D Gaussian fitting — set True to enable per-frame fitting
ENABLE_GAUSS_FIT = False

# Output settings
DPI = 600


# Matplotlib global font size
plt.rcParams.update({'font.size': 14})

# Colormaps — fixed per analysis type
COLORMAPS = {
    "intensity":       "viridis",
    "xcom":            "RdBu_r",
    "ycom":            "BrBG_r",
    # Gaussian fit parameter colormaps (used when ENABLE_GAUSS_FIT is True)
    "gauss_amplitude": "viridis",
    "gauss_mu_x":      "RdBu_r",
    "gauss_mu_y":      "BrBG_r",
    "gauss_sigma_x":   "managua_r",
    "gauss_sigma_y":   "berlin",
    "gauss_offset":    "RdYlBu_r",
}

# --- --- Data Loading --- ---


def load_scan_data(scan_num):
    """Load MDA + HDF5 data for a single scan number.

    Returns a data-bundle dict, or None if the scan cannot be loaded.
    """
    mda_path = Path(DATA_ROOT) / "data" / "mda" / f"26idbSOFT_{scan_num:04d}.mda"

    if not mda_path.exists():
        print(f"  ⚠ MDA file not found: {mda_path} — skipping scan {scan_num}")
        return None

    print(f"  Loading scan {scan_num}...")
    t0 = time.time()

    mda = load_scan(str(mda_path))

    if mda.ndim != 2:
        print(f"  ⚠ Scan {scan_num} is {mda.ndim}D, not 2D — skipping")
        return None

    if mda.h5 is None:
        print(f"  ⚠ No HDF5 data for scan {scan_num} — skipping")
        return None


    # --- ---


    # Load all frames into RAM via context manager
    h5_filename = mda.h5.filename
    mda.h5.close()  # Close the handle opened by load_scan

    with h5py.File(h5_filename, 'r') as f:
        all_frames = f['/entry/data/data'][()]
        sample_theta = float(f['/entry/instrument/26-ID-C/SAMTH'][()])


    # --- ---


    # Hot pixel masking
    all_frames[all_frames > 2**31] = 0


    # --- ---


    # Projective distortion correction: 1/sin(theta)
    if ENABLE_PROJECTIVE_CORRECTION:
        x_stretched = mda.positioner_x / np.sin(np.radians(sample_theta))
    else:
        x_stretched = mda.positioner_x.copy()
    dx = float(np.mean(np.abs(np.diff(x_stretched[0]))))

    elapsed = time.time() - t0
    proj_str = "ON" if ENABLE_PROJECTIVE_CORRECTION else "OFF"
    print(f"  Loaded {all_frames.shape[0]} frames in {elapsed:.1f}s "
          f"(θ={sample_theta}°, {mda.nrow}×{mda.ncol} grid, "
          f"proj. corr. {proj_str})")

    return {
        "all_frames": all_frames,
        "mda": mda,
        "x_stretched": x_stretched,
        "positioner_y": mda.positioner_y,
        "dx": dx,
        "sample_theta": sample_theta,
        "nrow": mda.nrow,
        "ncol": mda.ncol,
    }


# --- --- ROI Analysis --- ---


def analyze_roi(data_bundle, roi):
    """Compute Total Intensity, X-CoM, Y-CoM maps for a single ROI.

    Returns (intensity_map, xcom_map, ycom_map), each shaped (nrow, ncol).
    """
    all_frames = data_bundle["all_frames"]
    nrow = data_bundle["nrow"]
    ncol = data_bundle["ncol"]

    xmin, xmax = roi["xmin"], roi["xmax"]
    ymin, ymax = roi["ymin"], roi["ymax"]

    n_frames = nrow * ncol

    intensity_map = np.zeros((nrow, ncol))
    xcom_map = np.full((nrow, ncol), np.nan)
    ycom_map = np.full((nrow, ncol), np.nan)

    x_indices = np.arange(xmin, xmax, dtype=np.float64)
    y_indices = np.arange(ymin, ymax, dtype=np.float64)

    for i in range(n_frames):
        row = i // ncol
        col = i % ncol

        frame = all_frames[i]
        roi_frame = frame[ymin:ymax, xmin:xmax].astype(np.float64)

        total = np.sum(roi_frame)
        intensity_map[row, col] = total

        if total > 0:
            # X-CoM: sum ROI along y-axis (axis=0) → intensity per x-column
            col_sums = np.sum(roi_frame, axis=0)
            xcom_map[row, col] = np.sum(x_indices * col_sums) / total

            # Y-CoM: sum ROI along x-axis (axis=1) → intensity per y-row
            row_sums = np.sum(roi_frame, axis=1)
            ycom_map[row, col] = np.sum(y_indices * row_sums) / total

    # Boustrophedon (snake scan) correction — reverse every other row
    if ENABLE_BOUSTROPHEDON_CORRECTION:
        intensity_map[1::2, :] = intensity_map[1::2, ::-1]
        xcom_map[1::2, :] = xcom_map[1::2, ::-1]
        ycom_map[1::2, :] = ycom_map[1::2, ::-1]

    return intensity_map, xcom_map, ycom_map


# --- --- 2D Gaussian Fitting (optional) --- ---


def _gaussian_2d(xy, amplitude, mu_x, mu_y, sigma_x, sigma_y, offset):
    """2D Gaussian model (no rotation).

    f(x, y) = A * exp(-((x-mu_x)^2/(2*sigma_x^2) +
                         (y-mu_y)^2/(2*sigma_y^2))) + offset
    """
    x, y = xy
    return offset + amplitude * np.exp(
        -((x - mu_x) ** 2 / (2 * sigma_x ** 2)
          + (y - mu_y) ** 2 / (2 * sigma_y ** 2))
    )


def fit_gaussian_2d(data_bundle, roi):
    """Fit a 2D Gaussian to each frame's ROI sub-frame.

    Returns a dict of parameter maps {name: ndarray(nrow, ncol)}.
    Parameters: amplitude, mu_x, mu_y, sigma_x, sigma_y, offset.
    Failed fits are stored as np.nan.
    """
    if not SCIPY_AVAILABLE:
        print("  ⚠ scipy not available — skipping Gaussian fitting")
        return {}

    all_frames = data_bundle["all_frames"]
    nrow = data_bundle["nrow"]
    ncol = data_bundle["ncol"]

    xmin, xmax = roi["xmin"], roi["xmax"]
    ymin, ymax = roi["ymin"], roi["ymax"]

    n_frames = nrow * ncol

    param_names = ["amplitude", "mu_x", "mu_y", "sigma_x", "sigma_y", "offset"]
    param_maps = {name: np.full((nrow, ncol), np.nan) for name in param_names}

    # Meshgrid for ROI pixel coordinates
    x = np.arange(xmin, xmax, dtype=np.float64)
    y = np.arange(ymin, ymax, dtype=np.float64)
    xx, yy = np.meshgrid(x, y)
    xy = (xx.ravel(), yy.ravel())

    # Bounds
    bounds_lower = [0, xmin, ymin, 0.1, 0.1, 0]
    bounds_upper = [np.inf, xmax, ymax, xmax - xmin, ymax - ymin, np.inf]

    n_failed = 0
    print(f"  Fitting 2D Gaussians to {n_frames} frames...")
    t0 = time.time()

    for i in range(n_frames):
        row = i // ncol
        col = i % ncol

        frame = all_frames[i]
        roi_frame = frame[ymin:ymax, xmin:xmax].astype(np.float64)

        median_val = float(np.median(roi_frame))
        max_val = float(np.max(roi_frame))

        p0 = [
            max(max_val - median_val, 1.0),
            (xmin + xmax) / 2.0,
            (ymin + ymax) / 2.0,
            (xmax - xmin) / 6.0,
            (ymax - ymin) / 6.0,
            max(median_val, 0.0),
        ]

        try:
            popt, _ = curve_fit(
                _gaussian_2d, xy, roi_frame.ravel(),
                p0=p0, bounds=(bounds_lower, bounds_upper),
                maxfev=5000,
            )
            for j, name in enumerate(param_names):
                param_maps[name][row, col] = popt[j]
        except Exception:
            n_failed += 1

    elapsed = time.time() - t0
    print(f"  Gaussian fitting complete in {elapsed:.1f}s "
          f"({n_failed} failed fits out of {n_frames})")

    return param_maps


# --- --- Visualization --- ---


def plot_single_map(data_bundle, map_data, title, cbar_label, cmap,
                    output_path_base, colorscale):
    """Plot a single 2D map with colorbar + histogram, save as PNG and SVG."""
    x_stretched = data_bundle["x_stretched"]
    positioner_y = data_bundle["positioner_y"]
    dx = data_bundle["dx"]
    if colorscale == 'log':
        norm = 'log'
    else:
        norm = 'linear'

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_aspect('equal')
    ax.set_title(title)
    ax.set_xlabel(r"Sample X $\mu$m")
    ax.set_ylabel(r"Sample Y $\mu$m")

    img = ax.pcolormesh(x_stretched, positioner_y, map_data, cmap=cmap, norm=norm)

    # --- --- Colorbar with Superimposed Histogram --- ---


    cbar = fig.colorbar(img, ax=ax, label=cbar_label, orientation='horizontal')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')

    hist_ax = cbar.ax.twinx()

    valid_data = map_data[~np.isnan(map_data)].flatten()
    if len(valid_data) > 0:
        hist_ax.hist(valid_data, bins=100,
                     range=(np.nanmin(map_data), np.nanmax(map_data)),
                     color='black', alpha=0.3, density=True)

    hist_ax.set_yticks([])
    hist_ax.set_ylim(bottom=0)


    # --- --- Scalebar --- ---


    scalebar = ScaleBar(dx=dx, units="um", location="lower right",
                        color="black", box_alpha=0.0)
    ax.add_artist(scalebar)


    # --- --- Save --- ---


    output_path_base = Path(output_path_base)
    fig.savefig(str(output_path_base.with_suffix('.png')),
                dpi=DPI, bbox_inches='tight')
    fig.savefig(str(output_path_base.with_suffix('.svg')),
                bbox_inches='tight')
    plt.close(fig)

    print(f"    Saved {output_path_base.stem}.png / .svg")


# --- --- Data Saving --- ---


def save_results(output_dir, roi, results_dict, scan_num, data_bundle):
    """Save .npy arrays and JSON metadata for one scan-ROI combination."""
    output_dir = Path(output_dir)

    # Save .npy files
    maps_meta = {}
    for name, map_data in results_dict.items():
        npy_file = f"{name}.npy"
        np.save(str(output_dir / npy_file), map_data)
        maps_meta[name] = {
            "min": float(np.nanmin(map_data)),
            "max": float(np.nanmax(map_data)),
            "file": npy_file,
        }

    # Build metadata
    metadata = {
        "scan_number": scan_num,
        "roi": roi,
        "dimensions": {
            "nrow": int(data_bundle["nrow"]),
            "ncol": int(data_bundle["ncol"]),
        },
        "sample_theta": data_bundle["sample_theta"],
        "dx_um": data_bundle["dx"],
        "projective_correction": ENABLE_PROJECTIVE_CORRECTION,
        "boustrophedon_correction": ENABLE_BOUSTROPHEDON_CORRECTION,
        "gaussian_fit_enabled": ENABLE_GAUSS_FIT,
        "maps": maps_meta,
        "timestamp": datetime.now().isoformat(),
    }

    json_path = output_dir / "metadata.json"
    with open(str(json_path), 'w') as f:
        json.dump(metadata, f, indent=4)

    print(f"    Saved metadata.json + {len(results_dict)} .npy files")


def plot_detector_roi(data_bundle, roi, output_dir):
    """Save a detector image showing the ROI rectangle overlay.

    Picks a frame from the center of the scan for a representative view.
    """
    all_frames = data_bundle["all_frames"]
    nrow = data_bundle["nrow"]
    ncol = data_bundle["ncol"]

    # Pick a frame near the center of the scan
    center_row = nrow // 2
    center_col = ncol // 2
    center_idx = center_row * ncol + center_col
    frame = all_frames[center_idx].astype(np.float64)
    frame = np.ma.masked_less_equal(frame, 0)  # mask zeros → transparent

    xmin, xmax = roi["xmin"], roi["xmax"]
    ymin, ymax = roi["ymin"], roi["ymax"]
    roi_name = roi["name"]

    fig, ax = plt.subplots(figsize=(16, 12))
    cmap = plt.cm.viridis.copy()
    cmap.set_bad(color='white', alpha=0)  # masked pixels → white/transparent
    im = ax.imshow(frame, norm=LogNorm(), cmap=cmap)
    ax.set_xlim(50,550)
    ax.set_ylim(500,950)
    fig.colorbar(im, ax=ax, label='Intensity (Log Scale)')
    ax.set_title(f"Detector Image — {roi_name}  "
                 f"(Row {center_row}, Col {center_col})")

    # Draw ROI rectangle
    roi_rect = patches.Rectangle(
        (xmin, ymin), xmax - xmin, ymax - ymin,
        linewidth=1.5, edgecolor='red', facecolor='none',
        linestyle='--', label=roi_name,
    )
    ax.add_patch(roi_rect)
    ax.legend(loc='upper right', fontsize=12)

    # Save
    output_path = Path(output_dir) / "detector_roi"
    fig.savefig(str(output_path.with_suffix('.png')),
                dpi=DPI, bbox_inches='tight')
    fig.savefig(str(output_path.with_suffix('.svg')),
                bbox_inches='tight')
    plt.close(fig)

    print(f"    Saved detector_roi.png / .svg")


# --- --- Main Driver --- ---


def process_scans():
    """Main batch processing loop."""
    print(f"Nanodiffraction Batch Processor")
    print(f"Data root: {DATA_ROOT}")
    print(f"Scans to process: {SCAN_LIST}")
    print(f"ROIs: {[r['name'] for r in ROI_LIST]}")
    print(f"Projective correction: {'ON' if ENABLE_PROJECTIVE_CORRECTION else 'OFF'}")
    print(f"Boustrophedon correction: {'ON' if ENABLE_BOUSTROPHEDON_CORRECTION else 'OFF'}")
    print(f"Gaussian fitting: {'ENABLED' if ENABLE_GAUSS_FIT else 'DISABLED'}")
    print(f"Output: {OUTPUT_ROOT}")

    overall_t0 = time.time()

    for scan_num in SCAN_LIST:
        print(f"\n{'=' * 60}")
        print(f"Processing Scan {scan_num:04d}...")

        try:
            data_bundle = load_scan_data(scan_num)
            if data_bundle is None:
                continue

            for roi in ROI_LIST:
                roi_name = roi['name']
                output_dir = Path(OUTPUT_ROOT) / f"Scan_{scan_num:04d}" / roi_name
                output_dir.mkdir(parents=True, exist_ok=True)


                # --- --- Core Analysis --- ---


                print(f"  Analyzing {roi_name}...")
                t0 = time.time()
                intensity, xcom, ycom = analyze_roi(data_bundle, roi)
                print(f"  Core analysis done in {time.time() - t0:.1f}s")

                # Detector image with ROI overlay
                plot_detector_roi(data_bundle, roi, output_dir)


                # --- --- Plot each map separately --- ---


                plot_single_map(
                    data_bundle, intensity,
                    f"Total Intensity — {roi_name}",
                    "Total Intensity (counts)",
                    COLORMAPS["intensity"],
                    output_dir / "intensity",
                    'log',
                )

                plot_single_map(
                    data_bundle, xcom,
                    f"X-CoM — {roi_name}",
                    "X Centroid (Weighted, px)",
                    COLORMAPS["xcom"],
                    output_dir / "xcom",
                    'linear',
                )

                plot_single_map(
                    data_bundle, ycom,
                    f"Y-CoM — {roi_name}",
                    "Y Centroid (Weighted, px)",
                    COLORMAPS["ycom"],
                    output_dir / "ycom",
                    'linear',
                )

                results = {
                    "intensity": intensity,
                    "xcom": xcom,
                    "ycom": ycom,
                }


                # --- --- Optional Gaussian Fitting --- ---


                if ENABLE_GAUSS_FIT:
                    gauss_params = fit_gaussian_2d(data_bundle, roi)

                    gauss_labels = {
                        "amplitude": "Gaussian Amplitude",
                        "mu_x": "Gaussian X Center (px)",
                        "mu_y": "Gaussian Y Center (px)",
                        "sigma_x": r"Gaussian $\sigma_x$ (px)",
                        "sigma_y": r"Gaussian $\sigma_y$ (px)",
                        "offset": "Gaussian Offset (counts)",
                    }

                    for param_name, param_map in gauss_params.items():
                        cmap_key = f"gauss_{param_name}"
                        plot_single_map(
                            data_bundle, param_map,
                            f"Gauss {param_name} — {roi_name}",
                            gauss_labels.get(param_name, param_name),
                            COLORMAPS[cmap_key],
                            output_dir / f"gauss_{param_name}",
                            'linear',
                        )
                        results[f"gauss_{param_name}"] = param_map


                # --- --- Save --- ---


                save_results(output_dir, roi, results, scan_num, data_bundle)

                print(f"  ✓ {roi_name} complete")

        except Exception as e:
            print(f"  ✗ Error processing Scan {scan_num}: {e}")
            import traceback
            traceback.print_exc()
            continue

    total_elapsed = time.time() - overall_t0
    print(f"\n{'=' * 60}")
    print(f"Batch processing complete in {total_elapsed:.1f}s")


if __name__ == "__main__":
    process_scans()
