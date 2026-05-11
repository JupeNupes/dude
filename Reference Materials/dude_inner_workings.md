# DUDE — Diffraction User Data Explorer

## Reference Guide: Inner Workings

> **Version**: 2.7.0  
> **Author**: Daniel Z. Tao (danielzt12)  
> **Stack**: Python 3, GTK 3, Matplotlib, NumPy, SciPy, FabIO, h5py, netCDF4

---

## Table of Contents

1. [High-Level Overview](#1-high-level-overview)
2. [Repository Structure](#2-repository-structure)
3. [Application Entry Point & Startup](#3-application-entry-point--startup)
4. [The MDA Data Format](#4-the-mda-data-format)
5. [GUI Architecture](#5-gui-architecture)
6. [Core Data Flow: Opening a Folder & Loading Scans](#6-core-data-flow)
7. [Detector Image System](#7-detector-image-system)
8. [1D and 2D Scan Plotting](#8-1d-and-2d-scan-plotting)
9. [Region of Interest (ROI) System](#9-region-of-interest-roi-system)
10. [XRF Mode](#10-xrf-mode)
11. [Settings & Toggles](#11-settings--toggles)
12. [Add-on Module System](#12-add-on-module-system)
13. [Google Drive / Logbook Integration](#13-google-drive--logbook-integration)
14. [Python Console](#14-python-console)
15. [Helper Modules](#15-helper-modules)
16. [Key State Variables](#16-key-state-variables)
17. [Mouse Interaction Model](#17-mouse-interaction-model)

---

## 1. High-Level Overview

DUDE is a desktop GUI application for visualizing synchrotron diffraction data collected at APS Sector 26. It reads **MDA** (Multi-Dimensional Archive) scan files, displays associated **detector images** (TIFF or HDF5/Eiger), and provides tools for interactive ROI analysis, XRF spectroscopy, and 2D spatial mapping.

The application is a single-window GTK 3 app with embedded Matplotlib canvases. Everything lives inside one monolithic class: `MyMainWindow`.

### What it does at the highest level:

1. **Browse** a folder of MDA scan files (displayed in a sortable tree view)
2. **Load** scan metadata and detector channel data from MDA files
3. **Display** individual detector images (TIFF/HDF5) tied to each scan point
4. **Plot** 1D line scans or 2D spatial maps from detector channels or custom ROIs
5. **Analyze** data using ROI integration, hot-pixel masking, XRF spectra, and add-on tools

---

## 2. Repository Structure

| File | Role |
|---|---|
| `dude.py` | Main application (~2470 lines). Contains `MyMainWindow` class and `main()`. |
| `readMDA.py` | MDA file parser. Reads the XDR-encoded binary scan format into Python objects. |
| `misc_dude.py` | Small helper functions: coordinate transforms, motor name lookup, parallel loaders. |
| `Shortcuts.py` | Add-on: user-defined shortcut scripts window. |
| `PowderHelper.py` | Add-on: powder diffraction analysis tools. |
| `ShiftCorrection.py` | Add-on: image shift/drift correction. |
| `ShowMetadata.py` | Add-on: HDF5 metadata browser. |
| `Analyze5D.py` | Add-on: 5D dataset analysis tools. |
| `run_dude.sh` | Launch script (sets macOS dynamic library paths for GTK). |
| `requirements.txt` | Python dependencies. |

---

## 3. Application Entry Point & Startup

```python
# dude.py, lines 2466-2471
if __name__ == "__main__":
    MyMainWindow()
    main()          # calls Gtk.main() — enters the GTK event loop
```

### `__init__` Sequence (lines 1644–2455)

The constructor does the following in order:

1. **Google Auth** (lines 1647–1670): If running at the beamline (checks for `/home/sector26/.../token.json`), authenticates with Google Drive and Google Sheets for logbook integration. Sets `beamline = True/False`.

2. **Default State** (lines 1672–1682): Initializes key variables:
   - `self.cm = "viridis"` — default colormap
   - `self.dimY = 1062`, `self.dimX = 1028` — default Eiger detector dimensions
   - Mouse state flags: `Image_Rectangle_Drawing`, `Image_Zoomed`, `Image_Zooming`, `Image_Panning`

3. **MDA File Browser** (lines 1684–1711): Creates a `Gtk.TreeView` backed by a `ListStore` with 12 columns (scan number, motor names, ranges, number of points, count time, filepath, checkbox).

4. **ROI / Detector Controls** (lines 1714–1843): Builds the left-side toolbox with combo boxes for detector channel (DET), ROI, and monitor (MON) selection, plus spin buttons for ROI coordinates and hot-pixel masking.

5. **1D Plot Canvas** (lines 1845–1886): Creates a Matplotlib `Figure` + `FigureCanvas` for line plots.

6. **2D Plot Canvas** (lines 1888–1975): Creates a Matplotlib `Figure` + `FigureCanvas` for 2D image/map plots, with Log/Auto/Scale toggle buttons and Vmin/Vmax sliders.

7. **Detector Image Canvas** (lines 1977–2072): Creates the main detector image display with zoom, pan, auto-scale, log-scale controls, and Vmin/Vmax sliders.

8. **XRF Display** (lines 2074–2253): Creates dual-axis XRF spectrum display with per-element radio buttons (1–8) and ROI sliders.

9. **Console** (lines 2256–2265): A `Gtk.Entry` with auto-completion for executing arbitrary Python commands.

10. **Menu Bar** (lines 2267–2417): File (Open/Refresh/Quit), Settings (Detector/Colormap/Sparse/ShowAngle/DirtyFix/PumpProbe/XRF), Add-ons, Help.

11. **Window Assembly** (lines 2419–2455): Packs everything into the final layout and shows the window.

---

## 4. The MDA Data Format

MDA files are binary files using XDR (External Data Representation) encoding, produced by the EPICS `sscan` record at APS beamlines.

### `readMDA.py` — Key Classes

| Class | Purpose |
|---|---|
| `scanDim` | One dimension of scan data (rank, npts, positioners, detectors, triggers) |
| `scanPositioner` | Motor/positioner info: name, unit, step mode, readback, and data array |
| `scanDetector` | Detector channel info: name, unit, and data array |
| `scanTrigger` | Trigger info: name and command value |

### `readMDA()` Return Value

```python
data = readMDA(fname, verbose=0, maxdim=2)
# Returns a list:
#   data[0] = dict with scan metadata (rank, dimensions, filename, environment PVs)
#   data[1] = scanDim for dimension 1 (outer scan, e.g. slow axis)
#   data[2] = scanDim for dimension 2 (inner scan, e.g. fast axis) — if 2D
#   data[3] = scanDim for dimension 3 — if 3D (used for XRF 2048-channel data)
```

### How DUDE determines scan dimensionality

```python
# If the last dimension has 2048 points, it's XRF spectral data, not a scan axis
xrf_in_mda = data[0]["dimensions"][-1] == 2048
ndim = len(data) - 2 if xrf_in_mda else len(data) - 1
```

---

## 5. GUI Architecture

The main window layout (simplified):

```
┌──────────────────────────────────────────────────────────────────┐
│ Menu Bar: [File] [Settings] [Add-ons] [Help]                    │
├─────────────────────────────┬────────────────────────────────────┤
│ MDA File       │ Toolbox    │                                    │
│ TreeView       │ (DET/ROI/  │   Detector Image Canvas            │
│ (scan list)    │  MON/Mask) │   (or XRF Spectra in XRF mode)     │
│                │            │                                    │
├────────────────┴────────────┤                                    │
│ Plot Notebook               │                                    │
│ ┌─────────┬────────┐        │                                    │
│ │ 1D Plot │ 2D Map │        │                                    │
│ └─────────┴────────┘        │                                    │
├─────────────────────────────┴────────────────────────────────────┤
│ Python Console Entry                                             │
└──────────────────────────────────────────────────────────────────┘
```

**Key Notebooks** (tab containers with hidden tabs, switched programmatically):
- `Plot_Notebook`: Switches between 1D plot (page 0) and 2D map (page 1) based on scan dimensionality
- `Image_Notebook`: Switches between detector image view (page 0) and XRF spectra view (page 1)

---

## 6. Core Data Flow

### Step 1: Open Folder — `FileDialog_Construction(flag=0)`

When the user opens an MDA folder via File → Open (F2):

1. A GTK folder chooser dialog is shown
2. The selected path is stored in `self.MDA_folder`
3. Derived paths are set:
   - `self.Image_folder = .../Images/` (sibling of `mda/`)
   - `self.h5_folder = .../h5/` (sibling of `mda/`)
4. Attempts to open a Google Sheet named `S26_<experiment_name>` for logbook integration
5. Calls `Folder_Open()`

### Step 2: Scan Directory — `Folder_Open()` → `Folder_Scan()`

1. Lists all `.mda` files, sorts by scan number
2. Uses `multiprocessing.Pool` to call `mda_loader()` (from `misc_dude.py`) in parallel for each file
3. `mda_loader()` reads just the metadata from each MDA file and returns a list:
   ```
   [scan_number, motor1_name, min, max, npts, motor2_name, min, max, npts, count_time, filepath, False]
   ```
4. Each result is appended to `self.MDA_File_ListStore` (the TreeView's data model)

### Step 3: Select a Scan — `Scan_TreeView_Selection_Changed()`

When a row is clicked in the TreeView:

1. Gets the MDA file path from the selected row
2. Resets ROI and image state
3. Calls `Scan_Load(mdapath)`

### Step 4: Load Scan Data — `Scan_Load()`

This is the **central data loading method** (~200 lines). It:

1. **Reads MDA data**: `self.data = readMDA(mdapath, verbose=0, maxdim=...)`
2. **Populates detector list**: Fills `self.MDA_Det_store` with detector channel names
3. **Handles XRF data**: If in XRF mode, loads from MDA (3D dimension) or from netCDF files in a `fluo/` folder
4. **Locates image files**:
   - **Eiger (HDF5)**: Finds matching `scan_NNNN_*.h5` in `self.h5_folder`, opens with `h5py`
   - **TIFF**: Lists `.tif` files in `self.Image_folder/<scan_num>/`, sorts by index
5. **Handles file count mismatches** (lines 1442–1522): Contains extensive logic to deal with:
   - Interrupted scans (fewer images than expected)
   - Binned/repeated measurements (more images than expected)
   - Index offset bugs ("This bug is new in 2019...")
   - "Hiccups" where the file numbering skips
6. **Sets image slider range**: Updates `Image_Plot_HScale_Adjustment` bounds
7. **Triggers initial plot**: Via `Scan_ToolBox_Plot_Changed()`

---

## 7. Detector Image System

### Loading a Single Image — `Image_Plot_HScale_Changed()`

Triggered when the image slider moves. Loads one detector frame:

- **Eiger**: `self.h5["/entry/instrument/detector/data"][index]`
- **TIFF**: `fabio.open(filepath).data`
- Applies hot-pixel mask if enabled (`self.custom_mask`)
- Sets the image data on the Matplotlib `AxesImage`, updates normalization (linear or log)

### Loading All Images — `LoadAllImages()`

Loads the **entire image stack** into memory for ROI analysis:

- **Eiger + Sparse**: Reads HDF5 data in blocks of 1000 frames, converts each block to a `scipy.sparse.csr_matrix`, stacks vertically. This flattens each frame from (1062, 1028) to a 1D vector of length 1,093,736.
- **Eiger + Dense**: `self.h5["/entry/instrument/detector/data"][()]` — loads entire dataset
- **TIFF**: Uses `multiprocessing.Pool` with `image_loader()` to read all TIFFs in parallel

After loading, enables the ROI analysis buttons (Add, Sum, custom ROI entry).

### Intensity Scaling

- **Auto-scale** (`Image_AutoScale_toggled`): Sets Vmin/Vmax to the data min/max
- **Manual scale** (`Image_Vscale_Changed`): Uses slider values
- **Log scale** (`Image_Log_ToggleButton`): Switches between `colors.LogNorm` and `colors.Normalize`

---

## 8. 1D and 2D Scan Plotting

### `Scan_ToolBox_Plot_Changed()` — The Master Plot Refresh

This method is called whenever the detector channel, ROI, or monitor selection changes. It determines `ndim` and branches:

#### 1D Scans (`ndim == 1`)
- **X data**: `self.data[1].p[0].data` — the positioner (motor) values
- **Y data** (detector mode): `self.data[1].d[selected_detector].data`
- **Y data** (ROI mode): Sums the loaded image stack within the ROI bounds
- **Monitor normalization**: If MON toggle is active, divides Y by the selected monitor channel (normalized to its mean)
- **Dirty fix**: If enabled, shifts detector data by one index to correct a known first-point bug for Eiger/QMPX3 detectors
- Calls `Scan_Plot1D(xdata, ydata)` which plots with `-o` markers

#### 2D Scans (`ndim == 2`)
- **X axes**: Outer motor (`data[1].p[0]`) and inner motor (`data[2].p[0]`)
- **Y data**: Reshaped into a 2D array matching the scan grid
- **Spiral scans**: If motor name is `"ao01"`, uses scatter plot (`Scan_PlotSpiral`) instead of `imshow`
- Calls `Scan_Plot2D(xdata1, xdata2, ydata)` which uses `imshow` with the selected colormap

### Plot Interaction

- **1D plot**: Right-clicking or dragging over the plot highlights the nearest data point and updates the detector image slider to show that scan point's image
- **2D map**: Right-clicking on a pixel shows coordinates in the P0 label and jumps the image slider; left-clicking shows distance from P0 in the P1 label

---

## 9. Region of Interest (ROI) System

### Defining an ROI

ROIs can be defined two ways:

1. **Mouse drawing**: Right-click and drag on the detector image canvas creates a blue rectangle. Coordinates are automatically populated into the spin buttons.
2. **Spin buttons**: Manual entry of XMin, XMax, YMin, YMax coordinates.

### Saving an ROI — `Scan_ToolBox_CustomROI_Added()`

Stores the ROI name and coordinates in `self.CustomROI_store` (a `Gtk.ListStore`). ROIs persist during the session and can be selected from the ROI combo box.

### Summing an ROI — `Scan_ToolBox_CustomROI_Summed()`

Integrates pixel intensity within the ROI across all scan points:

- **Sparse matrix path**: Creates a binary mask, finds non-zero indices, sums columns of the sparse matrix
- **Dense array path**: Slices the 3D array `image[:, ymin:ymax+1, xmin:xmax+1]` and sums over detector axes
- For 1D scans: produces a line plot of integrated intensity vs. motor position
- For 2D scans: produces a spatial map of integrated intensity

### Pump-Probe ROI Handling

When pump-probe mode is active, the ROI sum is computed twice — once for the original ROI position and once mirrored by ±550 pixels vertically (the two halves of the Eiger detector correspond to pump and probe frames).

---

## 10. XRF Mode

Activated via Settings → XRF Mode. Switches the UI to show XRF spectra instead of detector images.

### Data Sources

- **MDA-embedded XRF**: When the last MDA dimension has 2048 channels, it's XRF data
- **netCDF files**: Reads from a `fluo/` directory, reshaping the fluorescence detector array into `(8_elements, n_spectra, 2048_channels)`

### Display

- Two XRF axes (XRF1 and XRF2) showing spectra for individual or averaged elements
- Radio buttons select which element(s) to display (1–7 average, individual 1–8)
- Vmin/Vmax sliders control the energy range displayed
- **Sum button**: Shows the summed spectrum across all scan points
- **ROI buttons**: Integrates within a selected energy window and plots vs. scan position (1D) or as a spatial map (2D)

---

## 11. Settings & Toggles

| Setting | Variable | Purpose |
|---|---|---|
| **Detector** | `self.dimX/Y`, `self.eiger_enabled`, `self.pilatus_enabled` | Sets detector pixel dimensions. Options: 1028×1062 (Eiger), 516×516, 515×515, 256×256, 487×195 (Pilatus), 1360×1024 |
| **Colormap** | `self.cm` | Matplotlib colormap for images/maps. Options: jet, binary, gray, hot, coolwarm, viridis, plasma, inferno, magma, cividis |
| **Sparse Matrix** | `self.sparse_enabled` | When ON (default), stores loaded image stacks as `scipy.sparse.csr_matrix` to save memory |
| **Show Angle** | `self.show_angle` | Replaces pixel X/Y coordinates with 2θ/γ angles (reads from `Analysis/gamma.csv` and `twotheta.csv`) |
| **Dirty Fix** | `self.dirty_fix` | Corrects a known index-offset bug where Eiger/QMPX3 detector data is shifted by one point relative to positioner data |
| **Pump Probe** | `self.pump_probe` | Enables specialized image handling: interleaves alternating frames and splits the detector into two halves |
| **XRF Mode** | `self.xrf_mode` | Switches between diffraction image display and XRF spectrum display |

---

## 12. Add-on Module System

Add-ons are separate Python modules, each defining a `MainWindow` class that receives the main `dude` instance (`self`) as a parent reference. This gives them full access to loaded data.

Each add-on follows the same lazy-init pattern:

```python
def Addon_XYZ(self, widget):
    try:
        self.XYZ.win.present()       # If window exists, bring to front
    except:
        self.XYZ = XYZ.MainWindow(self)  # Otherwise create it
```

| Add-on | Module | Purpose |
|---|---|---|
| User Shortcuts | `Shortcuts.py` | Custom user-defined analysis scripts |
| Show Metadata | `ShowMetadata.py` | Browse HDF5 file metadata keys |
| Run PtychoLib | `PtychoLib.py` | Ptychography reconstruction (optional import) |
| Probe Library | `ProbeLibrary.py` | Probe function management (requires PtychoLib) |
| Powder Helper | `PowderHelper.py` | Powder diffraction ring analysis |
| Shift Correction | `ShiftCorrection.py` | Image registration / drift correction |
| Analyze 5D | `Analyze5D.py` | 5D dataset analysis (scanning + diffraction) |

---

## 13. Google Drive / Logbook Integration

Only active when running at the beamline (`beamline = True`). Upload buttons ("+") appear next to each plot canvas.

### `Upload_To_Logbook(flag)`

1. Saves the current plot as a JPEG to `/tmp/`
2. Uploads to a specific Google Drive folder via the Drive API
3. Sets the file to public read access
4. Inserts a row in a Google Sheet with an `=IMAGE()` formula pointing to the uploaded file
5. Sets the row height to 200px for image visibility

---

## 14. Python Console

A `Gtk.Entry` at the bottom of the window that executes arbitrary Python code via `exec()`.

- **`clear`**: Clears command history
- **`history`**: Prints all previous commands
- **Any other text**: Executed as Python with access to `self` (the main window instance)
- **Shift+Up/Down**: Navigates command history
- Has auto-completion via `Gtk.EntryCompletion`

This is powerful — you can type things like:
```python
self.image.shape
self.data[1].p[0].name
np.save("/tmp/mydata.npy", self.image)
```

---

## 15. Helper Modules

### `misc_dude.py`

| Function | Purpose |
|---|---|
| `Display2Data(Axe, x, y)` | Converts display pixel coordinates to data coordinates using Matplotlib's inverse transform |
| `MotorMNE_Parser(longname)` | Maps EPICS PV names (e.g. `26idcnpi:m10.VAL`) to human-readable motor mnemonics (e.g. `fomx`) |
| `image_loader(fname)` | Loads a single TIFF file via `fabio.open()` — used as a worker for multiprocessing |
| `mda_loader(fname)` | Reads MDA metadata for the file browser — used as a worker for multiprocessing |
| `MotorMNE` | Dictionary mapping common EPICS PV names to short motor names for Sector 26 |

---

## 16. Key State Variables

### Data State

| Variable | Type | Description |
|---|---|---|
| `self.data` | list | Current MDA data: `[dict, scanDim, ...]` |
| `self.image` | ndarray or csr_matrix or None | Loaded image stack (all frames) |
| `self.h5` | h5py.File | Currently open HDF5 (Eiger) file |
| `self.xrf_data` | ndarray | XRF spectral data array |
| `self.AllImagesLoaded` | bool | Whether `LoadAllImages()` has been called |

### Path State

| Variable | Description |
|---|---|
| `self.MDA_folder` | Path to the `mda/` directory |
| `self.Image_folder` | Path to the `Images/` directory |
| `self.h5_folder` | Path to the `h5/` directory |
| `self.image_path` | Path to the current scan's TIFF subfolder |
| `self.image_list` | Sorted array of TIFF filenames for the current scan |

### Display State

| Variable | Description |
|---|---|
| `self.cm` | Current colormap string |
| `self.dimX`, `self.dimY` | Current detector dimensions |
| `self.nbin` | Number of repeated measurements per scan point |
| `self.image_index_min/max` | Range of valid image indices for slider |
| `self.custom_mask` | Boolean mask array for hot-pixel removal |

### Interaction Flags

| Variable | Description |
|---|---|
| `self.Image_Zoomed` | Whether the image is currently zoomed in |
| `self.Image_Zooming` | Whether a zoom drag is in progress |
| `self.Image_Panning` | Whether a pan drag is in progress |
| `self.Image_Rectangle_Drawing` | Whether an ROI rectangle drag is in progress |
| `self.non_animated_background` | Cached canvas for efficient rectangle animation (blitting) |

---

## 17. Mouse Interaction Model

The detector image canvas uses a state-machine pattern for mouse interactions:

### Button Press (`Image_Canvas_Button_Pressed`)

| Button | Condition | Action |
|---|---|---|
| Left (1) | Not zoomed | Start **zoom** mode |
| Left (1) | Already zoomed | Start **pan** mode (cursor changes to FLEUR) |
| Middle (2) | Zoomed | **Zoom out** to full view |
| Right (3) | Any | Start **ROI drawing** mode |

On press, the `button_press_event` handler is disconnected and replaced with a `button_release_event` handler. This prevents conflicting actions during drag.

### Mouse Move (`Image_Canvas_Mouse_Hover`)

- Updates the pixel coordinate / intensity label
- If zooming or drawing ROI: uses **blitting** for efficient rectangle animation
  - Saves the clean background once, restores it each frame, draws the rectangle artist, then blits
- If panning: shifts the axes limits to follow the mouse

### Button Release (`Image_Canvas_Button_Released`)

- **ROI drawing**: Finalizes the rectangle, updates spin button values
- **Zooming**: Sets new axis limits maintaining aspect ratio
- **Panning**: Resets cursor to crosshair
- Re-connects the `button_press_event` handler

### Scroll Wheel (`Image_Canvas_Button_Scrolled`)

Zooms in/out centered on the mouse position. Scales both X and Y ranges by `(1 - 0.1 * scroll_step)`.

---

## Appendix: Expected Directory Structure

DUDE expects the experiment data to be organized as:

```
<experiment_root>/
├── mda/                    ← MDA scan files (selected via File → Open)
│   ├── prefix_0001.mda
│   ├── prefix_0002.mda
│   └── ...
├── h5/                     ← Eiger HDF5 files (auto-discovered)
│   ├── scan_0001_data_000001.h5
│   └── ...
├── Images/                 ← TIFF image folders (auto-discovered)
│   ├── 1/                  ← TIFFs for scan 1
│   │   ├── prefix_00001.tif
│   │   └── ...
│   └── 2/
├── fluo/                   ← XRF netCDF files (optional)
│   ├── scan_0001.nc
│   └── ...
└── Analysis/               ← Pre-computed angle maps (optional)
    ├── gamma.csv
    └── twotheta.csv
```

---

## Appendix: Sparse Matrix Strategy

When `sparse_enabled = True` (default), `LoadAllImages()` stores the image stack as a **scipy CSR sparse matrix** with shape `(n_frames, 1062*1028)`. Each 2D frame is flattened to 1D.

This is effective because diffraction patterns are typically very sparse (most pixels are zero or near-zero).

**ROI extraction with sparse matrices**: Rather than slicing a 3D array, the code creates a binary index mask, finds the non-zero column indices, then sums across those columns:

```python
indice0 = np.zeros((1062, 1028))
indice0[ymin:ymax+1, xmin:xmax+1] = 1
indice0 = indice0.flatten().nonzero()[0]    # column indices within the sparse matrix
ydata = self.image[:, indice0].sum(1)       # sum across ROI for each frame
```

---

## Appendix: The "Dirty Fix"

The `dirty_fix` toggle (enabled by default) corrects a known synchronization bug where Eiger or QMPX3 detector data is recorded one scan point ahead of the positioner data:

```python
# For 1D scans: shift detector data backward by one
ydata[:-1] = ydata[1:]

# For 2D scans: flatten, shift, reshape
ydata_flat = ydata.flatten()
ydata_flat[:-1] = ydata_flat[1:]
ydata = ydata_flat.reshape(ydata.shape)
```

This only applies when the selected detector channel name contains "eiger" or "QMPX3".
