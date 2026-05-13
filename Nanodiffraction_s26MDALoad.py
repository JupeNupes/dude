"""
Nanodiffraction Script — standalone scan loader.

Extracts the core data-loading logic from dude.py's Scan_Load method,
free of any GTK/GUI dependencies, so it can be called from notebooks
or scripts.
"""

import os
import numpy as np
import h5py
import hdf5plugin
import fabio
from readMDA import readMDA


class ScanResult:
    """Container for everything Scan_Load produces.

    Attributes
    ----------
    data : list
        The raw readMDA output (list of scan dimensions + header dict).
    ndim : int
        Number of spatial scan dimensions (1 or 2).
    detector_names : list[str]
        Names of detector channels in the last MDA dimension.
    xrf_data : np.ndarray or list
        XRF spectra array (empty list when unavailable).
    xrf_in_mda : bool
    xrf_in_netcdf : bool
    h5 : h5py.File or None
        Open HDF5 handle for Eiger data (caller must close).
    image_list : np.ndarray or None
        Sorted array of TIFF filenames (non-Eiger detectors).
    image_path : str or None
        Directory containing the TIFF images.
    image_index_min : int
    image_index_max : int
    nbin : int
        Repeat-measurement stride.
    ntiff : int or None
        Expected number of TIFF images.
    tif_index : np.ndarray or None
        Numeric indices parsed from TIFF filenames.
    nrow : int or None
        Number of rows (y points) for 2D scans.
    ncol : int or None
        Number of columns (x points) for 2D scans.
    positioner_x : np.ndarray or None
        Spatial positions along the fast (x/column) scan axis.
        For 2D scans this is a 2D array of shape (nrow, ncol).
        For 1D scans this is a 1D array.
    positioner_y : np.ndarray or None
        Spatial positions along the slow (y/row) scan axis.
        For 2D scans this is a 2D array of shape (nrow, ncol).
        None for 1D scans.
    """

    def __init__(self):
        self.data = None
        self.ndim = 0
        self.detector_names = []
        self.xrf_data = []
        self.xrf_in_mda = False
        self.xrf_in_netcdf = False
        self.h5 = None
        self.image_list = None
        self.image_path = None
        self.image_index_min = 0
        self.image_index_max = -1
        self.nbin = 1
        self.ntiff = None
        self.tif_index = None
        self.nrow = None
        self.ncol = None
        self.positioner_x = None
        self.positioner_y = None


def load_scan(mdapath, *,
              scan_number=None,
              mda_folder=None,
              h5_folder=None,
              image_folder=None,
              eiger_enabled=True,
              pilatus_enabled=False,
              xrf_mode=False):
    """Load a nanodiffraction scan from an MDA file and its associated images.

    This is a standalone version of ``MyMainWindow.Scan_Load`` from dude.py.
    All GUI / GTK logic has been removed; the function returns a
    :class:`ScanResult` object instead of mutating widget state.

    Parameters
    ----------
    mdapath : str
        Full path to the ``.mda`` file **or** just the filename (in which
        case *mda_folder* must also be supplied so that sibling folders
        like ``h5/`` and ``Images/`` can be resolved).
    scan_number : str or int, optional
        The scan number used to locate companion HDF5 / fluorescence
        files (e.g. ``"0150"``).  If *None* it is extracted from
        *mdapath* automatically.
    mda_folder : str, optional
        Path to the ``mda/`` directory.  Required when *xrf_mode* is
        True with NetCDF data, or when folder-relative paths must be
        resolved.  If *None* it is inferred from *mdapath*.
    h5_folder : str, optional
        Path to the ``h5/`` directory containing Eiger HDF5 files.
        If *None* it is inferred as ``<parent of mda_folder>/h5``.
    image_folder : str, optional
        Path to the ``Images/`` directory containing per-scan TIFF
        sub-folders.  If *None* it is inferred as
        ``<parent of mda_folder>/Images``.
    eiger_enabled : bool
        Whether the Eiger detector (HDF5) path should be used.
    pilatus_enabled : bool
        Whether Pilatus naming conventions apply to TIFFs.
    xrf_mode : bool
        Whether to look for XRF fluorescence data (in-MDA or NetCDF).

    Returns
    -------
    ScanResult
        A plain object carrying all loaded data.  When *eiger_enabled*
        is True the ``.h5`` attribute holds an **open** ``h5py.File``
        that the caller is responsible for closing.
    """

    result = ScanResult()

    # --- resolve folders ---------------------------------------------------
    if os.path.isabs(mdapath):
        mdapath_full = mdapath
    else:
        if mda_folder is None:
            raise ValueError("mdapath is relative but mda_folder was not supplied")
        mdapath_full = os.path.join(mda_folder, mdapath)

    if mda_folder is None:
        mda_folder = os.path.dirname(mdapath_full)

    parent = os.path.abspath(os.path.join(mda_folder, os.pardir))
    if h5_folder is None:
        h5_folder = os.path.join(parent, "h5")
    if image_folder is None:
        image_folder = os.path.join(parent, "Images")

    # --- scan number -------------------------------------------------------
    if scan_number is None:
        scan_number = os.path.splitext(os.path.basename(mdapath_full))[0].split("_")[-1]
    scan_number = str(int(scan_number))  # strip leading zeros for matching

    # --- read MDA ----------------------------------------------------------
    result.data = readMDA(mdapath_full, verbose=0, maxdim=3 if xrf_mode else 2)
    data = result.data

    result.xrf_data = []
    result.xrf_in_mda = data[0]["dimensions"][-1] == 2048
    ndim = len(data) - 2 if result.xrf_in_mda else len(data) - 1
    result.ndim = ndim
    result.xrf_in_netcdf = False

    if ndim == 2:
        result.nrow = data[0]["dimensions"][0]
        result.ncol = data[0]["dimensions"][1]
        print("Scan {0} loaded: {1} points in x, {2} points in y".format(
            scan_number, result.ncol, result.nrow))

        # --- build 2D positioner grids (mirrors dude.py logic) ---
        # Fast axis (inner/column positions) — data[2].p[0].data is a 2D list
        pos_fast = np.zeros((result.nrow, result.ncol))
        datatmp = np.array(data[2].p[0].data)
        pos_fast[:datatmp.shape[0]] = datatmp  # handles interrupted scans
        result.positioner_x = pos_fast

        # Slow axis (outer/row positions) — data[1].p[0].data is a 1D list
        pos_slow = np.ones((result.nrow, result.ncol))
        pos_slow *= np.array(data[1].p[0].data)[:, np.newaxis]
        result.positioner_y = pos_slow

    elif ndim == 1:
        print("Scan {0} loaded: {1} points in x".format(
            scan_number, data[0]["dimensions"][0]))

        # --- 1D positioner ---
        result.positioner_x = np.array(data[1].p[0].data)

    # --- XRF ---------------------------------------------------------------
    if xrf_mode:
        if result.xrf_in_mda:
            nchan = data[ndim + 1].nd
            xrf_list = []
            for i in range(nchan):
                xrf_list.append(data[ndim + 1].d[i].data)
            result.xrf_data = np.array(xrf_list).reshape(nchan, -1, 2048)
        else:
            try:
                import netCDF4
            except ImportError:
                print("netCDF4 not available — skipping NetCDF XRF loading")
            else:
                fluo_folder = os.path.join(parent, "fluo")
                f_fluo = [f for f in os.listdir(fluo_folder)
                          if f.startswith("scan_{0}".format(scan_number))]
                if len(f_fluo):
                    f_fluo = f_fluo[0]
                    netcdffile = netCDF4.Dataset(os.path.join(fluo_folder, f_fluo), "r")
                    ny = data[0]["dimensions"][0]
                    nx = data[0]["dimensions"][1]
                    result.xrf_data = np.swapaxes(
                        np.swapaxes(
                            (netcdffile.variables['array_data'][:, :, 256:]
                             .reshape(ny, 2, 124, 256 + 2048 * 4)[:, :, :nx, 256:]),
                            1, 2
                        ).reshape(-1, 8, 2048),
                        0, 1
                    )
                    netcdffile.close()
                    result.xrf_in_netcdf = True

    # --- detector names ----------------------------------------------------
    result.detector_names = []
    for i in range(data[ndim].nd):
        dname = data[ndim].d[i].name
        result.detector_names.append(dname.replace("s26_eiger_cnm", "eiger"))

    # --- image loading (non-XRF branch) ------------------------------------
    if xrf_mode and (result.xrf_in_mda or result.xrf_in_netcdf):
        # Nothing more to do for XRF — caller can inspect xrf_data directly.
        pass
    else:
        if eiger_enabled:
            h5_filename = [f.name for f in os.scandir(h5_folder)
                           if "scan_{0}".format(scan_number) in f.name]
            if len(h5_filename) == 1:
                h5_filename = h5_filename[0]
            elif len(h5_filename) > 1:
                h5_index = np.array([int(f.split(".")[0].split("_")[-1])
                                     for f in h5_filename])
                h5_filename = np.take(h5_filename, h5_index.argsort())[0]
            else:
                print("No HDF5 file found for scan_{0} in {1}".format(
                    scan_number, h5_folder))
                result.image_index_max = -1
                return result

            try:
                result.h5 = h5py.File(os.path.join(h5_folder, h5_filename), 'r')
            except Exception as e:
                print("Failed to open HDF5: {0}".format(e))
                result.image_index_max = -1
            else:
                result.image_index_max = result.h5['entry/data/data'].shape[0] - 1
                result.image_index_min = 0
                result.nbin = 1
        else:
            # --- TIFF-based detectors ---
            scan_dir = os.path.splitext(os.path.basename(mdapath_full))[0].split("_")[-1].lstrip("0")
            result.image_path = os.path.join(image_folder, scan_dir)

            if pilatus_enabled:
                image_list = [f.name for f in os.scandir(result.image_path)
                              if '.tif' in f.name and 'il' in f.name]
            else:
                image_list = [f.name for f in os.scandir(result.image_path)
                              if '.tif' in f.name and 'il' not in f.name]

            tif_index = np.array([int(fn.split(".")[0].split("_")[-1])
                                  for fn in image_list])
            image_list = np.take(np.array(image_list), tif_index.argsort())
            tif_index = np.sort(tif_index)

            result.image_list = image_list
            result.tif_index = tif_index
            result.nbin = 1

            mda_index = None  # may be set below for 2-D scans

            if ndim == 2:
                ntiff = data[0]['dimensions'][0] * data[0]['dimensions'][1]
                for i in range(data[ndim].nd):
                    dname = data[ndim].d[i].name
                    if "FileNumber" in dname:
                        mda_index = np.array(data[2].d[i].data).flatten()
                if mda_index is not None:
                    print(mda_index.max(), tif_index.max(),
                          mda_index.min(), tif_index.min())
            else:
                ntiff = data[0]['dimensions'][0]

            result.ntiff = ntiff

            # ---- reconcile image count mismatches (mirrors dude.py logic) ----
            if len(image_list) != ntiff:
                print("Expecting {0} files but got {1} instead".format(
                    ntiff, len(image_list)))
                if ndim == 1:
                    if ntiff > len(image_list):
                        print("Scan appears interrupted — no correction attempted")
                    else:
                        if len(image_list) % ntiff == 0:
                            result.nbin = max(1, int(len(image_list) / ntiff))
                            print("Scan has {0} repeated measurement(s) per point".format(result.nbin))
                            result.image_list = image_list[0::result.nbin]
                        else:
                            print("Files may be in the wrong folder. Please correct manually 1!")
                else:
                    if len(image_list) > ntiff:
                        result.nbin = max(1, int(len(image_list) / ntiff))
                        if len(image_list) % ntiff:
                            result.nbin += 1
                        if result.nbin > 1:
                            if not len(image_list) % ntiff:
                                print("Scan has {0} repeated measurement(s) per point".format(result.nbin))
                                result.image_list = image_list[0::result.nbin]
                            else:
                                if mda_index is not None and mda_index.min() == tif_index.min():
                                    inc_index = mda_index[1:] - mda_index[:-1]
                                    hiccup_index = (inc_index == 2).nonzero()[0]
                                    if inc_index.max() == 2:
                                        print("Scan has {0} repeated measurement(s) per point, but {1} images missing".format(
                                            result.nbin, ntiff * result.nbin - len(image_list)))
                                        if inc_index[hiccup_index + 1].sum() == 0:
                                            print("Experiencing {0} hiccups. Using tif index.".format((inc_index == 2).sum()))
                                            mda_index[hiccup_index + 1] -= 1
                                            result.image_list = image_list[(mda_index - mda_index.min()).astype(int)]
                                            result.image_list = result.image_list[0::result.nbin]
                                        else:
                                            print("Experiencing {0} hiccups, unable to correct.".format((inc_index == 2).sum()))
                                    elif inc_index.min() < 0 or inc_index.max() > 2:
                                        print("This should not happen. Please correct manually 5!")
                                    elif inc_index.min() == 0 and (inc_index == 0).sum() == mda_index.shape[0] - tif_index.shape[0] + hiccup_index.shape[0]:
                                        print("Lost {0} images — trusting mda index.".format(
                                            (inc_index == 0).sum() - hiccup_index.shape[0]))
                                        result.image_list = image_list[0::result.nbin]
                                    else:
                                        print("Dunno what to do 7!")
                                else:
                                    print("Dunno what to do 6!")
                        else:
                            if mda_index is not None and mda_index.max() < tif_index.max():
                                if mda_index.min() == tif_index.min():
                                    print("Ignoring {0} extra images in the folder.".format(
                                        tif_index.max() - mda_index.max()))
                                    mask = tif_index <= mda_index.max()
                                    result.image_list = image_list[mask]
                                    result.tif_index = tif_index[mask]
                                    tif_index = result.tif_index
                                else:
                                    print("Files may be in the wrong folder. Please correct manually 2!")
                            else:
                                if mda_index is not None and mda_index.max() == tif_index.max() and mda_index.min() == tif_index.min() + 1:
                                    print("Known 2019 bug: mda index shifted by 1 vs tif index; extra point at start")
                                    result.image_list = image_list[:-1]
                                else:
                                    print("This should not happen. Please correct manually 3!")
                    else:
                        if mda_index is not None and mda_index.min() == tif_index.min():
                            inc_index = mda_index[1:] - mda_index[:-1]
                            hiccup_index = (inc_index == 2).nonzero()[0]
                            if inc_index.min() == 0 and (inc_index == 0).sum() == mda_index.shape[0] - tif_index.shape[0] + hiccup_index.shape[0]:
                                print("Lost {0} images — trusting mda index.".format(
                                    (inc_index == 0).sum() - hiccup_index.shape[0]))
                                result.image_list = image_list[(mda_index - mda_index.min()).astype(int)]
                            elif inc_index.max() == 2:
                                if inc_index[hiccup_index + 1].sum() == 0:
                                    print("Experiencing {0} hiccups. Using tif index.".format((inc_index == 2).sum()))
                                    mda_index[hiccup_index + 1] -= 1
                                    result.image_list = image_list[mda_index]
                                else:
                                    print("Experiencing {0} hiccups, unable to correct.".format((inc_index == 2).sum()))
                            elif inc_index.min() < 0 or inc_index.max() > 2:
                                print("This should not happen. Please correct manually 4!")
                            else:
                                print("Dunno what to do 8!")
                        else:
                            if mda_index is not None and mda_index.min() == tif_index.min() + 1:
                                print("Known 2019 bug: mda index shifted by 1 vs tif index (interrupted scan)")
                            else:
                                print("Files may be in the wrong folder. Please correct manually!")

            if tif_index.shape[0]:
                result.image_index_min = tif_index.min()
                result.image_index_max = tif_index.max()

    return result


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------

def get_image(result, index=None):
    """Read a single detector frame from an already-loaded scan.

    Parameters
    ----------
    result : ScanResult
        Object returned by :func:`load_scan`.
    index : int, optional
        Frame index.  Defaults to ``result.image_index_min``.

    Returns
    -------
    np.ndarray
        2-D detector image.
    """
    if index is None:
        index = result.image_index_min

    frame = int((index - result.image_index_min) / result.nbin)
    if result.h5 is not None:
        return result.h5['/entry/data/data'][frame].astype(np.int32)
    elif result.image_list is not None and result.image_path is not None:
        return fabio.open(os.path.join(result.image_path,
                                       result.image_list[frame])).data
    else:
        raise RuntimeError("No image data available in this ScanResult")


# ---------------------------------------------------------------------------
# Main — quick test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    file_path = filedialog.askopenfilename(
        title="Select MDA file",
        initialdir='/Users/smith.scott/Documents/Postdoc/Data/S26_20260211/',
        filetypes=[("MDA files", "*.mda"), ("All files", "*.*")]
    )
    root.destroy()

    if file_path:
        print(f"Selected file: {file_path}")
        result = load_scan(file_path)
        print(f"Loaded {result.ndim}-D scan with {len(result.detector_names)} detectors")
        print(f"Image index range: {result.image_index_min} — {result.image_index_max}")
        if result.h5 is not None:
            print(f"HDF5 dataset shape: {result.h5['entry/data/data'].shape}")
            result.h5.close()
    else:
        print("No file selected.")