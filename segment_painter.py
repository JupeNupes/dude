#!/usr/bin/env python3
"""
Segment Painter — Click and paint pixels into up to 3 segments (R/G/B).
Loads .npy 2D maps exported from DUDE, displays in grayscale,
and lets you manually assign pixels to segments by painting.

Usage:
    python segment_painter.py                         # launch, then click Load
    python segment_painter.py path/to/file.npy        # load directly
"""

import sys
import os
import json
import subprocess
import time
import numpy as np
from scipy import sparse
from scipy.optimize import curve_fit
import h5py
try:
    import hdf5plugin
except ImportError:
    pass
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors, gridspec
from matplotlib.widgets import Button, TextBox


# ─────────────────────────────────────────────────────────────────────
class SegmentPainter:

    COLORS = {
        1: (1.0, 0.2, 0.2, 0.55),   # Red
        2: (0.2, 0.5, 1.0, 0.55),   # Blue
        3: (0.2, 0.8, 0.3, 0.55),   # Green
    }
    LABELS = {1: 'Red', 2: 'Blue', 3: 'Green'}

    def __init__(self, npy_path=None):
        # State
        self.map_data = None
        self.segments = None          # 2D int array: 0=none, 1=R, 2=B, 3=G
        self.active_color = 1         # currently painting: 1, 2, or 3
        self.show_underlay = True
        self.show_overlay = True
        self.painting = False
        self._painted_pixels = set()  # pixels painted in this drag stroke
        self.meta = {}
        self.npy_path = None
        self.scan_folder = ''
        self.data_prefix = ''         # e.g. 'Orth_Mid'
        self.output_dir = ''          # e.g. scan_193/segmentation_Orth_Mid/
        self._h5_offset = 0           # +1 when Eiger sync bug is active

        # H5 / detector state
        self.image_stack = None
        self.dimY = 1062
        self.dimX = 1028
        self.det_vmin = 0.1
        self.det_vmax = 100.0
        self.det_use_sum = True       # True=sum, False=avg
        self._det_fig = None
        self._det_ims = []
        self._last_fit_results = {}   # seg_id -> popt list
        self._fit_fig = None
        self._fit_axes = ()

        self._build_ui()

        if npy_path and os.path.isfile(npy_path):
            self._load_file(os.path.abspath(npy_path))

        # Manual event loop — prevents macOS segfault on sub-window close
        try:
            while plt.fignum_exists(self.fig.number):
                self.fig.canvas.flush_events()
                time.sleep(0.05)
        except KeyboardInterrupt:
            pass
        sys.exit(0)

    # ─── UI ──────────────────────────────────────────────────────────
    def _build_ui(self):
        plt.ion()  # interactive mode — no blocking show()
        self.fig = plt.figure(figsize=(10, 8))
        self.fig.canvas.manager.set_window_title('Segment Painter')

        # Main image axes (leave room on right for buttons)
        self.ax = self.fig.add_axes([0.05, 0.05, 0.72, 0.90])
        self.ax.set_title('Segment Painter — click Load to begin')
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        # Placeholder images
        dummy = np.zeros((2, 2))
        self.im_under = self.ax.imshow(
            dummy, cmap='gray', origin='lower', aspect='equal',
            interpolation='nearest')
        self.im_overlay = self.ax.imshow(
            np.zeros((2, 2, 4)), origin='lower', aspect='equal',
            interpolation='nearest')

        # ── Buttons (right side) ─────────────────────────────────────
        bx = 0.82          # left edge of buttons
        bw, bh = 0.14, 0.05  # width, height
        gap = 0.015

        def btn(y, label, color, callback):
            ax_b = self.fig.add_axes([bx, y, bw, bh])
            b = Button(ax_b, label, color=color, hovercolor='#EEEEEE')
            b.on_clicked(callback)
            return b

        y = 0.90
        self._btn_load = btn(y, 'Load', '#E0E0E0', self._on_load)

        y -= bh + gap * 3
        self._btn_red   = btn(y, '● Red (1)', '#FFCDD2', self._sel_red)
        y -= bh + gap
        self._btn_blue  = btn(y, '● Blue (2)', '#BBDEFB', self._sel_blue)
        y -= bh + gap
        self._btn_green = btn(y, '● Green (3)', '#C8E6C9', self._sel_green)

        y -= bh + gap * 3
        self._btn_toggle = btn(y, 'Toggle Underlay', '#E0E0E0',
                               self._on_toggle)
        y -= bh + gap
        self._btn_toggle_ov = btn(y, 'Toggle Overlay', '#E0E0E0',
                                  self._on_toggle_overlay)

        y -= bh + gap * 3
        self._btn_detector = btn(y, 'Show Detector', '#D1C4E9',
                                 self._on_show_detector)

        y -= bh + gap * 3
        self._btn_save = btn(y, 'Save Segments', '#FFF9C4',
                             self._on_save)
        y -= bh + gap
        self._btn_load_seg = btn(y, 'Load Segments', '#E0E0E0',
                                 self._on_load_segments)

        # ── vmin / vmax text boxes ────────────────────────────────────
        y -= bh + gap * 3
        self.fig.text(bx, y + 0.015, 'vmin:', fontsize=8, color='#555')
        ax_vmin = self.fig.add_axes([bx + 0.04, y, bw - 0.04, 0.03])
        self._tb_vmin = TextBox(ax_vmin, '', initial='0.1')
        self._tb_vmin.on_submit(self._on_vmin_submit)

        y -= 0.045
        self.fig.text(bx, y + 0.015, 'vmax:', fontsize=8, color='#555')
        ax_vmax = self.fig.add_axes([bx + 0.04, y, bw - 0.04, 0.03])
        self._tb_vmax = TextBox(ax_vmax, '', initial='100')
        self._tb_vmax.on_submit(self._on_vmax_submit)

        # Active color indicator
        self._color_text = self.fig.text(
            bx + bw / 2, 0.90 + bh + 0.01,
            'Brush: Red', ha='center', fontsize=10, fontweight='bold',
            color='red')

        # Status text
        self._status = self.fig.text(
            bx + bw / 2, 0.05,
            '', ha='center', fontsize=8, color='#555555',
            style='italic')

        # ── Mouse events ─────────────────────────────────────────────
        c = self.fig.canvas
        c.mpl_connect('button_press_event', self._on_press)
        c.mpl_connect('motion_notify_event', self._on_motion)
        c.mpl_connect('button_release_event', self._on_release)
        c.mpl_connect('key_press_event', self._on_key)

    # ─── File loading ────────────────────────────────────────────────
    def _on_load(self, event=None):
        script = (
            'set f to POSIX path of (choose file '
            'with prompt "Select .npy file")\nreturn f'
        )
        try:
            result = subprocess.run(
                ['osascript', '-e', script],
                capture_output=True, text=True, timeout=120)
            path = result.stdout.strip()
            if path and os.path.isfile(path):
                self._load_file(os.path.abspath(path))
        except Exception as e:
            print(f"File dialog error: {e}")

    def _on_load_segments(self, event=None):
        """Import pre-existing segment masks from another folder.

        Only replaces self.segments — does NOT change data_prefix,
        output_dir, or any data-channel state.  Saving still uses
        the current data's context.
        """
        if self.segments is None:
            print("Load a data file first before importing segments.")
            return
        # Pick a folder containing mask_red/blue/green.npy
        try:
            res = subprocess.run(
                ['osascript', '-e',
                 'set f to POSIX path of (choose folder '
                 'with prompt "Select folder with segment masks")'
                 '\nreturn f'],
                capture_output=True, text=True, timeout=120)
            folder = res.stdout.strip()
            if not folder:
                return
        except Exception as e:
            print(f"Folder dialog error: {e}")
            return

        # Try loading individual masks
        loaded = False
        new_seg = np.zeros(self.segments.shape, dtype=int)
        for seg_id, label in self.LABELS.items():
            mask_path = os.path.join(folder, f'mask_{label.lower()}.npy')
            if os.path.exists(mask_path):
                mask = np.load(mask_path).astype(bool)
                if mask.shape == self.segments.shape:
                    new_seg[mask] = seg_id
                    loaded = True
                    print(f"  Loaded mask_{label.lower()}.npy  "
                          f"({mask.sum()} px)")
                else:
                    print(f"  ⚠ mask_{label.lower()}.npy shape mismatch: "
                          f"{mask.shape} vs {self.segments.shape}")

        # Also try segment_map.npy as fallback
        if not loaded:
            map_path = os.path.join(folder, 'segment_map.npy')
            if os.path.exists(map_path):
                smap = np.load(map_path).astype(int)
                if smap.shape == self.segments.shape:
                    new_seg = smap
                    loaded = True
                    print(f"  Loaded segment_map.npy")

        if loaded:
            self.segments = new_seg
            self._last_fit_results = {}   # clear old fits
            self._refresh_overlay()
            self._update_status()
            self.fig.canvas.draw_idle()
            src = os.path.basename(folder.rstrip('/'))
            print(f"Segments imported from {src}/")
            print(f"  Saves will still go to: "
                  f"{os.path.basename(self.output_dir)}/")
        else:
            print(f"No valid mask files found in {folder}")

    @staticmethod
    def _extract_prefix(basename):
        """Extract the user prefix from a DUDE-exported .npy filename.

        Examples:
            'Orth_Mid_intensity_user_roi.npy'   → 'Orth_Mid'
            'FluoPt_intensity_mda_channel.npy'  → 'FluoPt'
            'custom_name.npy'                   → 'custom_name'
        """
        stem = basename.replace('.npy', '')
        # Known DUDE suffixes (longest first so greedy match works)
        suffixes = [
            '_intensity_user_roi',
            '_intensity_mda_channel',
            '_intensity_full_frame',
            '_com_x', '_com_y',
        ]
        for suf in suffixes:
            if stem.endswith(suf):
                return stem[:-len(suf)]
        return stem

    def _load_file(self, npy_path):
        self.npy_path = npy_path
        self.scan_folder = os.path.dirname(npy_path)
        basename = os.path.basename(npy_path)

        # Derive prefix and output directory
        self.data_prefix = self._extract_prefix(basename)
        self.output_dir = os.path.join(
            self.scan_folder, f'segmentation_{self.data_prefix}')

        # Load map
        self.map_data = np.load(npy_path).astype(float)
        print(f"Loaded: {basename}  shape={self.map_data.shape}")
        print(f"  Prefix: {self.data_prefix}")
        print(f"  Output: {os.path.basename(self.scan_folder)}/"
              f"segmentation_{self.data_prefix}/")

        # Load metadata
        meta_name = basename.replace('.npy', '_meta.json')
        meta_path = os.path.join(self.scan_folder, meta_name)
        self.meta = {}
        if os.path.exists(meta_path):
            with open(meta_path) as f:
                self.meta = json.load(f)

        # Detector dims from metadata
        dt = self.meta.get('detector_type', {})
        self.dimY = dt.get('dimY', 1062)
        self.dimX = dt.get('dimX', 1028)

        # ── Dirty fix / sync-bug detection ─────────────────────────
        dirty_fix = self.meta.get('dirty_fix', False)
        map_type = self.meta.get('map_type', '')
        det_ch = self.meta.get('detector_channel') or ''

        # Was the dirty fix already applied to this .npy?
        is_eiger_mda = ('MDA Channel' in map_type and
                        ('eiger' in det_ch.lower() or
                         'qmpx3' in det_ch.lower()))
        map_already_fixed = dirty_fix and is_eiger_mda

        # Is the map derived from H5 (ROI, CoM) and NOT fixed?
        h5_derived = ('User ROI' in map_type or
                      'CoM' in map_type or
                      'Full Frame' in map_type)
        needs_fix = dirty_fix and h5_derived

        if needs_fix:
            # Apply the same shift DUDE uses for Eiger MDA data
            flat = self.map_data.ravel()
            flat[:-1] = flat[1:]
            self.map_data = flat.reshape(self.map_data.shape)
            print(f"  ⚡ Applied dirty fix to H5-derived map "
                  f"(Eiger sync correction)")

        # H5 frame offset: +1 when Eiger sync bug is present
        if dirty_fix and (h5_derived or is_eiger_mda):
            self._h5_offset = 1
            print(f"  ⚡ H5 frame offset = +1 (Eiger sync bug active)")
        else:
            self._h5_offset = 0

        # Load the H5 diffraction stack
        self._load_h5()

        # Reset segment mask and fit results
        self.segments = np.zeros(self.map_data.shape, dtype=int)
        self._last_fit_results = {}

        # Close any open detector or fit windows to avoid stale ROI views
        if self._det_fig is not None and plt.fignum_exists(self._det_fig.number):
            plt.close(self._det_fig.number)
            self._det_fig = None
            self._det_ims = []
            self._det_axes = None
            
        if self._fit_fig is not None and plt.fignum_exists(self._fit_fig.number):
            plt.close(self._fit_fig.number)
            self._fit_fig = None
            self._fit_axes = ()

        # Redirect matplotlib's native save dialog to our output folder
        os.makedirs(self.output_dir, exist_ok=True)
        matplotlib.rcParams['savefig.directory'] = self.output_dir

        # Update underlay
        self.im_under.set_data(self.map_data)
        self.im_under.set_extent(
            [-0.5, self.map_data.shape[1] - 0.5,
             -0.5, self.map_data.shape[0] - 0.5])
        self.im_under.autoscale()
        self.im_under.set_visible(self.show_underlay)

        # Update overlay
        self._refresh_overlay()

        self.ax.set_xlim(-0.5, self.map_data.shape[1] - 0.5)
        self.ax.set_ylim(-0.5, self.map_data.shape[0] - 0.5)
        self.ax.set_title(f'{self.data_prefix}  [{basename}]')
        self.fig.canvas.manager.set_window_title(
            f'Segment Painter — {self.data_prefix}')
        self._update_status()
        self.fig.canvas.draw_idle()

    def _load_h5(self):
        """Load diffraction images from the H5 file as sparse matrices."""
        self.image_stack = None
        h5_file = self.meta.get('h5_file')
        mda_file = self.meta.get('mda_file')
        if not h5_file or not mda_file:
            print("No H5 info in metadata — detector view disabled")
            return
        exp_root = os.path.dirname(os.path.dirname(mda_file))
        h5_path = os.path.join(exp_root, 'h5', h5_file)
        if not os.path.exists(h5_path):
            print(f"H5 not found: {h5_path}")
            return
        print(f"Loading H5: {h5_file} …")
        h5 = h5py.File(h5_path, 'r')
        dd = h5["/entry/instrument/detector/data"]
        n = dd.shape[0]
        blk = 1000
        parts = []
        for i in range(0, n, blk):
            end = min(i + blk, n)
            b = dd[i:end][()]
            parts.append(sparse.csr_matrix(b.reshape(b.shape[0], -1)))
            print(f"  frames {i}–{end - 1} / {n - 1}")
        self.image_stack = sparse.vstack(parts)
        mp = self.map_data.size
        if self.image_stack.shape[0] > mp:
            self.image_stack = self.image_stack[:mp]
        print(f"  Stack: {self.image_stack.shape}")

    # ─── Overlay ─────────────────────────────────────────────────────
    def _refresh_overlay(self):
        if self.segments is None:
            return
        h, w = self.segments.shape
        rgba = np.zeros((h, w, 4), dtype=float)
        for seg_id, col in self.COLORS.items():
            mask = self.segments == seg_id
            rgba[mask] = col
        self.im_overlay.set_data(rgba)
        self.im_overlay.set_extent(
            [-0.5, w - 0.5, -0.5, h - 0.5])

    # ─── Painting ────────────────────────────────────────────────────
    def _paint_pixel(self, row, col, allow_toggle=False):
        """Paint a single pixel. Only toggles off on initial click, not drag."""
        if (row, col) in self._painted_pixels:
            return  # already handled in this stroke
        if self.segments[row, col] == self.active_color and allow_toggle:
            self.segments[row, col] = 0
        else:
            self.segments[row, col] = self.active_color
        self._painted_pixels.add((row, col))

    def _on_press(self, event):
        if event.inaxes == self.ax and event.button == 1:
            if self.segments is None:
                return
            self.painting = True
            self._painted_pixels = set()
            col = int(round(event.xdata))
            row = int(round(event.ydata))
            if 0 <= row < self.segments.shape[0] and \
               0 <= col < self.segments.shape[1]:
                self._paint_pixel(row, col, allow_toggle=True)
                self._refresh_overlay()
                self._update_status()
                self.fig.canvas.draw_idle()

    def _on_motion(self, event):
        if not self.painting or self.segments is None:
            return
        if event.inaxes != self.ax:
            return
        col = int(round(event.xdata))
        row = int(round(event.ydata))
        if 0 <= row < self.segments.shape[0] and \
           0 <= col < self.segments.shape[1]:
            self._paint_pixel(row, col, allow_toggle=False)
            self._refresh_overlay()
            self._update_status()
            self.fig.canvas.draw_idle()

    def _on_release(self, event):
        self.painting = False
        self._painted_pixels = set()

    # ─── Button callbacks ────────────────────────────────────────────
    def _sel_red(self, event):
        self.active_color = 1
        self._color_text.set_text('Brush: Red')
        self._color_text.set_color('red')
        self.fig.canvas.draw_idle()

    def _sel_blue(self, event):
        self.active_color = 2
        self._color_text.set_text('Brush: Blue')
        self._color_text.set_color('#2196F3')
        self.fig.canvas.draw_idle()

    def _sel_green(self, event):
        self.active_color = 3
        self._color_text.set_text('Brush: Green')
        self._color_text.set_color('green')
        self.fig.canvas.draw_idle()

    def _on_toggle(self, event):
        self.show_underlay = not self.show_underlay
        self.im_under.set_visible(self.show_underlay)
        self.fig.canvas.draw_idle()

    def _on_toggle_overlay(self, event):
        self.show_overlay = not self.show_overlay
        self.im_overlay.set_visible(self.show_overlay)
        self.fig.canvas.draw_idle()

    def _on_key(self, event):
        if event.key == '1':
            self._sel_red(None)
        elif event.key == '2':
            self._sel_blue(None)
        elif event.key == '3':
            self._sel_green(None)
        elif event.key == 't':
            self._on_toggle(None)
        elif event.key == 's':
            self._on_save(None)
        elif event.key == 'f':
            self._on_load(None)
        elif event.key == 'd':
            self._on_show_detector(None)
        elif event.key == 'm':
            self._on_load_segments(None)

    def _on_vmin_submit(self, text):
        try:
            self.det_vmin = float(text)
            self._update_detector()
        except ValueError:
            pass

    def _on_vmax_submit(self, text):
        try:
            self.det_vmax = float(text)
            self._update_detector()
        except ValueError:
            pass

    # ─── Detector viewer ─────────────────────────────────────────────
    def _on_show_detector(self, event):
        """Compute and show sum/avg detector images for each segment."""
        if self.image_stack is None:
            print("No H5 data loaded — cannot show detector images.")
            return
        if self.segments is None:
            print("No segments painted yet.")
            return
        self._build_detector_window()
        self._update_detector()
        self._det_fig.show()

    def _build_detector_window(self):
        """Create (or reshow) the detector popup with 3 linked panels."""
        if self._det_fig is not None and plt.fignum_exists(self._det_fig.number):
            return  # already open

        self._det_fig, axes = plt.subplots(
            1, 3, figsize=(18, 6), num='Detector — Segments',
            sharex=True, sharey=True)
        self._det_fig.canvas.mpl_connect(
            'close_event', self._on_det_close)
        self._det_fig.subplots_adjust(
            left=0.03, right=0.97, top=0.90, bottom=0.05,
            wspace=0.08)

        dummy = np.zeros((self.dimY, self.dimX))
        norm = colors.LogNorm(vmin=max(self.det_vmin, 0.1),
                              vmax=self.det_vmax)
        self._det_axes = axes
        self._det_ims = []
        seg_colors = ['red', '#2196F3', 'green']
        for i, (ax, label, sc) in enumerate(zip(
                axes, ['Red', 'Blue', 'Green'], seg_colors)):
            im = ax.imshow(dummy, cmap='viridis', origin='lower',
                           aspect='equal', norm=norm,
                           interpolation='nearest')
            ax.set_title(label, color=sc, fontweight='bold')
            self._det_ims.append(im)

        # Auto-zoom to ROI
        roi = self.meta.get('roi', {'xmin': 0, 'xmax': self.dimX - 1,
                                    'ymin': 0, 'ymax': self.dimY - 1})
        dx = max(roi['xmax'] - roi['xmin'], 1)
        dy = max(roi['ymax'] - roi['ymin'], 1)
        m = 0.20
        axes[0].set_xlim(roi['xmin'] - m * dx, roi['xmax'] + m * dx)
        axes[0].set_ylim(roi['ymax'] + m * dy, roi['ymin'] - m * dy)

        # Sum/Avg toggle button
        ax_btn = self._det_fig.add_axes([0.44, 0.93, 0.12, 0.05])
        self._det_btn_mode = Button(ax_btn, 'Mode: SUM',
                                    color='#E0E0E0',
                                    hovercolor='#BBDEFB')
        self._det_btn_mode.on_clicked(self._toggle_det_mode)

        # Fit Gaussian button
        ax_fit = self._det_fig.add_axes([0.58, 0.93, 0.14, 0.05])
        self._det_btn_fit = Button(ax_fit, 'Fit Gaussian',
                                   color='#C8E6C9',
                                   hovercolor='#A5D6A7')
        self._det_btn_fit.on_clicked(self._on_fit_gaussian)

        self._det_fig.canvas.mpl_connect('key_press_event',
                                         self._det_on_key)

    def _on_det_close(self, event):
        """When detector window is closed, just reset the reference."""
        self._det_fig = None
        self._det_ims = []
        self._det_axes = None

    def _toggle_det_mode(self, event=None):
        self.det_use_sum = not self.det_use_sum
        label = 'SUM' if self.det_use_sum else 'AVG'
        self._det_btn_mode.label.set_text(f'Mode: {label}')
        self._update_detector()

    def _det_on_key(self, event):
        if event.key == 'a':
            self._toggle_det_mode()

    def _update_detector(self):
        """Recompute summed/averaged detector images and update panels."""
        if self._det_fig is None or not plt.fignum_exists(self._det_fig.number):
            return
        if self.image_stack is None or self.segments is None:
            return

        n_fr = self.image_stack.shape[0]
        seg_flat = self._get_h5_indices()
        mode = 'sum' if self.det_use_sum else 'avg'
        vmin = max(self.det_vmin, 0.001)
        vmax = max(self.det_vmax, vmin + 0.1)

        for seg_id, im, ax in zip(
                [1, 2, 3], self._det_ims, self._det_axes):
            mask = seg_flat == seg_id
            npx = mask.sum()
            if npx > 0:
                d = self.image_stack[mask].sum(axis=0)
                d = np.asarray(d).reshape(self.dimY, self.dimX).astype(float)
                if not self.det_use_sum:
                    d = d / npx
            else:
                d = np.zeros((self.dimY, self.dimX))
            im.set_data(d)
            im.set_norm(colors.LogNorm(vmin=vmin, vmax=vmax))
            label = self.LABELS[seg_id]
            ax.set_title(f"{label} ({npx} px, {mode})",
                         fontweight='bold')

        self._det_fig.canvas.draw_idle()

    # ─── 2D Gaussian fitting ─────────────────────────────────────────
    @staticmethod
    def _gaussian2d(xy, amp, x0, y0, sx, sy, theta, offset):
        """Rotated 2D Gaussian."""
        x, y = xy
        a = np.cos(theta)**2 / (2*sx**2) + np.sin(theta)**2 / (2*sy**2)
        b = -np.sin(2*theta) / (4*sx**2) + np.sin(2*theta) / (4*sy**2)
        c = np.sin(theta)**2 / (2*sx**2) + np.cos(theta)**2 / (2*sy**2)
        return (offset + amp * np.exp(
            -(a*(x - x0)**2 + 2*b*(x - x0)*(y - y0) + c*(y - y0)**2)))

    def _on_fit_gaussian(self, event=None):
        """Fit 2D Gaussians to each segment's detector data within ROI."""
        if self._det_fig is None or self.image_stack is None:
            return

        # Use ROI bounds from metadata for fitting, fallback to view limits if missing
        roi = self.meta.get('roi')
        if roi:
            x0i = int(max(roi['xmin'], 0))
            x1i = int(min(roi['xmax'], self.dimX - 1))
            y0i = int(max(roi['ymin'], 0))
            y1i = int(min(roi['ymax'], self.dimY - 1))
        else:
            ax0 = self._det_axes[0]
            xlo, xhi = sorted(ax0.get_xlim())
            ylo, yhi = sorted(ax0.get_ylim())
            x0i, x1i = int(max(xlo, 0)), int(min(xhi, self.dimX - 1))
            y0i, y1i = int(max(ylo, 0)), int(min(yhi, self.dimY - 1))

        n_fr = self.image_stack.shape[0]
        seg_flat = self._get_h5_indices()

        # Coordinate grids for the ROI sub-image
        yy, xx = np.mgrid[y0i:y1i+1, x0i:x1i+1]
        xflat = xx.ravel().astype(float)
        yflat = yy.ravel().astype(float)

        fit_results = {}  # seg_id -> (popt, data_roi)
        seg_colors_3d = {1: 'red', 2: 'blue', 3: 'green'}

        for seg_id in [1, 2, 3]:
            mask = seg_flat == seg_id
            if mask.sum() == 0:
                continue
            # Get the detector image for this segment
            d = self.image_stack[mask].sum(axis=0)
            d = np.asarray(d).reshape(self.dimY, self.dimX).astype(float)
            if not self.det_use_sum:
                d = d / mask.sum()

            roi_data = d[y0i:y1i+1, x0i:x1i+1]
            zflat = roi_data.ravel()

            # Initial guesses
            total = zflat.sum()
            if total <= 0:
                continue
            cx = np.sum(xflat * zflat) / total
            cy = np.sum(yflat * zflat) / total
            amp0 = zflat.max()
            sx0 = max((x1i - x0i) / 4, 1)
            sy0 = max((y1i - y0i) / 4, 1)
            off0 = max(np.median(zflat), 0)
            p0 = [amp0, cx, cy, sx0, sy0, 0.0, off0]
            bounds_lo = [0, x0i, y0i, 0.5, 0.5, -np.pi, 0]
            bounds_hi = [amp0*10, x1i, y1i, (x1i-x0i), (y1i-y0i), np.pi, amp0]

            try:
                popt, _ = curve_fit(
                    self._gaussian2d, (xflat, yflat), zflat,
                    p0=p0, bounds=(bounds_lo, bounds_hi),
                    maxfev=10000)
                fit_results[seg_id] = (popt, roi_data)
                print(f"  {self.LABELS[seg_id]} fit: amp={popt[0]:.1f} "
                      f"x0={popt[1]:.1f} y0={popt[2]:.1f} "
                      f"sx={popt[3]:.2f} sy={popt[4]:.2f} "
                      f"theta={np.degrees(popt[5]):.1f}° off={popt[6]:.1f}")
            except Exception as e:
                print(f"  {self.LABELS[seg_id]} fit failed: {e}")

        if fit_results:
            # Store for saving later
            self._last_fit_results = {
                sid: popt.tolist() for sid, (popt, _) in fit_results.items()}
            # Auto-save fit params
            self._save_json()
            self._show_fit_3d(fit_results, xx, yy, seg_colors_3d)

    def _show_fit_3d(self, fit_results, xx, yy, seg_colors_3d):
        """Show 3D plots: raw data (left) and Gaussian fits (right)."""
        fig = plt.figure(figsize=(18, 8), num='2D Gaussian Fits')
        fig.clf()

        ax_data = fig.add_subplot(121, projection='3d')
        ax_fit  = fig.add_subplot(122, projection='3d')
        ax_data.set_title('Raw Data', fontsize=11)
        ax_fit.set_title('Gaussian Fits', fontsize=11)

        text_lines = []
        mode = 'SUM' if self.det_use_sum else 'AVG'

        for seg_id, (popt, roi_data) in fit_results.items():
            label = self.LABELS[seg_id]
            col = seg_colors_3d[seg_id]
            amp, x0, y0, sx, sy, theta, offset = popt
            peak_z = amp + offset

            # ── Left panel: raw data as wireframe ──────────────────
            ax_data.plot_wireframe(
                xx, yy, roi_data, color=col, alpha=0.35,
                rstride=2, cstride=2, linewidth=0.5,
                label=label)

            # ── Right panel: fit surface ───────────────────────────
            z_fit = self._gaussian2d(
                (xx.astype(float), yy.astype(float)), *popt
            ).reshape(xx.shape)
            ax_fit.plot_surface(
                xx, yy, z_fit, color=col, alpha=0.3,
                rstride=2, cstride=2, edgecolor='none',
                label=label)

            # ── Vertical center line (from 0 up to peak) ──────────
            ax_fit.plot([x0, x0], [y0, y0], [offset, peak_z],
                        color=col, linewidth=2, linestyle='-',
                        alpha=0.9)
            # Small marker at the peak
            ax_fit.scatter([x0], [y0], [peak_z],
                           color=col, s=30, zorder=5)

            # ── Major-axis indicator at peak top ───────────────────
            # Line length = max(sx, sy), oriented along theta
            major_s = max(sx, sy)
            # If sy > sx, the major axis is at theta+90°
            if sy > sx:
                axis_angle = theta + np.pi / 2
            else:
                axis_angle = theta
            dx = major_s * np.cos(axis_angle)
            dy = major_s * np.sin(axis_angle)
            ax_fit.plot(
                [x0 - dx, x0 + dx],
                [y0 - dy, y0 + dy],
                [peak_z, peak_z],
                color=col, linewidth=3, alpha=0.8,
                linestyle='-')

            # Parameter text
            text_lines.append(
                f"{label}: A={amp:.1f}  x₀={x0:.1f}  y₀={y0:.1f}  "
                f"σx={sx:.2f}  σy={sy:.2f}  θ={np.degrees(theta):.1f}°  "
                f"off={offset:.1f}")

        # Axis labels for both panels
        for ax in (ax_data, ax_fit):
            ax.set_xlabel('X (px)')
            ax.set_ylabel('Y (px)')
            ax.set_zlabel('Intensity')

        ax_data.legend(fontsize=8, loc='upper left')

        fig.suptitle(f'2D Gaussian Fits ({mode})', fontsize=13, y=0.98)

        # Parameter annotation
        param_text = '\n'.join(text_lines)
        fig.text(0.02, 0.02, param_text, fontsize=8,
                 family='monospace',
                 bbox=dict(boxstyle='round', fc='white', alpha=0.85))

        # ── View angle controls ────────────────────────────────────
        self._fit_fig = fig
        self._fit_axes = (ax_data, ax_fit)

        fig.text(0.72, 0.02, 'Elev:', fontsize=8, color='#555')
        ax_elev = fig.add_axes([0.76, 0.01, 0.06, 0.03])
        tb_elev = TextBox(ax_elev, '', initial='30')
        tb_elev.on_submit(lambda t: self._set_fit_view(elev=t))

        fig.text(0.83, 0.02, 'Azim:', fontsize=8, color='#555')
        ax_azim = fig.add_axes([0.87, 0.01, 0.06, 0.03])
        tb_azim = TextBox(ax_azim, '', initial='-60')
        tb_azim.on_submit(lambda t: self._set_fit_view(azim=t))

        fig.text(0.94, 0.02, 'Roll:', fontsize=8, color='#555')
        ax_roll = fig.add_axes([0.96, 0.01, 0.035, 0.03])
        tb_roll = TextBox(ax_roll, '', initial='0')
        tb_roll.on_submit(lambda t: self._set_fit_view(roll=t))

        # Store references to keep widgets alive
        self._fit_tb = (tb_elev, tb_azim, tb_roll)
        self._fit_syncing = False  # guard against recursive sync

        fig.subplots_adjust(left=0.01, right=0.99, top=0.93,
                            bottom=0.10, wspace=0.02)

        # Link camera between the two 3D panels
        fig.canvas.mpl_connect('draw_event', self._sync_fit_views)

    def _sync_fit_views(self, event=None):
        """Sync view angles between the two 3D axes on every draw."""
        if self._fit_syncing:
            return
        if self._fit_fig is None or not plt.fignum_exists(self._fit_fig.number):
            return
        ax_data, ax_fit = self._fit_axes

        # Detect which axis the user interacted with (the one that changed)
        # and copy its view to the other
        src, dst = ax_data, ax_fit
        if hasattr(self, '_last_fit_elev'):
            # Check which one changed
            if (ax_fit.elev != self._last_fit_elev or
                    ax_fit.azim != self._last_fit_azim):
                src, dst = ax_fit, ax_data

        self._fit_syncing = True
        dst.view_init(elev=src.elev, azim=src.azim, roll=src.roll)
        self._last_fit_elev = src.elev
        self._last_fit_azim = src.azim
        self._fit_syncing = False

    def _set_fit_view(self, elev=None, azim=None, roll=None):
        """Update view angle on both 3D axes from text box input."""
        if self._fit_fig is None or not plt.fignum_exists(self._fit_fig.number):
            return
        self._fit_syncing = True
        for ax in self._fit_axes:
            if elev is not None:
                try:
                    ax.elev = float(elev)
                except ValueError:
                    pass
            if azim is not None:
                try:
                    ax.azim = float(azim)
                except ValueError:
                    pass
            if roll is not None:
                try:
                    ax.roll = float(roll)
                except ValueError:
                    pass
        self._fit_syncing = False
        self._fit_fig.canvas.draw_idle()

    # ─── Status ──────────────────────────────────────────────────────
    def _update_status(self):
        if self.segments is None:
            self._status.set_text('')
            return
        nr = (self.segments == 1).sum()
        nb = (self.segments == 2).sum()
        ng = (self.segments == 3).sum()
        total = self.segments.size
        out = os.path.basename(self.output_dir) if self.output_dir else ''
        self._status.set_text(
            f'R:{nr}  B:{nb}  G:{ng}  free:{total-nr-nb-ng}\n'
            f'→ {out}/')

    # ─── Save ────────────────────────────────────────────────────────
    def _get_h5_indices(self):
        """Map segment pixels to H5 frame indices, accounting for sync offset.

        Returns a 1D array the same length as image_stack, where each
        element is the segment ID (0/1/2/3) for the corresponding H5 frame.
        When _h5_offset is +1, pixel (r,c) maps to H5 frame (r*ncols+c)+1.
        """
        n_fr = self.image_stack.shape[0]
        seg_flat = self.segments.ravel()

        if self._h5_offset == 0:
            return seg_flat[:n_fr]

        # Build offset mapping: pixel i → H5 frame i + offset
        n_px = len(seg_flat)
        h5_seg = np.zeros(n_fr, dtype=int)
        for i in range(n_px):
            h5_idx = i + self._h5_offset
            if 0 <= h5_idx < n_fr:
                h5_seg[h5_idx] = seg_flat[i]
        return h5_seg

    def _build_segment_info(self):
        """Build the segment_info dict (shared by save and auto-save)."""
        ncols = self.segments.shape[1]
        info = {
            'prefix': self.data_prefix,
            'source_npy': os.path.basename(self.npy_path),
            'scan_number': self.meta.get('scan_number'),
            'h5_file': self.meta.get('h5_file'),
            'mda_file': self.meta.get('mda_file'),
            'map_type': self.meta.get('map_type'),
            'dirty_fix': self.meta.get('dirty_fix', False),
            'h5_frame_offset': self._h5_offset,
            'map_shape': list(self.segments.shape),
            'detector_mode': 'sum' if self.det_use_sum else 'avg',
            'segments': {}
        }
        for seg_id, label in self.LABELS.items():
            rows, cols = np.where(self.segments == seg_id)
            frame_indices = (rows * ncols + cols + self._h5_offset).tolist()
            seg_info = {
                'pixel_count': int(len(frame_indices)),
                'frame_indices': frame_indices,
                'pixel_coords': list(zip(rows.tolist(), cols.tolist()))
            }
            if seg_id in self._last_fit_results:
                p = self._last_fit_results[seg_id]
                seg_info['gaussian_fit'] = {
                    'amplitude': p[0],
                    'x0': p[1],
                    'y0': p[2],
                    'sigma_x': p[3],
                    'sigma_y': p[4],
                    'theta_deg': float(np.degrees(p[5])),
                    'offset': p[6]
                }
            info['segments'][label.lower()] = seg_info
        return info

    def _save_json(self):
        """Save segment_info.json to the output directory."""
        if not self.output_dir or self.segments is None:
            return
        os.makedirs(self.output_dir, exist_ok=True)
        info = self._build_segment_info()
        path = os.path.join(self.output_dir, 'segment_info.json')
        with open(path, 'w') as f:
            json.dump(info, f, indent=2)
        print(f"  → segment_info.json saved")

    def _on_save(self, event):
        """Save all segment data to the auto-derived output folder."""
        if self.segments is None:
            print("Nothing to save — load a file first.")
            return

        os.makedirs(self.output_dir, exist_ok=True)

        # Masks
        for seg_id, label in self.LABELS.items():
            mask = self.segments == seg_id
            np.save(os.path.join(self.output_dir,
                                 f'mask_{label.lower()}.npy'), mask)

        # Full segment array
        np.save(os.path.join(self.output_dir, 'segment_map.npy'),
                self.segments)

        # JSON with pixel lists + fit params
        self._save_json()

        rel = os.path.join(os.path.basename(self.scan_folder),
                           os.path.basename(self.output_dir))
        print(f"\nSaved to {rel}/")
        print(f"  mask_red.npy, mask_blue.npy, mask_green.npy")
        print(f"  segment_map.npy")
        print(f"  segment_info.json")


# ─────────────────────────────────────────────────────────────────────
def main():
    npy_path = None
    if len(sys.argv) > 1:
        p = sys.argv[1]
        if os.path.isfile(p) and p.endswith('.npy'):
            npy_path = p
    SegmentPainter(npy_path)


if __name__ == '__main__':
    main()
