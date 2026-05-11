# Dirty Fix & Data Alignment Guide

## Quick Reference: What Happens Automatically

When you load a `.npy` file in Segment Painter, it reads the metadata and applies corrections automatically:

| You load... | Map fixed on load? | H5 offset | Notes |
|---|---|---|---|
| **User ROI** (dirty_fix=true) | ✅ Yes — shifted to align | +1 | Map shifted to match motor positions |
| **MDA Eiger channel** (dirty_fix=true) | No (already fixed by DUDE) | +1 | Already correct |
| **MDA fluorescence** (dirty_fix=true) | No (no bug for this detector) | 0 | No sync bug for non-Eiger |
| **Any map** (dirty_fix=false) | No | 0 | No bug present |

**You don't need to do anything special** — the tool handles it.

---

## How to Reuse Segmentation Masks Across Data Types

### Safe workflow:
1. **Load any map** (ROI or MDA) for a given scan
2. **Paint your segments** or **Load Segments** from a previous session
3. **Save** — the tool records the correct H5 frame indices for that data type
4. **Load a different map** for the same scan
5. **Load the same segments** — the mask applies to the same pixel grid
6. The tool re-detects the data type and adjusts H5 frame indices automatically

### Why this works:
After correction, **all maps from the same scan are aligned to the same pixel grid** — features appear at the same (row, col) position regardless of whether the source was ROI or MDA. The H5 frame offset is tracked separately.

---

## Pitfalls to Avoid

### 1. Don't mix dirty_fix ON and OFF exports for the same scan
If you export the same scan twice from DUDE — once with dirty_fix ON, once OFF — the MDA Eiger maps will be shifted relative to each other. Stick to one setting per scan.

### 2. Fluorescence maps have NO H5 detector correspondence
`FluoPt` (XSP3 channels) comes from the MDA file, not the Eiger H5. When you paint segments on a fluorescence map and click "Show Detector", you get the Eiger frames at those same scan positions — but the fluorescence and diffraction signals come from different detectors and may not correlate spatially in the same way.

### 3. The last pixel in a dirty-fixed map is stale
The dirty fix shifts data left by 1: `data[:-1] = data[1:]`. The very last pixel in the raster retains the value from before the shift. This is usually the bottom-right pixel. It's a known limitation — don't rely on that single pixel for analysis.

### 4. Check `segment_info.json` for audit trail
Every save records:
```json
{
  "dirty_fix": true,
  "h5_frame_offset": 1,
  "map_type": "Intensity (User ROI)"
}
```
If you ever need to verify which correction was applied, this is your audit trail.

### 5. Single-pixel features matter
With `h5_offset = 1`, painting pixel `(r, c)` retrieves H5 frame `r*ncols + c + 1`. For a 1-pixel-wide feature, this offset is the difference between the right frame and the wrong one. The tool handles this correctly now, but if you ever manually index into the H5 file outside this tool, remember to apply the +1 offset when dirty_fix was active.

---

## Technical Details

### The Eiger sync bug
During scanning at APS 26-ID, the Eiger detector records each frame **one trigger ahead** of the MDA positioner data. This means:
- MDA position index `i` → Eiger actually recorded the image for position `i-1`
- H5 frame `i` contains the detector image for motor position `i-1`

### The dirty fix correction
```python
ydata_flat[:-1] = ydata_flat[1:]  # shift data left by 1
```
After this, MDA index `i` shows the value from raw index `i+1`, which IS the correct data for position `i`.

### Segment Painter's correction
1. **ROI maps**: applies the same shift on load, so ROI maps align with dirty-fixed MDA maps
2. **H5 frame lookup**: uses `frame_index = pixel_index + 1` to retrieve the correct Eiger image
3. **Saved indices**: `frame_indices` in JSON already include the +1 offset
