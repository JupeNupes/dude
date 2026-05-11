# 2D Gaussian Fit — Parameter Reference

Reference for interpreting the rotated 2D Gaussian fit used in `segment_painter.py` for analyzing nano-diffraction peak profiles.

---

## The Model

A **rotated 2D Gaussian** with 7 free parameters:

```
f(x, y) = A · exp(−[a(x − x₀)² + 2b(x − x₀)(y − y₀) + c(y − y₀)²]) + offset
```

where the rotation coefficients are:

```
a = cos²θ / (2σx²) + sin²θ / (2σy²)
b = −sin(2θ) / (4σx²) + sin(2θ) / (4σy²)
c = sin²θ / (2σx²) + cos²θ / (2σy²)
```

This parameterization allows the Gaussian ellipse to be oriented at an arbitrary angle `θ` relative to the detector axes.

---

## Parameter Definitions

| Parameter | Symbol | Units | Description |
|-----------|--------|-------|-------------|
| `amplitude` | A | counts | Height of the peak above the background |
| `x0` | x₀ | pixels | Horizontal center position on the detector |
| `y0` | y₀ | pixels | Vertical center position on the detector |
| `sigma_x` | σx | pixels | Standard deviation along the first principal axis |
| `sigma_y` | σy | pixels | Standard deviation along the second principal axis |
| `theta_deg` | θ | degrees | Rotation angle of the ellipse (CCW from x-axis) |
| `offset` | bg | counts | Constant background level |

---

## Derived Quantities

### Full Width at Half Maximum (FWHM)

The FWHM is the most commonly reported measure of peak width:

```
FWHM = 2√(2 ln 2) × σ  ≈  2.3548 × σ
```

So `FWHMx ≈ 2.355 × σx` and `FWHMy ≈ 2.355 × σy`.

### Integrated Intensity

The total integrated intensity (volume under the Gaussian, excluding background):

```
I_total = 2π × A × σx × σy
```

This is proportional to the total diffracted signal from that segment — more physically meaningful than `A` alone since it accounts for peak broadening. A broad, weak peak and a sharp, strong peak can have the same integrated intensity.

### Eccentricity (Peak Shape)

Ratio of the principal widths:

```
ε = σ_max / σ_min
```

- **ε ≈ 1** → circular peak (isotropic broadening)
- **ε > 1** → elongated peak (anisotropic broadening)
- **ε ≫ 1** → streak-like (strong preferred orientation or strain gradient)

### Aspect Ratio Convention

By convention, `σx` is the width along the first principal axis before rotation. After rotation by `θ`, the major axis of the ellipse points at angle `θ` from the horizontal. If `σx > σy`, the elongation direction is along `θ`. If `σx < σy`, it's along `θ + 90°`.

---

## Physical Interpretation for Diffraction

### Peak Position (x₀, y₀)

The center of the Bragg peak on the detector. Combined with the experimental geometry (sample-detector distance, detector tilt, beam center), this maps to the scattering vector **Q** and ultimately to:

- **d-spacing**: `d = 2π / |Q|` — the interplanar spacing of the diffracting planes
- **Bragg angle**: `2θ_B = 2 arcsin(λ / 2d)` where λ is the X-ray wavelength
- **Lattice parameter**: through the crystal structure's relationship between d-spacing and Miller indices (hkl)

**When comparing segments**: A shift in x₀ or y₀ between segments indicates a **strain gradient** across the sample — different regions have slightly different d-spacings, meaning the lattice is more compressed or expanded in one area vs. another.

### Peak Width (σx, σy)

Peak broadening encodes microstructural information via the **Scherrer equation** and its generalizations:

```
β = Kλ / (D cos θ)
```

where β is the FWHM (in radians of 2θ), K ≈ 0.9 is the shape factor, and D is the crystallite/domain size.

Sources of broadening (they add in quadrature):

| Source | Effect on σ | Signature |
|--------|-------------|-----------|
| **Small domain size** | Increases σ isotropically | σx ≈ σy, both large |
| **Microstrain** | Increases σ, often anisotropically | σx ≠ σy |
| **Mosaicity** | Increases σ along specific directions | Elongated peak |
| **Instrumental** | Baseline broadening (always present) | Minimum σ even for perfect crystal |
| **Stacking faults** | Asymmetric broadening | Peak may be non-Gaussian |

**When comparing segments**: Broader peaks (larger σ) in one region indicate smaller crystalline domains, more strain disorder, or greater mosaicity. If one segment shows significantly broader peaks than another, that region of the sample is more disordered.

### Rotation Angle (θ)

The orientation of the broadening ellipse on the detector. This relates to the **anisotropy direction** of the microstructural broadening:

- **θ aligned with radial direction** (away from beam center) → radial broadening → typically **strain broadening** (d-spacing variation)
- **θ aligned with azimuthal direction** (tangent to Debye ring) → azimuthal broadening → typically **mosaicity** or **domain misorientation**
- **θ at intermediate angle** → mixed contribution or a tilted strain field

**When comparing segments**: A change in θ between regions may indicate:
- Different crystallographic texture / preferred orientation
- Different strain field geometry (e.g., near a grain boundary vs. grain interior)
- Rotation of the crystal lattice itself

### Amplitude (A)

The peak height above background. In isolation, it depends on:

- Structure factor of the reflection (fixed for a given hkl)
- Illuminated volume of crystalline material
- Texture / orientation distribution
- Absorption

**When comparing segments**: Differences in amplitude (at the same hkl) suggest differences in **diffracting volume** or **crystallographic orientation**. Higher amplitude = more material oriented to satisfy the Bragg condition for this reflection.

### Background Offset

The constant pedestal under the peak. Contributions include:

- **Thermal diffuse scattering** (TDS) — phonon-related
- **Compton scattering** — incoherent
- **Fluorescence** — element-specific X-ray emission
- **Detector dark current / readout noise**
- **Air scatter / parasitic scattering**

A higher offset in one segment could indicate more fluorescence (if measuring near an absorption edge) or more diffuse scattering (more disorder).

---

## Comparing Segments: What to Look For

| Observation | Likely Physical Cause |
|-------------|----------------------|
| Δx₀ or Δy₀ between segments | Strain gradient (lattice parameter variation) |
| Larger σ in one segment | Smaller domains or more disorder in that region |
| σx ≠ σy (elongated peak) | Anisotropic strain, mosaicity, or preferred orientation |
| θ differs between segments | Different texture or strain field orientation |
| Higher A in one segment | More diffracting volume or better orientation alignment |
| Higher offset in one segment | More fluorescence, diffuse scattering, or disorder |
| Same x₀/y₀ but different σ | Same phase, different domain size / quality |
| Different x₀/y₀, similar σ | Strained but similarly ordered regions |

---

## Caveats & Limitations

1. **Gaussian is an approximation.** Real diffraction peaks may be better described by:
   - **Lorentzian** (Cauchy) — for strain-dominated broadening
   - **Voigt** or **Pseudo-Voigt** — mixed Gaussian + Lorentzian (most physical)
   - **Pearson VII** — flexible shape parameter
   
   If the Gaussian fit residuals are large or systematic, consider whether a different profile function is needed.

2. **Single-peak assumption.** The fit assumes one peak in the ROI. Overlapping reflections, split peaks (due to phase coexistence or twinning), or satellite reflections will produce unreliable fits.

3. **Pixel discretization.** For very sharp peaks (σ < 2 pixels), the fit becomes unreliable because the peak is under-sampled. The FWHM should span at least ~5 pixels for a reliable fit.

4. **Sum vs. Average mode.** Summed images have better statistics (higher counts) but the amplitude scales with the number of pixels in the segment. Average mode normalizes this out, making amplitudes directly comparable between segments of different sizes.

5. **Background subtraction.** The flat `offset` parameter assumes spatially uniform background across the ROI. If the background has a gradient (e.g., near a beam stop shadow), the fit may be biased. Consider fitting with a tilted plane background for more accuracy (not yet implemented).

6. **Correlation between parameters.** Amplitude and offset are correlated — a higher offset can suppress the fitted amplitude. Similarly, σ and A are correlated: broadening reduces peak height at constant integrated intensity.

---

## Converting Pixel Units to Physical Units

To convert σ from pixels to reciprocal space (Å⁻¹) or angular units (°):

```
σ_Q = σ_pixels × (pixel_size / sample_detector_distance) × (2π / λ)
```

Or more practically, if you have the Q-map or 2θ-map from the geometry calibration (the CSV files exported by DUDE: `gamma.csv`, `twotheta.csv`, etc.):

```python
# Load the 2θ map for the detector
twotheta = np.loadtxt('twotheta.csv', delimiter=',')

# Get the 2θ value at the peak center
twotheta_center = twotheta[int(y0), int(x0)]

# Get the 2θ gradient (degrees per pixel) at the peak
d2t_dx = np.gradient(twotheta, axis=1)[int(y0), int(x0)]
d2t_dy = np.gradient(twotheta, axis=0)[int(y0), int(x0)]

# Convert σ from pixels to degrees
sigma_2theta = sigma_pixels * sqrt(d2t_dx**2 + d2t_dy**2)
```
