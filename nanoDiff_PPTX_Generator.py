import os
import json
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt

# --- Configuration ---
ANALYSIS_ROOT = "/Users/smith.scott/Documents/Postdoc/Data/S26_20260211/analysis"

# Specify which scans to include in the presentation
SCAN_LIST = list(range(78, 150))

OUTPUT_PPTX = "Nanodiffraction_Summary.pptx"

# Override list — scans in this set will ALWAYS use the raster (non-boustrophedon)
# template, regardless of what the metadata.json says.
# Leave empty if no overrides are needed.
OVERRIDE_RASTER_SCANS = {122,123,124,125,126}  # e.g. {89, 126, 140}

    # 85,86,87,88,92,98,100,
    # 103,104,105,106,107,108,109,110,
    # 112,113,114,115,116,117,118,119,120,
    # 121,122,123,124,125,

# Section heading slides — insert a title slide before each scan's ROI slides
ENABLE_SECTION_SLIDES = True

def add_slide_header(slide, scan_num, roi_name, detector_path, is_boustrophedon):
    """Add the standardized title and detector image to the slide.

    Template selection is driven by the boustrophedon_correction flag
    from the scan's metadata.json:
        True  → Template 1 (snake/boustrophedon scans)
        False → Template 2 (raster scans)
    If the flag is missing (unknown), falls back to a default layout.
    """
    # Add Title text box (Top Left)
    txBox = slide.shapes.add_textbox(Inches(0.5), Inches(0.2), Inches(6.0), Inches(1.0))
    tf = txBox.text_frame
    p = tf.paragraphs[0]
    p.text = f"Scan #{scan_num:04d} — {roi_name}"
    p.font.size = Pt(32)
    p.font.bold = True

    # Add Detector Image (Top Right) — position depends on template
    if detector_path.exists():
        if is_boustrophedon is True:
            # Template 1: boustrophedon/snake scans
            slide.shapes.add_picture(str(detector_path), left=Inches(7.30), top=Inches(0.00), width=Inches(5.45))
        elif is_boustrophedon is False:
            # Template 2: raster (non-boustrophedon) scans
            slide.shapes.add_picture(str(detector_path), left=Inches(8.73), top=Inches(0.00), height=Inches(3.5))
        else:
            # Default fallback (unknown / key missing)
            slide.shapes.add_picture(str(detector_path), left=Inches(7.30), top=Inches(0.00), width=Inches(5.45))


def add_core_images(slide, inten_img, xcom_img, ycom_img, is_boustrophedon):
    """Place the 3 main analysis images using the appropriate template."""
    if is_boustrophedon is True:
        # Template 1: boustrophedon/snake scans
        slide.shapes.add_picture(str(inten_img), left=Inches(0.2), top=Inches(1.2), width=Inches(6.3))
        slide.shapes.add_picture(str(xcom_img),  left=Inches(0.2), top=Inches(4.42), width=Inches(6.3))
        slide.shapes.add_picture(str(ycom_img),  left=Inches(6.83), top=Inches(4.42), width=Inches(6.3))
    else: #is_boustrophedon is False:
        # Template 2: raster (non-boustrophedon) scans
        slide.shapes.add_picture(str(ycom_img),  left=Inches(8.14), top=Inches(3.5), height=Inches(4.00))
        slide.shapes.add_picture(str(xcom_img),  left=Inches(4.07), top=Inches(2.18), height=Inches(4.00))
        slide.shapes.add_picture(str(inten_img), left=Inches(0), top=Inches(0.92), height=Inches(4.00))
    # else:
    #     # Default fallback (unknown / key missing)
        # slide.shapes.add_picture(str(inten_img), left=Inches(0.2), top=Inches(1.2), height=Inches(2.9))
        # slide.shapes.add_picture(str(xcom_img),  left=Inches(0.2), top=Inches(4.42), height=Inches(2.9))
        # slide.shapes.add_picture(str(ycom_img),  left=Inches(6.83), top=Inches(4.42), height=Inches(2.9))


def create_presentation():
    # Initialize presentation
    prs = Presentation()

    # Set to 16:9 wide format (13.333 x 7.5 inches)
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    # Use a blank slide layout (typically index 6)
    blank_slide_layout = prs.slide_layouts[6]

    print(f"Generating PowerPoint for Scans: {SCAN_LIST}")

    for scan_num in SCAN_LIST:
        scan_dir = Path(ANALYSIS_ROOT) / f"Scan_{scan_num:04d}"
        if not scan_dir.exists():
            print(f"  ⚠ Directory not found for Scan {scan_num}, skipping.")
            continue

        # Section heading slide
        if ENABLE_SECTION_SLIDES:
            section_slide = prs.slides.add_slide(blank_slide_layout)
            txBox = section_slide.shapes.add_textbox(
                Inches(0), Inches(2.5), prs.slide_width, Inches(2.5))
            tf = txBox.text_frame
            p = tf.paragraphs[0]
            p.text = f"Scan #{scan_num:04d}"
            p.font.size = Pt(54)
            p.font.bold = True
            p.alignment = 1  # PP_ALIGN.CENTER

        # Get all ROIs in the scan directory
        roi_dirs = [d for d in scan_dir.iterdir() if d.is_dir()]
        roi_dirs.sort()  # Sort alphabetically

        for roi_dir in roi_dirs:
            roi_name = roi_dir.name
            metadata_path = roi_dir / "metadata.json"

            if not metadata_path.exists():
                continue

            with open(metadata_path, 'r') as f:
                metadata = json.load(f)

            # Determine template from metadata — None means "unknown/default"
            is_boustrophedon = metadata.get("boustrophedon_correction", None)

            # Apply override if this scan is in the override list
            if scan_num in OVERRIDE_RASTER_SCANS:
                is_boustrophedon = False

            # File paths
            detector_img = roi_dir / "detector_roi.png"
            inten_img = roi_dir / "intensity.png"
            xcom_img = roi_dir / "xcom.png"
            ycom_img = roi_dir / "ycom.png"

            # --- Slide 1: Main Analysis (Inten, XCoM, YCoM) ---
            if inten_img.exists() and xcom_img.exists() and ycom_img.exists():
                slide1 = prs.slides.add_slide(blank_slide_layout)
                add_slide_header(slide1, scan_num, roi_name, detector_img, is_boustrophedon)
                add_core_images(slide1, inten_img, xcom_img, ycom_img, is_boustrophedon)

                bost_str = "Bost" if is_boustrophedon else ("Raster" if is_boustrophedon is False else "Default")
                print(f"  Scan {scan_num:04d} / {roi_name} → {bost_str} template")

            # --- Slide 2: Gaussian Fits (if enabled) ---
            if metadata.get("gaussian_fit_enabled", False):
                # Check if Gaussian images exist
                gauss_amp = roi_dir / "gauss_amplitude.png"
                gauss_offset = roi_dir / "gauss_offset.png"
                gauss_mux = roi_dir / "gauss_mu_x.png"
                gauss_muy = roi_dir / "gauss_mu_y.png"
                gauss_sigx = roi_dir / "gauss_sigma_x.png"
                gauss_sigy = roi_dir / "gauss_sigma_y.png"

                if gauss_amp.exists() and gauss_mux.exists():
                    slide2 = prs.slides.add_slide(blank_slide_layout)
                    add_slide_header(slide2, scan_num, f"{roi_name} (Gaussian Fits)",
                                     detector_img, is_boustrophedon)

                    # 2x3 Grid for Gaussian images
                    row1_top = Inches(1.5)
                    row2_top = Inches(4.4)
                    g_width = Inches(4.2)

                    # Top Row: Amplitude, Mu X, Sigma X
                    slide2.shapes.add_picture(str(gauss_amp),  left=Inches(0.2), top=row1_top, width=g_width)
                    slide2.shapes.add_picture(str(gauss_mux),  left=Inches(4.6), top=row1_top, width=g_width)
                    slide2.shapes.add_picture(str(gauss_sigx), left=Inches(9.0), top=row1_top, width=g_width)

                    # Bottom Row: Offset, Mu Y, Sigma Y
                    slide2.shapes.add_picture(str(gauss_offset), left=Inches(0.2), top=row2_top, width=g_width)
                    slide2.shapes.add_picture(str(gauss_muy),    left=Inches(4.6), top=row2_top, width=g_width)
                    slide2.shapes.add_picture(str(gauss_sigy),   left=Inches(9.0), top=row2_top, width=g_width)

    # Save final presentation
    print(f"Saving presentation to {OUTPUT_PPTX}...")
    prs.save(OUTPUT_PPTX)
    print("Done!")

if __name__ == "__main__":
    create_presentation()
