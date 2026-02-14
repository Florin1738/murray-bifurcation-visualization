# Murray Exponent Drivers Analysis - Installation & Usage Guide

## Overview

This script (`murray_exponent_drivers_analysis.py`) analyzes what factors drive variation in the fitted Murray exponent (alpha) at vascular bifurcations. It generates diagnostic plots to test whether exponent variation is due to vessel scale, branching geometry, or other factors.

---

## Installation: scikit-learn

The script requires **scikit-learn** for the permutation importance analysis (Plot 3). The script will run without it, but will skip the machine learning analysis.

### Safe Installation via Conda (Recommended)

1. **Open Anaconda Prompt** (search for it in Windows Start menu)

2. **Check your current environment:**
   ```bash
   conda info --envs
   ```
   The active environment has an asterisk (*) next to it.

3. **Activate your working environment** (if needed):
   ```bash
   conda activate your_environment_name
   ```

4. **Install scikit-learn:**
   ```bash
   conda install scikit-learn
   ```

   Or using conda-forge (often more up-to-date):
   ```bash
   conda install -c conda-forge scikit-learn
   ```

5. **Verify installation:**
   ```bash
   python -c "import sklearn; print('scikit-learn version:', sklearn.__version__)"
   ```

   You should see something like: `scikit-learn version: 1.3.0`

### Why Conda Instead of pip?

- Conda handles dependency compatibility automatically
- Ensures NumPy, SciPy are compatible versions
- Less likely to break your environment

---

## What the Script Does

### Generated Plots

The script generates **two sets** of plots:

#### 1. Combined Analysis (All Species Together)
Located in: `plots/`

- `murray_drivers_1_alpha_vs_radius.png` - Tests for scale-dependent effects
- `murray_drivers_2_alpha_vs_asymmetry.png` - Tests for geometry-dependent effects
- `murray_drivers_3_permutation_importance.png` - Quantifies feature importance
- `murray_drivers_summary_3panel.png` - Three-panel summary figure

#### 2. Individual Sample Analysis
Located in: `plots/individual_samples/[Sample_Name]/`

For each sample (Alligator, GS40_Monkey, GS55_Monkey, Chicken, Turtle):
- `[Sample]_1_alpha_vs_radius.png`
- `[Sample]_2_alpha_vs_asymmetry.png`
- `[Sample]_3_permutation_importance.png` (only if ≥100 bifurcations)
- `[Sample]_summary_3panel.png`

---

## Usage

### Running the Script

**Option 1: From Command Line**
```bash
cd "c:\Users\Florin\OneDrive - Johns Hopkins\Documents\kimimaro_skeletonization"
python murray_exponent_drivers_analysis.py
```

**Option 2: In Python/IPython**
```python
exec(open('murray_exponent_drivers_analysis.py').read())
```

**Option 3: In Jupyter Notebook**
```python
%run murray_exponent_drivers_analysis.py
```

### Expected Output

The script will:
1. Load all 5 skeleton samples
2. Extract bifurcation features (radii, angles, tortuosity, etc.)
3. Fit local Murray exponent (alpha) for each bifurcation
4. Generate combined diagnostic plots
5. Generate individual plots for each sample

**Progress indicators:**
- Real-time console output showing processing status
- File sizes and confirmation messages when plots are saved

**Runtime:** Approximately 2-5 minutes depending on:
- Number of bifurcations
- Whether scikit-learn is available (permutation importance takes longer)
- Your CPU speed

---

## Configuration

You can modify parameters at the top of the script:

### Analysis Parameters
```python
JUNCTION_RADIUS_POINTS = 5    # Points to average for local radius
MIN_PARENT_RADIUS = 2.0        # Minimum parent radius (filter noise)
ALPHA_BRACKET = (0.5, 6.0)     # Search range for alpha fitting
BRANCH_DIRECTION_POINTS = 3    # Points for estimating direction
TORTUOSITY_POINTS = 10         # Points for tortuosity calculation
```

### Output Settings
```python
OUTPUT_DIR = r'path\to\output\directory'
FIGURE_DPI = 300               # Resolution (300 for publication)
FIGURE_FORMAT = 'png'          # 'png', 'pdf', or 'svg'
```

### Sample Files
Edit `SKELETON_FILES` dictionary to add/remove samples:
```python
SKELETON_FILES = {
    'Sample Name': r'path\to\skeleton.pkl',
    # ... add more samples
}
```

---

## Interpreting the Plots

### Plot 1: Alpha vs Parent Radius (Hexbin)
**Question:** Is exponent variation scale-dependent?

**What to look for:**
- **Horizontal trend** → Alpha is independent of vessel size (good!)
- **Upward slope** → Larger vessels have higher alpha
- **Downward slope** → Smaller vessels have higher alpha (may indicate resolution limits)
- **Wide spread at small radii** → Measurement noise dominates at small scales

### Plot 2: Alpha vs Asymmetry (Hexbin)
**Question:** Is exponent variation driven by branching geometry?

**Asymmetry index:** a = min(r₁, r₂) / max(r₁, r₂)
- a ≈ 1 = symmetric split (equal daughters)
- a ≈ 0 = highly asymmetric split (one large, one small)

**What to look for:**
- **Horizontal trend** → Alpha independent of asymmetry
- **Trend with asymmetry** → Different branching modes have different exponents
- **Clustering** → Discrete branching types (trunk splits vs side branches)

### Plot 3: Permutation Importance (Bar Chart)
**Question:** Which features actually predict alpha?

**What to look for:**
- **High importance of log(r_parent)** → Scale-dependent
- **High importance of asymmetry** → Geometry-dependent
- **High importance of angle** → Hemodynamic optimization
- **Low overall R²** → Alpha is noisy or driven by unmeasured factors

---

## Troubleshooting

### "scikit-learn not available"
- Install scikit-learn using instructions above
- Script will still run and generate Plots 1, 2, and 4

### "Insufficient data for regression"
- Need ≥100 bifurcations for reliable permutation importance
- Check that skeletons loaded correctly
- Try reducing `MIN_PARENT_RADIUS` to include more junctions

### "No valid data for [sample]"
- Check skeleton file path is correct
- Verify skeleton contains bifurcations (degree-3 junctions)
- Ensure radius data exists in skeleton

### Plots look empty
- Check if `MIN_PARENT_RADIUS` is too high (excluding too many junctions)
- Verify voxel spacing is correct
- Check console output for number of valid bifurcations

### Memory errors
- Reduce `FIGURE_DPI` from 300 to 150
- Process samples one at a time instead of all together

---

## Output Structure

```
plots/
├── murray_drivers_1_alpha_vs_radius.png
├── murray_drivers_2_alpha_vs_asymmetry.png
├── murray_drivers_3_permutation_importance.png
├── murray_drivers_summary_3panel.png
└── individual_samples/
    ├── Alligator/
    │   ├── Alligator_1_alpha_vs_radius.png
    │   ├── Alligator_2_alpha_vs_asymmetry.png
    │   ├── Alligator_3_permutation_importance.png
    │   └── Alligator_summary_3panel.png
    ├── GS40_Monkey/
    │   └── ...
    ├── GS55_Monkey/
    │   └── ...
    ├── Chicken/
    │   └── ...
    └── Turtle/
        └── ...
```

---

## Mathematical Details

### Murray Exponent Fitting
For each bifurcation, solves:
```
r_parent^α = r_child1^α + r_child2^α
```

Uses log-space root finding for numerical stability:
```
f(α) = log(r_child1^α + r_child2^α) - α·log(r_parent) = 0
```

Solved using Brent's method (bracketed root finder).

### Feature Engineering
- **log(r_parent):** Vessel scale (log-transformed for scale invariance)
- **asymmetry:** min(r₁, r₂) / max(r₁, r₂) ∈ (0, 1]
- **branch_angle:** Angle between daughter direction vectors (degrees)
- **tortuosity:** path_length / straight_line_distance
- **radius_taper:** (r_near - r_far) / distance along branch

### Permutation Importance
- Trains Random Forest regressor on 70% of data
- Tests on held-out 30%
- Shuffles each feature and measures drop in R²
- Repeats 20-30 times for stability estimation

---

## Citation

If you use this analysis in publications, consider citing:
- Murray, C. D. (1926). The physiological principle of minimum work. PNAS.
- Sherman, T. F. (1981). On connecting large vessels to small. J Gen Physiol.

---

## Contact

For questions about the script, check the inline comments or refer to:
- `MURRAY_EXPONENT_FITTING_NOTES.md` - Mathematical background
- `KIMIMARO_DATA_STRUCTURES_GUIDE.md` - Data structure reference

---

**Last Updated:** 2026-01-22
