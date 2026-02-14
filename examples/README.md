# Examples

This directory contains example scripts demonstrating how to use the Murray Bifurcation Visualization tool.

## Available Examples

### basic_example.py - Complete Workflow Demonstration

Comprehensive examples covering all major features of the visualizer.

**Functions included:**

1. **`basic_example()`** - Basic Murray ratio visualization
   - Load skeleton and volume data
   - Detect bifurcations automatically
   - Compute Murray metrics with default cube law (γ=3)
   - Visualize with bifurcation-localized coloring

2. **`metric_comparison()`** - Compare different Murray metrics
   - Murray ratio (φ) - compliance measure
   - Murray residual (ε) - deviation from ideal
   - Murray exponent (γ) - best-fit exponent per bifurcation
   - Side-by-side visualization of all three metrics

3. **`custom_exponent()`** - Test different Murray exponents
   - Compare γ=2.5, 3.0, 3.5
   - Analyze bifurcation compliance across exponent values
   - Print statistics for each exponent

4. **`custom_configuration()`** - Customize visualization appearance
   - Adjust bifurcation sphere sizes
   - Change surface opacity
   - Select different colormaps
   - Modify skeleton rendering color
   - Configure local radius estimation

5. **`export_statistics()`** - Export bifurcation data to CSV
   - Save all Murray metrics to CSV file
   - Include bifurcation positions and radii
   - Generate summary statistics
   - Useful for downstream analysis in R, Python, Excel

6. **`skeleton_only()`** - Visualize without volume data
   - Works with skeleton files only
   - Manual bifurcation marker creation
   - Custom PyVista visualization setup
   - Useful when volume data is unavailable

7. **`quick_analysis()`** - Batch processing without visualization
   - Compute metrics without rendering
   - Print bifurcation summary statistics
   - Categorize Murray compliance
   - Fast analysis for multiple specimens

## Running Examples

### Setup

1. **Place data files** in `../data/` directory:
   ```bash
   data/
   ├── chicken_liver_skeleton.pkl
   ├── GS40_liver_skeleton.pkl
   ├── GS55_skeleton2.pkl
   └── ...
   ```

2. **Update paths** in the example file:
   ```python
   skeleton_path="../data/your_skeleton.pkl"
   volume_path="../data/your_volume.mat"  # Optional
   ```

3. **Install dependencies**:
   ```bash
   pip install -r ../requirements.txt
   ```

### Run Example

```bash
cd examples/
python basic_example.py
```

### Select Specific Function

Edit the `if __name__ == "__main__":` section at the bottom of `basic_example.py`:

```python
if __name__ == "__main__":
    # Uncomment the function you want to run:

    basic_example()           # ← Currently active
    # metric_comparison()
    # custom_exponent()
    # custom_configuration()
    # export_statistics()
    # skeleton_only()
    # quick_analysis()
```

## Example Outputs

### Basic Example
- Interactive 3D window with rotatable vessel surface
- Bifurcations colored by Murray ratio (φ)
- Console output: bifurcation summary statistics

### Export Statistics
- **File**: `bifurcation_murray_analysis.csv`
- **Columns**: bifurcation_id, x, y, z, r_parent, r_daughter1, r_daughter2, local_radius, murray_phi, murray_residual, murray_gamma
- **Rows**: One per bifurcation

### Quick Analysis
- Console output only:
  - Murray consistency categories (consistent, under-perfused, over-perfused)
  - Cube law compliance statistics
  - No visualization window

## Data Requirements

### Minimum Requirements
- **Skeleton file** (.pkl) with vertices, edges, and radii
- No volume file needed for skeleton-only mode

### Full Visualization Requirements
- **Skeleton file** (.pkl)
- **Volume file** (.mat) - binary vessel mask
- **Voxel spacing** - (x, y, z) tuple in microns

See [../data/README.md](../data/README.md) for format specifications.

## Customization Tips

### Change Colormap
```python
viz.config['colormap'] = 'viridis'  # or 'coolwarm', 'plasma', 'jet', etc.
```

### Adjust Bifurcation Sphere Size
```python
viz.config['bifurcation_sphere_scale'] = 3.0  # Larger spheres
```

### Change Opacity
```python
viz.config['surface_opacity'] = 0.1       # More transparent non-bifurcation regions
viz.config['colored_opacity'] = 0.9       # More opaque bifurcation regions
```

### Window Size
```python
viz.config['window_size'] = (1920, 1080)  # Full HD
```

## Troubleshooting

### ImportError: No module named 'murray_viz'
- Make sure you've installed the package or added src/ to PYTHONPATH
- The example file automatically adds `../src` to sys.path

### FileNotFoundError: skeleton file not found
- Check that your skeleton file is in `../data/`
- Verify the path in the example script matches your filename

### Volume loading fails
- Verify .mat file format (v7.2 or v7.3)
- Check variable name matches `volume_var` parameter
- For skeleton-only visualization, set `volume_path=None`

### PyVista rendering issues
- Ensure PyQt5 is installed: `pip install PyQt5`
- Try setting environment variable: `export QT_QPA_PLATFORM=offscreen`
- Check OpenGL support on your system

## Additional Resources

- **[Quick Start Guide](../docs/QUICKSTART.md)** - 5-minute introduction
- **[Visualization Guide](../docs/VISUALIZATION_GUIDE.md)** - Complete documentation
- **[Implementation Details](../docs/IMPLEMENTATION.md)** - Technical architecture
- **[Main README](../README.md)** - Project overview

## Creating Your Own Examples

Template for new analysis:

```python
#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from murray_viz import MurrayBifurcationVisualizer

# Your analysis code here
viz = MurrayBifurcationVisualizer(
    skeleton_path="../data/your_skeleton.pkl"
)

viz.detect_bifurcations()
viz.compute_murray_metrics()
# ... your custom analysis ...
```
