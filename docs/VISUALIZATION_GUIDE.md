# Murray's Law Bifurcation Visualization - User Guide

## Overview

The Murray Bifurcation Visualizer creates an interactive 3D visualization showing:
- **Semi-translucent vessel surface** extracted from binary 3D volume
- **Opaque skeleton** rendered as tubes inside the surface
- **Localized bifurcation coloring** using spheres around branch points
- **Gray coloring** everywhere else (using NaN)
- **Switchable metrics** for analyzing Murray's Law consistency

## Key Concepts

### Murray's Law

Murray's Law states that at vascular bifurcations, the cube of the parent radius equals the sum of cubes of the daughter radii:

```
r_parent³ = r_daughter1³ + r_daughter2³
```

### Metrics Available

1. **Murray Ratio (φ)** - Default metric
   ```
   φ = (r₁³ + r₂³) / rₚ³
   ```
   - φ = 1.0 means perfect Murray consistency
   - φ > 1.0 means daughters larger than expected
   - φ < 1.0 means daughters smaller than expected
   - **This is the recommended default metric**

2. **Murray Residual (ε)**
   ```
   ε = 1 - φ
   ```
   - ε = 0.0 means perfect Murray consistency
   - Useful for highlighting deviations

3. **Murray Exponent (γ)**
   - Estimates the best-fit exponent for each bifurcation
   - Solves: rₚ^γ = r₁^γ + r₂^γ
   - Only valid when rₚ > max(r₁, r₂)
   - Shows deviation from cube law (γ=3)

### Generalized Exponent

You can compute Murray ratio with any exponent:
```python
viz.compute_murray_metrics(gamma=2.5)  # Use γ=2.5 instead of 3.0
```

## Installation

Required packages:
```bash
pip install pyvista numpy scipy h5py scikit-image
```

## Basic Usage

### Simple Example

```python
from murray_bifurcation_visualizer import MurrayBifurcationVisualizer

# Create visualizer
viz = MurrayBifurcationVisualizer(
    skeleton_path="GS55_skeleton2.pkl",
    volume_path="path/to/volume.mat",
    volume_var='rescaled_vol',
    volume_spacing=(16, 16, 16)  # Physical units
)

# Detect bifurcations
viz.detect_bifurcations()

# Compute Murray metrics (default γ=3.0)
viz.compute_murray_metrics()

# Print summary statistics
viz.print_bifurcation_summary()

# Extract surface
viz.extract_surface()

# Assign metrics to surface
viz.assign_metrics_to_surface(metric_name='murray_phi')

# Create skeleton
viz.create_skeleton_mesh()

# Visualize
viz.visualize(metric_name='murray_phi')
```

### Switching Metrics

```python
# View Murray ratio (φ) - RECOMMENDED DEFAULT
viz.assign_metrics_to_surface('murray_phi')
viz.visualize(metric_name='murray_phi')

# View Murray residual (ε)
viz.assign_metrics_to_surface('murray_residual')
viz.visualize(metric_name='murray_residual')

# View estimated exponent (γ)
viz.assign_metrics_to_surface('murray_gamma')
viz.visualize(metric_name='murray_gamma')
```

### Custom Exponent

```python
# Compute with different exponent (e.g., γ=2.5)
viz.visualize(metric_name='murray_phi', gamma=2.5)
```

## Configuration Options

### Bifurcation Detection

```python
viz.config['local_radius_neighbors'] = 5  # More neighbors for smoother local radius
```

### Sphere Visualization

```python
viz.config['bifurcation_sphere_scale'] = 2.0  # Larger spheres around bifurcations
viz.config['kdtree_k_candidates'] = 24  # More candidates for surface assignment
```

### Visual Appearance

```python
viz.config['surface_opacity'] = 0.3  # More/less transparent surface
viz.config['nan_color'] = 'lightgray'  # Color for non-bifurcation regions
viz.config['skeleton_color'] = 'darkred'  # Change skeleton color
viz.config['colormap'] = 'coolwarm'  # Different colormap
```

### Color Maps

Recommended colormaps for different metrics:

- **murray_phi**: 'viridis' (diverging around 1.0)
- **murray_residual**: 'coolwarm' (diverging around 0.0)
- **murray_gamma**: 'plasma' (sequential for exponent values)

```python
viz.config['colormap'] = 'coolwarm'
```

## Understanding the Visualization

### Surface Coloring

- **Colored regions**: Near bifurcations (within sphere radius)
- **Gray regions**: Far from bifurcations (NaN values)
- **Sphere radius**: `scale × local_radius_estimate`
  - Default scale: 1.5
  - Ensures localized highlighting

### Color Interpretation

For **Murray Ratio (φ)**:
- **Green/Yellow (φ ≈ 1.0)**: Consistent with Murray's Law
- **Blue (φ < 1.0)**: Daughters thinner than expected
- **Red (φ > 1.0)**: Daughters thicker than expected

For **Murray Residual (ε)**:
- **White (ε ≈ 0.0)**: Consistent with Murray's Law
- **Blue (ε < 0.0)**: Same as φ < 1.0
- **Red (ε > 0.0)**: Same as φ > 1.0

For **Murray Exponent (γ)**:
- **Yellow (γ ≈ 3.0)**: Cube law
- **Purple (γ < 3.0)**: Sub-cubic relationship
- **Red (γ > 3.0)**: Super-cubic relationship

### Interactive Controls

**Mouse:**
- Left-drag: Rotate view
- Right-drag: Pan view
- Scroll: Zoom in/out

**Keyboard:**
- `q`: Quit visualization
- Click surface point: Print metric value at that location

## Volume File Format

### MATLAB Files

The visualizer supports both MATLAB formats:

**v7.2 and earlier** (scipy.io.loadmat):
```python
volume_path = "volume.mat"
volume_var = 'rescaled_vol'
```

**v7.3 HDF5** (h5py):
```python
# Automatically detected and handled
# Data is transposed to correct for MATLAB column-major order
```

### Binary Volume Requirements

- Values: 1 = vessel, 0 = background
- Any non-zero value interpreted as vessel
- Converted to boolean internally

### Coordinate System Alignment

The skeleton and volume must be in the same coordinate system:

```python
# If skeleton is in voxel coordinates, convert to physical:
skeleton_physical = skeleton_voxels * spacing

# Or ensure volume spacing matches skeleton units
volume_spacing = (16, 16, 16)  # e.g., microns
```

## Advanced Usage

### Batch Analysis

Analyze multiple metrics at once:

```python
metrics = ['murray_phi', 'murray_residual', 'murray_gamma']

for metric in metrics:
    viz.assign_metrics_to_surface(metric)
    viz.visualize(metric_name=metric)
```

### Export Statistics

```python
import pandas as pd

# Extract bifurcation data
data = []
for i, bif in enumerate(viz.bifurcations):
    row = {
        'bif_id': i,
        'x': bif['position'][0],
        'y': bif['position'][1],
        'z': bif['position'][2],
        'r_parent': bif['r_parent'],
        'r_daughter1': bif['r_daughter1'],
        'r_daughter2': bif['r_daughter2'],
    }

    # Add metrics
    for metric_name, values in viz.bifurcation_metrics.items():
        row[metric_name] = values[i]

    data.append(row)

df = pd.DataFrame(data)
df.to_csv('bifurcation_analysis.csv', index=False)
```

### Custom Metric

Define your own bifurcation metric:

```python
# Example: Area-based metric
def compute_area_ratio(bif):
    rp = bif['r_parent']
    r1 = bif['r_daughter1']
    r2 = bif['r_daughter2']

    area_p = np.pi * rp**2
    area_daughters = np.pi * (r1**2 + r2**2)

    return area_daughters / area_p

# Apply to all bifurcations
custom_metric = np.array([compute_area_ratio(b) for b in viz.bifurcations])
viz.bifurcation_metrics['area_ratio'] = custom_metric

# Visualize
viz.assign_metrics_to_surface('area_ratio')
viz.visualize(metric_name='area_ratio')
```

### Surface-Only Visualization

If you don't have a volume file:

```python
# Create visualizer without volume
viz = MurrayBifurcationVisualizer(
    skeleton_path="skeleton.pkl",
    volume_path=None  # No volume
)

# Skip surface extraction
viz.surface_mesh = None

# Visualize skeleton with bifurcation markers
# (requires custom implementation)
```

## Troubleshooting

### No bifurcations detected

**Possible causes:**
- Skeleton has no degree-3 vertices
- Check: `print(np.bincount(degrees))`

**Solutions:**
- Verify skeleton has branching structure
- Check if skeleton was over-simplified

### Surface and skeleton misaligned

**Possible causes:**
- Coordinate system mismatch
- Wrong volume spacing

**Solutions:**
```python
# Check skeleton bounds
print(viz.skeleton.vertices.min(axis=0))
print(viz.skeleton.vertices.max(axis=0))

# Check volume bounds in physical space
vol_bounds = np.array(viz.volume_data.shape) * viz.volume_spacing
print(vol_bounds)

# Adjust spacing if needed
viz.volume_spacing = (new_sx, new_sy, new_sz)
```

### Few surface points colored

**Possible causes:**
- Bifurcation spheres too small
- Bifurcations far from surface

**Solutions:**
```python
# Increase sphere scale
viz.config['bifurcation_sphere_scale'] = 2.5

# Increase KD-tree candidates
viz.config['kdtree_k_candidates'] = 32

# Re-assign metrics
viz.assign_metrics_to_surface('murray_phi')
```

### NaN values for Murray exponent

**Expected behavior:**
- γ estimation requires rₚ > max(r₁, r₂)
- Invalid when parent thinner than daughters
- Optimization may fail to converge

**Not a bug** - these bifurcations don't satisfy Murray's Law assumptions

## Performance Tips

### Large Volumes

For large volumes (>500³ voxels):

```python
# Use flying_edges (faster than marching_cubes)
viz.extract_surface(method='flying_edges')

# Reduce surface resolution
viz.surface_mesh = viz.surface_mesh.decimate(0.5)  # Keep 50%
```

### Many Bifurcations

For skeletons with >10,000 bifurcations:

```python
# Increase KD-tree candidates for better accuracy
viz.config['kdtree_k_candidates'] = 32

# Or reduce for speed
viz.config['kdtree_k_candidates'] = 8
```

## Validation

### Sanity Checks

```python
# Check metric consistency
phi = viz.bifurcation_metrics['murray_phi']
residual = viz.bifurcation_metrics['murray_residual']

# Should be true: residual ≈ 1 - phi
print(np.allclose(residual, 1 - phi, rtol=1e-5))

# Check valid bifurcations
valid_bifs = ~np.isnan(phi)
print(f"Valid bifurcations: {np.sum(valid_bifs)} / {len(phi)}")
```

### Compare Metrics

```python
# Scatter plot of φ vs γ
import matplotlib.pyplot as plt

phi = viz.bifurcation_metrics['murray_phi']
gamma = viz.bifurcation_metrics['murray_gamma']

valid = ~np.isnan(gamma)

plt.scatter(gamma[valid], phi[valid], alpha=0.5)
plt.xlabel('Estimated Exponent γ')
plt.ylabel('Murray Ratio φ (γ=3)')
plt.axhline(1.0, color='r', linestyle='--', label='Perfect φ')
plt.axvline(3.0, color='b', linestyle='--', label='Cube law γ')
plt.legend()
plt.show()
```

## Example Output

### Terminal Output

```
Murray's Law Bifurcation Visualizer
============================================================

Initializing visualizer...
Loading skeleton from GS55_skeleton2.pkl...
  Loaded skeleton 1
  Vertices: 45,832
  Edges: 45,831
  Radius range: 8.45 to 287.32

Loading volume from volume.mat...
  Loaded via h5py (v7.3) with transpose
  Volume shape: (512, 512, 256)
  Non-zero voxels: 8,432,156

Detecting bifurcations...
  Found 3,421 degree-3 bifurcations
  Analyzed 3,421 bifurcations
  Parent radius range: 12.34 to 245.67

Computing Murray metrics (gamma=3.0)...
  Murray ratio φ (γ=3.0):
    Range: 0.456 to 1.834
    Mean: 1.023, Median: 0.987
    Near 1.0 (consistent): 2,156 / 3,421
  Murray exponent γ:
    Valid estimates: 2,834 / 3,421
    Range: 2.123 to 4.567
    Mean: 2.987, Median: 2.945

Extracting surface using flying_edges...
  Contouring at iso-value 0.5...
  Extracted surface: 1,245,678 points, 2,491,234 cells

Assigning metric 'murray_phi' to surface via sphere transfer...
  Surface points: 1,245,678
  Bifurcations: 3,421
  Sphere radius range: 18.51 to 368.51
  Querying 16 nearest bifurcations per surface point...
  Assigned metric to 234,567 / 1,245,678 surface points (18.8%)

Creating skeleton mesh...
  Created skeleton with 45,831 edges as tubes

Creating visualization for metric 'murray_phi'...
  Color limits: 0.678 to 1.432

Interactive controls:
  Mouse: Left-drag=rotate, Right-drag=pan, Scroll=zoom
  Click surface point to see metric value
  Press 'q' to quit
```

## Citation

If you use this visualization in your research, please cite:

Murray CD (1926) The physiological principle of minimum work: I. The vascular system and the cost of blood volume. PNAS 12(3):207-214.

## Support

For issues or questions:
- Check this guide first
- Review the code comments in `murray_bifurcation_visualizer.py`
- Ensure all dependencies are correctly installed
- Verify coordinate system alignment between skeleton and volume
