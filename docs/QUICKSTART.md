# Murray Bifurcation Visualizer - Quick Start

## 30-Second Start

```python
from murray_bifurcation_visualizer import MurrayBifurcationVisualizer

viz = MurrayBifurcationVisualizer(
    skeleton_path="GS55_skeleton2.pkl",
    volume_path="volume.mat",
    volume_var='rescaled_vol',
    volume_spacing=(16, 16, 16)
)

viz.detect_bifurcations()
viz.compute_murray_metrics()
viz.extract_surface()
viz.create_skeleton_mesh()
viz.assign_metrics_to_surface('murray_phi')
viz.visualize(metric_name='murray_phi')
```

## What You'll See

![Expected Visualization]
- **Gray semi-transparent vessel surface** (most of it)
- **Colored patches around bifurcations** (Murray metric)
- **Black opaque skeleton** (centerline structure)
- **Interactive 3D viewer** (rotate, zoom, click)

## Files You Need

### Required
- [murray_bifurcation_visualizer.py](murray_bifurcation_visualizer.py) - Main code
- Your skeleton file (`.pkl`)
- Your volume file (`.mat`) - optional but recommended

### Documentation
- [MURRAY_QUICKSTART.md](MURRAY_QUICKSTART.md) - This file
- [MURRAY_VISUALIZATION_GUIDE.md](MURRAY_VISUALIZATION_GUIDE.md) - Complete guide
- [MURRAY_IMPLEMENTATION_SUMMARY.md](MURRAY_IMPLEMENTATION_SUMMARY.md) - Technical details

### Examples
- [murray_example.py](murray_example.py) - 7 example scripts

## Installation

```bash
pip install pyvista numpy scipy h5py scikit-image
```

**That's it!** No compilation, no complex setup.

## Three Usage Patterns

### 1. Quick Visualization (Default Script)

```bash
# Edit murray_bifurcation_visualizer.py main() with your paths
python murray_bifurcation_visualizer.py
```

### 2. Interactive Python

```python
from murray_bifurcation_visualizer import MurrayBifurcationVisualizer

viz = MurrayBifurcationVisualizer("skeleton.pkl", "volume.mat")
viz.detect_bifurcations()
viz.compute_murray_metrics()
viz.print_bifurcation_summary()  # See statistics
viz.extract_surface()
viz.create_skeleton_mesh()
viz.visualize(metric_name='murray_phi')
```

### 3. Example Scripts

```bash
python murray_example.py
```

Then uncomment the example you want in the `__main__` block.

## Understanding the Metrics

### Murray Ratio (φ) - RECOMMENDED

```
φ = (r₁³ + r₂³) / rₚ³
```

- **φ = 1.0** → Perfect Murray consistency (green/yellow)
- **φ < 1.0** → Daughters thinner than expected (blue)
- **φ > 1.0** → Daughters thicker than expected (red)

### Murray Residual (ε)

```
ε = 1 - φ
```

- **ε = 0.0** → Perfect Murray consistency
- Useful for highlighting deviations from ideal

### Murray Exponent (γ)

Estimates best-fit exponent per bifurcation:
```
rₚ^γ = r₁^γ + r₂^γ
```

- **γ ≈ 3.0** → Cube law (Murray's original)
- **γ < 3.0** → Sub-cubic relationship
- **γ > 3.0** → Super-cubic relationship

## Switching Metrics

```python
# View Murray ratio (default, recommended)
viz.visualize(metric_name='murray_phi')

# View Murray residual
viz.visualize(metric_name='murray_residual')

# View estimated exponent
viz.visualize(metric_name='murray_gamma')

# Custom exponent (e.g., square law)
viz.visualize(metric_name='murray_phi', gamma=2.0)
```

## Common Customizations

### Larger Bifurcation Spheres

```python
viz.config['bifurcation_sphere_scale'] = 2.5  # Default: 1.5
```

### More Transparent Surface

```python
viz.config['surface_opacity'] = 0.15  # Default: 0.25
```

### Different Colormap

```python
viz.config['colormap'] = 'coolwarm'  # Default: 'viridis'
# Options: viridis, plasma, inferno, coolwarm, seismic, turbo
```

### Change Skeleton Color

```python
viz.config['skeleton_color'] = 'darkred'  # Default: 'black'
```

## Interactive Controls

**Mouse:**
- **Left-drag**: Rotate view
- **Right-drag**: Pan view
- **Scroll**: Zoom in/out

**Keyboard:**
- **q**: Quit visualization

**Click:**
- **Click surface point**: Print metric value at that location

## Export Data

```python
import pandas as pd

# After running analysis
data = []
for i, bif in enumerate(viz.bifurcations):
    row = {
        'id': i,
        'x': bif['position'][0],
        'y': bif['position'][1],
        'z': bif['position'][2],
        'r_parent': bif['r_parent'],
        'r_daughter1': bif['r_daughter1'],
        'r_daughter2': bif['r_daughter2'],
        'murray_phi': viz.bifurcation_metrics['murray_phi'][i],
    }
    data.append(row)

df = pd.DataFrame(data)
df.to_csv('bifurcations.csv', index=False)
print(f"Exported {len(df)} bifurcations")
```

## Troubleshooting

### "No bifurcations detected"

**Cause:** Skeleton has no degree-3 vertices

**Fix:** Check your skeleton has branching:
```python
# Count vertex degrees
adjacency = defaultdict(list)
for v1, v2 in viz.skeleton.edges:
    adjacency[v1].append(v2)
    adjacency[v2].append(v1)

degrees = [len(adjacency[i]) for i in range(len(viz.skeleton.vertices))]
print(np.bincount(degrees))  # Should see degree-3 vertices
```

### Surface and skeleton misaligned

**Cause:** Coordinate system mismatch

**Fix:** Check bounds match:
```python
# Skeleton bounds
print("Skeleton:", viz.skeleton.vertices.min(axis=0), "to", viz.skeleton.vertices.max(axis=0))

# Volume bounds (in physical units)
vol_bounds = np.array(viz.volume_data.shape) * viz.volume_spacing
print("Volume:", [0,0,0], "to", vol_bounds)
```

Adjust `volume_spacing` if needed.

### Few colored regions

**Cause:** Bifurcation spheres too small or far from surface

**Fix:**
```python
# Increase sphere size
viz.config['bifurcation_sphere_scale'] = 2.0

# Increase KD-tree search
viz.config['kdtree_k_candidates'] = 24

# Re-assign
viz.assign_metrics_to_surface('murray_phi')
```

### Slow performance

**Cause:** Large volume or dense surface

**Fix:**
```python
# Use faster contouring
viz.extract_surface(method='flying_edges')

# Decimate surface
viz.surface_mesh = viz.surface_mesh.decimate(0.5)  # Keep 50%
```

## File Paths

### Update These in Your Script

```python
# In murray_bifurcation_visualizer.py main() or murray_example.py

skeleton_path = "GS55_skeleton2.pkl"  # ← Your skeleton file

volume_path = r"path\to\your\volume.mat"  # ← Your volume file

volume_var = 'rescaled_vol'  # ← Variable name in .mat file

volume_spacing = (16, 16, 16)  # ← Voxel spacing (sx, sy, sz)
```

### Default Paths in Code

The main script uses:
```python
skeleton_path = "GS55_skeleton2.pkl"
volume_path = r"\\10.162.80.11\Andre_kit\...\02_brain_smooth_4xy.mat"
```

**Change these to your actual file paths!**

## Workflow Summary

```
1. Load skeleton + volume
   ↓
2. Detect bifurcations (degree-3 junctions)
   ↓
3. Compute Murray metrics (φ, ε, γ)
   ↓
4. Extract vessel surface
   ↓
5. Transfer metrics to surface (sphere-based)
   ↓
6. Render (gray + colored bifurcations + skeleton)
```

## Next Steps

### Just Starting?
→ Run `python murray_example.py` with `quick_analysis()` uncommented
→ See statistics without waiting for visualization

### Want to Visualize?
→ Update paths in `main()` of `murray_bifurcation_visualizer.py`
→ Run `python murray_bifurcation_visualizer.py`

### Want to Customize?
→ Read [MURRAY_VISUALIZATION_GUIDE.md](MURRAY_VISUALIZATION_GUIDE.md)
→ Try examples in `murray_example.py`

### Want Technical Details?
→ Read [MURRAY_IMPLEMENTATION_SUMMARY.md](MURRAY_IMPLEMENTATION_SUMMARY.md)
→ Review code comments in `murray_bifurcation_visualizer.py`

## Key Insight

**Murray Ratio φ = 1.0 means consistency**, not ε = 0

- Use `murray_phi` as your default metric
- φ near 1.0 → bifurcation follows Murray's Law
- φ far from 1.0 → deviation from expected flow optimization

The visualization **only colors regions near bifurcations** (within spheres), leaving the rest gray. This makes deviations from Murray's Law immediately visible while providing anatomical context.

## Support

If something doesn't work:

1. Check file paths are correct
2. Verify dependencies installed: `pip list | grep pyvista`
3. Check skeleton has bifurcations: `viz.detect_bifurcations()`
4. Review troubleshooting section above
5. Read full guide: [MURRAY_VISUALIZATION_GUIDE.md](MURRAY_VISUALIZATION_GUIDE.md)

---

**That's it! You're ready to visualize Murray's Law at bifurcations.**

Happy analyzing! 🎨🔬
