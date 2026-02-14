# Murray Bifurcation 3D Visualization

Interactive 3D visualization of Murray's Law compliance at vascular bifurcations using PyVista.

## Overview

This tool analyzes and visualizes Murray's Law at vascular bifurcations, a fundamental principle in biological fluid dynamics that states:

**r_parent^γ = r_daughter1^γ + r_daughter2^γ**

where γ ≈ 3 (the "cube law") represents the optimal balance between blood volume and flow resistance.

## Features

- **3D Vessel Surface Rendering**: Semi-transparent volume visualization with PyVista
- **Bifurcation-Localized Coloring**: Murray metrics displayed only near branch points
- **Multiple Metrics**:
  - Murray ratio (φ) - compliance measure
  - Murray residual (ε) - deviation from ideal
  - Murray exponent (γ) - best-fit exponent per bifurcation
- **Interactive Visualization**: Rotate, zoom, pick points, switch metrics
- **Driver Analysis**: ML-based analysis of factors affecting Murray exponent variation
- **Multi-Species Support**: Comparative analysis across different specimens

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/murray-bifurcation-visualization.git
cd murray-bifurcation-visualization

# Install dependencies
pip install -r requirements.txt

# For advanced analysis features (optional)
pip install -r requirements-analysis.txt
pip install git+https://github.com/gcskoenig/fippy.git
```

### Basic Usage

```python
from murray_viz import MurrayBifurcationVisualizer

# Initialize visualizer
viz = MurrayBifurcationVisualizer(
    skeleton_path="data/your_skeleton.pkl",
    volume_path="data/your_volume.mat",
    volume_var='rescaled_vol',
    volume_spacing=(16, 16, 8)  # voxel spacing in microns
)

# Analysis pipeline
viz.detect_bifurcations()
viz.compute_murray_metrics(gamma=3.0)
viz.print_bifurcation_summary()

# Visualization
viz.extract_surface()
viz.create_skeleton_mesh()
viz.assign_metrics_to_surface('murray_phi')
viz.visualize(metric_name='murray_phi')
```

### Skeleton-Only Mode

```python
# Visualize without volume data
viz = MurrayBifurcationVisualizer(
    skeleton_path="data/your_skeleton.pkl",
    volume_path=None  # No volume required
)

viz.detect_bifurcations()
viz.compute_murray_metrics()
viz.create_skeleton_mesh()
# Custom visualization code (see examples/basic_example.py)
```

## Data Requirements

### Skeleton Files (Required)

Pickle files containing kimimaro skeleton objects with:
- `vertices`: Nx3 array of 3D coordinates
- `edges`: Mx2 array of vertex indices
- `radius` or `radii`: N-length array of vessel radii

Place skeleton files in the `data/` directory.

### Volume Files (Optional)

MATLAB `.mat` files with binary 3D volumes:
- Format: MATLAB v7.2 or v7.3
- Variable: typically `rescaled_vol`
- Content: Binary vessel mask
- Spacing: Known voxel dimensions

See [data/README.md](data/README.md) for detailed format specifications.

## Documentation

- **[Quick Start Guide](docs/QUICKSTART.md)** - Get started in 5 minutes
- **[Visualization Guide](docs/VISUALIZATION_GUIDE.md)** - Complete usage documentation
- **[Implementation Details](docs/IMPLEMENTATION.md)** - Technical architecture
- **[Exponent Fitting](docs/EXPONENT_FITTING.md)** - Murray exponent calculation
- **[Driver Analysis](docs/DRIVERS_ANALYSIS.md)** - ML-based feature importance
- **[Data Structures](docs/KIMIMARO_STRUCTURES.md)** - Kimimaro skeleton format

## Examples

See the [examples/](examples/) directory for:
- **basic_example.py** - Complete workflow demonstration
  - Basic visualization with default settings
  - Metric comparison (φ, ε, γ)
  - Custom exponent testing
  - Configuration options
  - Statistics export
  - Skeleton-only visualization

## Configuration

The visualizer is highly customizable:

```python
viz.config['bifurcation_sphere_scale'] = 2.5  # Sphere size around bifurcations
viz.config['surface_opacity'] = 0.25           # Non-bifurcation regions
viz.config['colored_opacity'] = 0.65           # Bifurcation regions
viz.config['colormap'] = 'coolwarm'            # Color scheme
viz.config['skeleton_color'] = 'black'         # Skeleton rendering
viz.config['window_size'] = (1400, 900)        # Window dimensions
```

## Dependencies

### Core Requirements

```
Python >= 3.9
pyvista >= 0.40.0
numpy >= 1.20.0
scipy >= 1.7.0
scikit-image >= 0.19.0
h5py >= 3.0.0
vtk >= 9.0.0
PyQt5 >= 5.15.0
PyOpenGL >= 3.1.0
```

### Analysis Requirements (Optional)

```
scikit-learn >= 1.0.0
pandas >= 1.3.0
matplotlib >= 3.5.0
torch >= 1.10.0
fippy (install via: pip install git+https://github.com/gcskoenig/fippy.git)
```

## Citation

If you use this tool in your research, please cite:

```
[Citation information to be added]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/murray-bifurcation-visualization/issues)
- **Documentation**: [docs/](docs/)
- **Examples**: [examples/](examples/)

## Acknowledgments

- **Kimimaro**: Skeleton extraction library
- **PyVista**: 3D visualization framework
- **fippy**: Conditional feature importance analysis

## Project Structure

```
murray-bifurcation-visualization/
├── src/murray_viz/          # Core source code
├── examples/                # Usage examples
├── docs/                    # Documentation
├── data/                    # User data (gitignored)
└── tests/                   # Unit tests
```
