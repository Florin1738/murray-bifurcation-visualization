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

### 1. Install Dependencies

```bash
# Clone the repository
git clone https://github.com/Florin1738/murray-bifurcation-visualization.git
cd murray-bifurcation-visualization

# Install core dependencies
pip install -r requirements.txt

# Optional: Install ML analysis dependencies
pip install -r requirements-analysis.txt
```

### 2. Add Your Data

Place your skeleton file(s) in the `data/` directory:
```bash
data/
└── your_skeleton.pkl  # Required: kimimaro skeleton file
└── your_volume.mat    # Optional: binary volume for surface rendering
```

### 3. Run the Example

**That's it!** Just run the example script:

```bash
cd examples
python basic_example.py
```

The script includes several pre-configured examples. At the bottom of `basic_example.py`, **uncomment the function you want to run**:

```python
if __name__ == "__main__":
    basic_example()           # ← Default: Basic Murray visualization
    # metric_comparison()     # Compare φ, ε, and γ metrics
    # custom_exponent()       # Test different exponents
    # custom_configuration()  # Customize appearance
    # export_statistics()     # Export data to CSV
    # skeleton_only()         # Visualize without volume
    # quick_analysis()        # Fast batch processing
```

**Before running**, update the centralized configuration file `src/murray_viz/config.py`:
```python
# In config.py, add your dataset:
DATASETS = {
    'your_data': {
        'skeleton_path': 'path/to/your_skeleton.pkl',
        'volume_path': 'path/to/your_volume.mat',  # or None
        'volume_var': 'rescaled_vol',
        'volume_spacing': (16, 16, 8)  # IMPORTANT: Set correct voxel spacing!
    }
}
```

⚠️ **Important**: Set the correct `volume_spacing` (x, y, z) in microns for your data! If the skeleton and volume don't overlap when visualized, this is usually a sign that the voxel spacing is incorrect.

That's it! The visualization will open in a 3D interactive window.

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
- **Spacing**: Must specify correct voxel dimensions (x, y, z) in microns

⚠️ **Voxel Spacing is Critical**: If skeleton and volume don't overlap in the visualization, the voxel spacing is likely incorrect. Double-check your imaging parameters!

See [data/README.md](data/README.md) for detailed format specifications.

## What the Examples Include

The `basic_example.py` script includes 7 different analyses:

1. **`basic_example()`** - Simple Murray ratio visualization
2. **`metric_comparison()`** - Compare φ, ε, and γ side-by-side
3. **`custom_exponent()`** - Test different Murray exponents
4. **`custom_configuration()`** - Customize colors and appearance
5. **`export_statistics()`** - Save data to CSV for further analysis
6. **`skeleton_only()`** - Visualize without volume data
7. **`quick_analysis()`** - Fast batch processing without visualization

Just uncomment the one you want at the bottom of the file and run it!

---

## For Advanced Users

### Custom Python Scripts

You can write your own analysis scripts using the package:

```python
from murray_viz import MurrayBifurcationVisualizer, get_dataset_config

# Option 1: Use centralized configuration
config = get_dataset_config('your_data')
viz = MurrayBifurcationVisualizer(
    skeleton_path=config['skeleton_path'],
    volume_path=config['volume_path'],
    volume_var=config['volume_var'],
    volume_spacing=config['volume_spacing']
)

# Option 2: Direct paths (for one-off scripts)
viz = MurrayBifurcationVisualizer(
    skeleton_path="data/your_skeleton.pkl",
    volume_path="data/your_volume.mat",  # Optional
    volume_var='rescaled_vol',
    volume_spacing=(16, 16, 8)
)

viz.detect_bifurcations()
viz.compute_murray_metrics(gamma=3.0)
viz.visualize(metric_name='murray_phi')
```

### Customization Options

```python
# Customize visualization appearance
viz.config['bifurcation_sphere_scale'] = 2.5  # Sphere size
viz.config['surface_opacity'] = 0.25           # Transparency
viz.config['colormap'] = 'coolwarm'            # Color scheme
viz.config['window_size'] = (1400, 900)        # Window size
```

### Detailed Documentation

For in-depth technical details, see the [docs/](docs/) directory:
- **[Visualization Guide](docs/VISUALIZATION_GUIDE.md)** - Complete API documentation
- **[Implementation Details](docs/IMPLEMENTATION.md)** - Technical architecture
- **[Exponent Fitting](docs/EXPONENT_FITTING.md)** - Murray exponent calculation
- **[Data Structures](docs/KIMIMARO_STRUCTURES.md)** - Skeleton file format

## Requirements

**Python 3.9+** and the packages in `requirements.txt` (installed automatically with `pip install -r requirements.txt`)

Main packages:
- **PyVista** - 3D visualization
- **NumPy/SciPy** - Scientific computing
- **scikit-image** - Image processing

Optional (for ML analysis):
- **scikit-learn** - Machine learning
- **pandas** - Data analysis
- **PyTorch** - Deep learning

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

## Project Structure

```
murray-bifurcation-visualization/
├── src/murray_viz/          # Core source code
├── examples/                # Usage examples
├── docs/                    # Documentation
├── data/                    # User data (gitignored)
└── tests/                   # Unit tests
```
