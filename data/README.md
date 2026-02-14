# Data Directory

Place your skeleton and volume files in this directory.

## File Organization

```
data/
├── your_skeleton.pkl          # Skeleton files (required)
├── your_volume.mat            # Volume files (optional)
└── README.md                  # This file
```

## Skeleton File Format (Required)

Pickle files containing kimimaro skeleton objects with the following structure:

### Required Fields

- **`vertices`**: Nx3 numpy array of 3D coordinates (float)
  - Each row represents a vertex position (x, y, z)
  - Physical coordinates in micrometers or voxels

- **`edges`**: Mx2 numpy array of vertex indices (int)
  - Each row represents an edge connecting two vertices
  - Indices reference rows in the `vertices` array

- **`radius` or `radii`**: N-length numpy array of vessel radii (float)
  - One radius value per vertex
  - Units: typically micrometers
  - Used for Murray Law analysis

### Optional Fields

- **`branching_points`**: Indices of bifurcation/branch points
- **`edge_radii`**: Radii associated with edges
- **`metadata`**: Additional skeleton information

### Example Structure

```python
skeleton = {
    'vertices': np.array([[0, 0, 0], [1, 0, 0], ...]),  # Nx3
    'edges': np.array([[0, 1], [1, 2], ...]),           # Mx2
    'radius': np.array([1.5, 2.0, ...]),                # N-length
}
```

## Volume File Format (Optional)

MATLAB `.mat` files with binary 3D volumes for surface visualization.

### Supported Formats

- **MATLAB v7.2** (scipy.io.loadmat)
- **MATLAB v7.3** (h5py, HDF5-based)

### Required Contents

- **Variable name**: typically `rescaled_vol` (configurable)
- **Data type**: Binary 3D array (uint8 or logical)
- **Content**: 1 = vessel, 0 = background
- **Voxel spacing**: Must provide spacing as (x, y, z) tuple in code

### Example Usage

```python
from murray_viz import MurrayBifurcationVisualizer

viz = MurrayBifurcationVisualizer(
    skeleton_path="data/your_skeleton.pkl",
    volume_path="data/your_volume.mat",
    volume_var='rescaled_vol',       # Variable name in .mat file
    volume_spacing=(16, 16, 8)        # Voxel spacing in microns (x, y, z)
)
```

## Creating Skeleton Files

Skeletons can be generated using the [kimimaro](https://github.com/seung-lab/kimimaro) library:

```python
import kimimaro
import pickle

# Generate skeleton from binary volume
skeletons = kimimaro.skeletonize(
    labels=binary_volume,
    teasar_params={...},
    ...
)

# Save to pickle
with open('data/my_skeleton.pkl', 'wb') as f:
    pickle.dump(skeletons[1], f)  # Save first skeleton
```

See [KIMIMARO_STRUCTURES.md](../docs/KIMIMARO_STRUCTURES.md) for detailed format specifications.

## Example Datasets

Sample datasets can be obtained from:
- [Contact repository maintainers for example data]
- [Link to data repository or DOI - if available]

## Git Ignore

This directory is configured in `.gitignore` to exclude:
- All `.pkl` files (skeleton data)
- All `.mat` files (volume data)
- All `.h5`/`.hdf5` files
- All `.npy`/`.npz` files

Your data files remain **private** and are not committed to version control.

## Data Size Considerations

- Skeleton files: typically 100 KB - 50 MB
- Volume files: can be 100 MB - several GB
- Consider compressing large files or using file compression

## File Naming Conventions

Suggested naming:
- `{species}_{organ}_{specimen}_skeleton.pkl`
- `{species}_{organ}_{specimen}_volume.mat`

Examples:
- `GS40_liver_skeleton.pkl`
- `chicken_liver_vessels_skeleton.pkl`
- `alligator_liver_vessels.pkl`
