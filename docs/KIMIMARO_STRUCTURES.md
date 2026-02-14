# Kimimaro Skeleton Data Structures - Complete Guide

## Table of Contents
1. [Overview](#overview)
2. [Loading Skeleton Data](#loading-skeleton-data)
3. [Skeleton Object Structure](#skeleton-object-structure)
4. [Accessing Skeleton Components](#accessing-skeleton-components)
5. [Radius/Thickness Information](#radiusthickness-information)
6. [Working with Edges](#working-with-edges)
7. [Complete Code Examples](#complete-code-examples)
8. [Visualization Examples](#visualization-examples)
9. [Advanced Operations](#advanced-operations)

---

## Overview

Kimimaro produces skeleton data structures from 3D volumetric images using the TEASAR (Tree-structure Extraction Algorithm for Accurate and Robust skeletons) algorithm. The output is a dictionary of skeleton objects, where each key represents a label ID from the original segmentation.

### Key Data Structure
```python
skeletons = {
    label_id_1: Skeleton_object_1,
    label_id_2: Skeleton_object_2,
    ...
}
```

---

## Loading Skeleton Data

### From Pickle Files
Skeleton data is typically saved as Python pickle files (`.pkl`).

```python
import pickle

# Load skeleton data from pickle file
with open('skeleton_file.pkl', 'rb') as f:
    skeletons = pickle.load(f)

# skeletons is a dictionary: {label_id: skeleton_object}
print(f"Loaded {len(skeletons)} skeletons")
print(f"Label IDs: {list(skeletons.keys())}")
```

### Accessing Individual Skeletons
```python
# Get first skeleton
first_label_id = list(skeletons.keys())[0]
first_skeleton = skeletons[first_label_id]

# Or iterate through all skeletons
for label_id, skeleton in skeletons.items():
    print(f"Processing skeleton {label_id}")
    # Work with skeleton here
```

---

## Skeleton Object Structure

Each `Skeleton` object contains the following key attributes:

### Core Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `vertices` | numpy.ndarray (N×3) | 3D coordinates of skeleton vertices/nodes |
| `edges` | numpy.ndarray (M×2) | Pairs of vertex indices defining connections |
| `radius` or `radii` | numpy.ndarray (N,) | Thickness/radius at each vertex point |
| `vertex_types` | numpy.ndarray (N,) | Type classification of each vertex |

### Example Inspection
```python
skeleton = skeletons[1]  # Get skeleton for label 1

print(f"Number of vertices: {len(skeleton.vertices)}")
print(f"Number of edges: {len(skeleton.edges)}")
print(f"Vertices shape: {skeleton.vertices.shape}")  # (N, 3)
print(f"Edges shape: {skeleton.edges.shape}")        # (M, 2)

# Check which radius attribute exists
if hasattr(skeleton, 'radius') and skeleton.radius is not None:
    print(f"Radius data length: {len(skeleton.radius)}")
elif hasattr(skeleton, 'radii') and skeleton.radii is not None:
    print(f"Radii data length: {len(skeleton.radii)}")
```

---

## Accessing Skeleton Components

### 1. Vertices (3D Points)

Vertices are stored as an N×3 numpy array where each row is [x, y, z] coordinates.

```python
# Access all vertices
vertices = skeleton.vertices
print(f"Vertices array shape: {vertices.shape}")  # (N, 3)

# Access individual vertex
vertex_index = 0
x, y, z = vertices[vertex_index]
print(f"Vertex {vertex_index}: ({x}, {y}, {z})")

# Get all X, Y, Z coordinates separately
all_x = vertices[:, 0]
all_y = vertices[:, 1]
all_z = vertices[:, 2]

# Find bounding box
min_coords = vertices.min(axis=0)  # [min_x, min_y, min_z]
max_coords = vertices.max(axis=0)  # [max_x, max_y, max_z]
```

### 2. Edges (Connections)

Edges are stored as an M×2 numpy array where each row contains indices of connected vertices.

```python
# Access all edges
edges = skeleton.edges
print(f"Edges array shape: {edges.shape}")  # (M, 2)

# Access individual edge
edge_index = 0
vertex_1_idx, vertex_2_idx = edges[edge_index]
print(f"Edge {edge_index} connects vertices {vertex_1_idx} and {vertex_2_idx}")

# Get the 3D coordinates of an edge's endpoints
vertex_1_pos = vertices[vertex_1_idx]  # [x1, y1, z1]
vertex_2_pos = vertices[vertex_2_idx]  # [x2, y2, z2]

# Calculate edge length
import numpy as np
edge_length = np.linalg.norm(vertex_2_pos - vertex_1_pos)
```

---

## Radius/Thickness Information

### Accessing Radius Data

Kimimaro stores radius/thickness information at each vertex. The attribute name may be either `radius` or `radii` depending on the version.

```python
# Safe way to access radius data
def get_radius_data(skeleton):
    """Safely extract radius/thickness data from skeleton."""
    if hasattr(skeleton, 'radius') and skeleton.radius is not None:
        return skeleton.radius
    elif hasattr(skeleton, 'radii') and skeleton.radii is not None:
        return skeleton.radii
    else:
        raise ValueError("No radius data found in skeleton")

# Usage
radius_data = get_radius_data(skeleton)
print(f"Radius data shape: {radius_data.shape}")  # (N,)
```

### Working with Radius

```python
# Get radius at specific vertex
vertex_idx = 100
radius_at_vertex = radius_data[vertex_idx]
print(f"Radius at vertex {vertex_idx}: {radius_at_vertex}")

# Get statistics
print(f"Min radius: {radius_data.min()}")
print(f"Max radius: {radius_data.max()}")
print(f"Mean radius: {radius_data.mean()}")
print(f"Median radius: {np.median(radius_data)}")

# Find vertices with specific radius range
thick_branches_mask = radius_data > 50
thick_branch_indices = np.where(thick_branches_mask)[0]
thick_branch_vertices = vertices[thick_branch_indices]

# Create vertex-to-radius mapping
vertex_radius_dict = {i: radius_data[i] for i in range(len(radius_data))}
```

### Calculating Edge Thickness

```python
# Calculate average radius for each edge
edge_radii = []
for edge in edges:
    v1_idx, v2_idx = edge
    avg_radius = (radius_data[v1_idx] + radius_data[v2_idx]) / 2
    edge_radii.append(avg_radius)

edge_radii = np.array(edge_radii)

# Or as a dictionary
edge_radius_dict = {}
for edge in edges:
    v1_idx, v2_idx = edge
    edge_key = (min(v1_idx, v2_idx), max(v1_idx, v2_idx))
    edge_radius_dict[edge_key] = (radius_data[v1_idx] + radius_data[v2_idx]) / 2
```

---

## Working with Edges

### Analyzing Connectivity

```python
# Calculate vertex degrees (number of connections)
num_vertices = len(vertices)
vertex_degrees = np.zeros(num_vertices, dtype=int)

for edge in edges:
    vertex_degrees[edge[0]] += 1
    vertex_degrees[edge[1]] += 1

# Identify different vertex types
endpoints = np.where(vertex_degrees == 1)[0]      # Terminal points
branch_points = np.where(vertex_degrees > 2)[0]   # Branch/junction points
continuation = np.where(vertex_degrees == 2)[0]   # Path continuation points

print(f"Endpoints: {len(endpoints)}")
print(f"Branch points: {len(branch_points)}")
print(f"Continuation points: {len(continuation)}")
```

### Building Adjacency Information

```python
# Create adjacency list
from collections import defaultdict

adjacency = defaultdict(list)
for v1_idx, v2_idx in edges:
    adjacency[v1_idx].append(v2_idx)
    adjacency[v2_idx].append(v1_idx)

# Get neighbors of a vertex
vertex_idx = 50
neighbors = adjacency[vertex_idx]
print(f"Vertex {vertex_idx} has {len(neighbors)} neighbors: {neighbors}")
```

---

## Complete Code Examples

### Example 1: Load and Analyze Skeleton

```python
import pickle
import numpy as np

def analyze_skeleton(pkl_path):
    """Complete skeleton analysis."""
    # Load data
    with open(pkl_path, 'rb') as f:
        skeletons = pickle.load(f)

    print(f"Loaded {len(skeletons)} skeletons")

    # Analyze each skeleton
    for label_id, skeleton in skeletons.items():
        print(f"\n=== Skeleton {label_id} ===")

        # Basic info
        n_vertices = len(skeleton.vertices)
        n_edges = len(skeleton.edges)
        print(f"Vertices: {n_vertices:,}")
        print(f"Edges: {n_edges:,}")

        # Spatial extent
        min_coords = skeleton.vertices.min(axis=0)
        max_coords = skeleton.vertices.max(axis=0)
        extent = max_coords - min_coords
        print(f"Spatial extent: {extent}")

        # Radius info
        if hasattr(skeleton, 'radius') and skeleton.radius is not None:
            radius = skeleton.radius
        elif hasattr(skeleton, 'radii') and skeleton.radii is not None:
            radius = skeleton.radii
        else:
            print("No radius data available")
            continue

        print(f"Radius range: {radius.min():.2f} to {radius.max():.2f}")
        print(f"Mean radius: {radius.mean():.2f}")

        # Vertex degree analysis
        vertex_degrees = np.zeros(n_vertices)
        for edge in skeleton.edges:
            vertex_degrees[edge[0]] += 1
            vertex_degrees[edge[1]] += 1

        endpoints = np.sum(vertex_degrees == 1)
        branches = np.sum(vertex_degrees > 2)
        print(f"Endpoints: {endpoints}")
        print(f"Branch points: {branches}")

# Usage
analyze_skeleton('GS55_skeleton2.pkl')
```

### Example 2: Extract Specific Data for Plotting

```python
import pickle
import numpy as np
import matplotlib.pyplot as plt

def extract_plotting_data(pkl_path, label_id=None):
    """Extract data for creating plots."""
    with open(pkl_path, 'rb') as f:
        skeletons = pickle.load(f)

    # Get specific skeleton or first one
    if label_id is not None:
        skeleton = skeletons[label_id]
    else:
        skeleton = list(skeletons.values())[0]

    # Extract vertices
    vertices = skeleton.vertices
    x_coords = vertices[:, 0]
    y_coords = vertices[:, 1]
    z_coords = vertices[:, 2]

    # Extract radius
    if hasattr(skeleton, 'radius'):
        radius = skeleton.radius
    else:
        radius = skeleton.radii

    # Create edge lines for plotting
    edge_lines_start = []
    edge_lines_end = []
    edge_colors = []  # Based on average radius

    for edge in skeleton.edges:
        v1_idx, v2_idx = edge
        edge_lines_start.append(vertices[v1_idx])
        edge_lines_end.append(vertices[v2_idx])
        avg_radius = (radius[v1_idx] + radius[v2_idx]) / 2
        edge_colors.append(avg_radius)

    edge_lines_start = np.array(edge_lines_start)
    edge_lines_end = np.array(edge_lines_end)
    edge_colors = np.array(edge_colors)

    return {
        'vertices': vertices,
        'x': x_coords,
        'y': y_coords,
        'z': z_coords,
        'radius': radius,
        'edge_starts': edge_lines_start,
        'edge_ends': edge_lines_end,
        'edge_colors': edge_colors
    }

# Usage
data = extract_plotting_data('skeleton.pkl')

# Example: 3D scatter plot colored by radius
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(data['x'], data['y'], data['z'],
                     c=data['radius'], cmap='viridis', s=10)
plt.colorbar(scatter, label='Radius')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Skeleton colored by radius')
plt.show()
```

### Example 3: Filter by Radius

```python
def filter_skeleton_by_radius(skeleton, min_radius=None, max_radius=None):
    """Extract only parts of skeleton within radius range."""
    # Get radius data
    if hasattr(skeleton, 'radius'):
        radius = skeleton.radius
    else:
        radius = skeleton.radii

    # Create mask
    mask = np.ones(len(radius), dtype=bool)
    if min_radius is not None:
        mask &= (radius >= min_radius)
    if max_radius is not None:
        mask &= (radius <= max_radius)

    # Filter vertices
    filtered_vertex_indices = np.where(mask)[0]
    filtered_vertices = skeleton.vertices[mask]
    filtered_radius = radius[mask]

    # Create index mapping
    old_to_new_idx = {old_idx: new_idx
                      for new_idx, old_idx in enumerate(filtered_vertex_indices)}

    # Filter edges (keep only if both vertices are in filtered set)
    filtered_edges = []
    for edge in skeleton.edges:
        v1, v2 = edge
        if v1 in old_to_new_idx and v2 in old_to_new_idx:
            new_edge = [old_to_new_idx[v1], old_to_new_idx[v2]]
            filtered_edges.append(new_edge)

    filtered_edges = np.array(filtered_edges)

    return {
        'vertices': filtered_vertices,
        'edges': filtered_edges,
        'radius': filtered_radius
    }

# Usage: Get only thick branches (radius > 50)
thick_branches = filter_skeleton_by_radius(skeleton, min_radius=50)
print(f"Thick branches: {len(thick_branches['vertices'])} vertices")
```

---

## Visualization Examples

### Using Vedo (from skeleton_visualizer.py)

```python
from vedo import Points, Lines, Plotter

def visualize_skeleton_vedo(skeleton):
    """Visualize skeleton using vedo library."""
    # Get radius data
    if hasattr(skeleton, 'radius'):
        radius = skeleton.radius
    else:
        radius = skeleton.radii

    # Create edge lines
    edge_lines = Lines(
        [skeleton.vertices[e[0]] for e in skeleton.edges],
        [skeleton.vertices[e[1]] for e in skeleton.edges],
        c='red',
        lw=4
    )
    edge_lines.alpha(0.9)

    # Create vertex points (only show important ones)
    vertex_degrees = np.zeros(len(skeleton.vertices))
    for edge in skeleton.edges:
        vertex_degrees[edge[0]] += 1
        vertex_degrees[edge[1]] += 1

    # Show branch points and endpoints
    important_mask = (vertex_degrees != 2)
    important_vertices = skeleton.vertices[important_mask]

    vertex_points = Points(important_vertices, r=4, c='yellow')
    vertex_points.alpha(0.8)

    # Create plotter and show
    plotter = Plotter(bg='black', axes=1)
    plotter.add(edge_lines)
    plotter.add(vertex_points)
    plotter.show()

# Usage
visualize_skeleton_vedo(skeleton)
```

### Using PyVista (from thickness_based_visualization.py)

```python
import pyvista as pv
import numpy as np

def visualize_with_thickness_pyvista(skeleton):
    """Visualize skeleton with thickness-based coloring."""
    # Get radius
    if hasattr(skeleton, 'radius'):
        radius = skeleton.radius
    else:
        radius = skeleton.radii

    # Create edge mesh
    lines = []
    edge_thickness = []

    for edge in skeleton.edges:
        v1, v2 = edge
        lines.extend([2, v1, v2])
        avg_radius = (radius[v1] + radius[v2]) / 2
        edge_thickness.append(avg_radius)

    lines_array = np.array(lines, dtype=np.int32)
    edge_mesh = pv.PolyData(skeleton.vertices, lines=lines_array)
    edge_mesh.cell_data['thickness'] = np.array(edge_thickness)

    # Create plotter
    plotter = pv.Plotter()
    plotter.add_mesh(
        edge_mesh,
        scalars='thickness',
        cmap='plasma',
        line_width=5,
        scalar_bar_args={'title': 'Branch Thickness'}
    )
    plotter.show()

# Usage
visualize_with_thickness_pyvista(skeleton)
```

---

## Advanced Operations

### Calculate Total Cable Length

```python
def calculate_total_cable_length(skeleton):
    """Calculate total length of all skeleton edges."""
    total_length = 0

    for edge in skeleton.edges:
        v1_idx, v2_idx = edge
        v1_pos = skeleton.vertices[v1_idx]
        v2_pos = skeleton.vertices[v2_idx]
        edge_length = np.linalg.norm(v2_pos - v1_pos)
        total_length += edge_length

    return total_length

# Or use built-in method if available
if hasattr(skeleton, 'cable_length'):
    total_length = np.sum(skeleton.cable_length())
```

### Find Paths Between Points

```python
from collections import deque

def find_path_bfs(skeleton, start_vertex, end_vertex):
    """Find shortest path between two vertices using BFS."""
    # Build adjacency list
    adjacency = defaultdict(list)
    for v1, v2 in skeleton.edges:
        adjacency[v1].append(v2)
        adjacency[v2].append(v1)

    # BFS
    queue = deque([(start_vertex, [start_vertex])])
    visited = {start_vertex}

    while queue:
        current, path = queue.popleft()

        if current == end_vertex:
            return path

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return None  # No path found

# Usage
path = find_path_bfs(skeleton, start_vertex=0, end_vertex=100)
if path:
    path_vertices = skeleton.vertices[path]
    print(f"Found path with {len(path)} vertices")
```

### Calculate Branch Statistics

```python
def calculate_branch_statistics(skeleton):
    """Calculate statistics for each branch segment."""
    # Get radius
    if hasattr(skeleton, 'radius'):
        radius = skeleton.radius
    else:
        radius = skeleton.radii

    branch_stats = []

    for edge in skeleton.edges:
        v1_idx, v2_idx = edge
        v1_pos = skeleton.vertices[v1_idx]
        v2_pos = skeleton.vertices[v2_idx]

        # Length
        length = np.linalg.norm(v2_pos - v1_pos)

        # Radius statistics
        r1, r2 = radius[v1_idx], radius[v2_idx]
        avg_radius = (r1 + r2) / 2
        radius_change = abs(r2 - r1)

        # Volume approximation (as cylinder)
        volume = np.pi * (avg_radius ** 2) * length

        branch_stats.append({
            'edge': (v1_idx, v2_idx),
            'length': length,
            'avg_radius': avg_radius,
            'radius_change': radius_change,
            'volume': volume,
            'start_pos': v1_pos,
            'end_pos': v2_pos
        })

    return branch_stats

# Usage
stats = calculate_branch_statistics(skeleton)
total_volume = sum(s['volume'] for s in stats)
avg_radius = np.mean([s['avg_radius'] for s in stats])
```

---

## Summary of Key Access Patterns

```python
# Quick reference for common operations:

# 1. Load skeletons
with open('skeleton.pkl', 'rb') as f:
    skeletons = pickle.load(f)

# 2. Get a skeleton
skeleton = skeletons[label_id]

# 3. Access vertices (N×3 array)
vertices = skeleton.vertices
x, y, z = vertices[i]  # Get vertex i

# 4. Access edges (M×2 array)
edges = skeleton.edges
v1_idx, v2_idx = edges[i]  # Get edge i

# 5. Access radius (N array)
radius = skeleton.radius if hasattr(skeleton, 'radius') else skeleton.radii
r = radius[i]  # Get radius at vertex i

# 6. Calculate vertex degrees
degrees = np.zeros(len(vertices))
for edge in edges:
    degrees[edge[0]] += 1
    degrees[edge[1]] += 1

# 7. Get edge positions
edge_start = vertices[edges[:, 0]]  # All edge starting points
edge_end = vertices[edges[:, 1]]    # All edge ending points

# 8. Find specific vertex types
endpoints = np.where(degrees == 1)[0]
branches = np.where(degrees > 2)[0]
```

---

## Data Flow Diagram

```
Pickle File (.pkl)
        ↓
    pickle.load()
        ↓
Dictionary {label_id: Skeleton}
        ↓
    skeletons[label_id]
        ↓
    Skeleton Object
        ├── .vertices (N×3)  → [[x1,y1,z1], [x2,y2,z2], ...]
        ├── .edges (M×2)     → [[v1,v2], [v3,v4], ...]
        ├── .radius (N,)     → [r1, r2, r3, ...]
        └── .vertex_types    → [type1, type2, ...]
```

---

## Common Gotchas and Tips

1. **Radius Attribute Name**: Always check for both `radius` and `radii`:
   ```python
   radius = skeleton.radius if hasattr(skeleton, 'radius') else skeleton.radii
   ```

2. **Edge Indices**: Remember edges contain vertex indices, not coordinates:
   ```python
   edge = edges[0]  # Returns [vertex_idx_1, vertex_idx_2]
   pos1 = vertices[edge[0]]  # Get actual 3D position
   ```

3. **Vertex Indexing**: Vertex indices in edges correspond to rows in the vertices array.

4. **Filtering**: When filtering vertices, remember to update edge indices accordingly.

5. **Memory**: Large skeletons can be memory-intensive. Consider downsampling or filtering.

---

This guide provides everything needed to work with kimimaro skeleton data, from basic loading to advanced analysis and visualization!
