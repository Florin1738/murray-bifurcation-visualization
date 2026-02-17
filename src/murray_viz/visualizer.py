#!/usr/bin/env python3
"""
Murray's Law Bifurcation Visualizer

Creates an interactive 3D visualization showing vessel surfaces with localized
coloring at bifurcations based on Murray's Law consistency metrics.

Features:
- Semi-translucent vessel surface from binary volume
- Opaque skeleton rendering
- Bifurcation-localized coloring (spheres around branch points)
- Gray coloring everywhere else (using NaN)
- Switchable metrics (Murray ratio φ, Murray exponent γ)
- PyVista-based interactive rendering

Author: Implementation based on rigorous Murray Law analysis
Dependencies: pyvista, numpy, scipy, h5py, scikit-image
"""

import numpy as np
import pyvista as pv
import pickle
from pathlib import Path
from typing import Dict, Tuple, Optional, List, Any
from collections import defaultdict
from scipy.spatial import cKDTree
from scipy.optimize import brentq
from scipy.io import loadmat
import h5py
from skimage import measure


class MurrayBifurcationVisualizer:
    """
    Visualizer for Murray's Law analysis at vascular bifurcations.

    Displays vessel surface colored by Murray metrics only near bifurcations,
    with neutral gray coloring elsewhere.
    """

    def __init__(self, skeleton_path: str, volume_path: Optional[str] = None,
                 volume_var: str = 'rescaled_vol', volume_spacing: Tuple = (16, 16, 16)):
        """
        Initialize the Murray bifurcation visualizer.

        Parameters
        ----------
        skeleton_path : str
            Path to skeleton pickle file
        volume_path : str, optional
            Path to binary volume MATLAB file
        volume_var : str
            Variable name in MATLAB file
        volume_spacing : tuple
            Voxel spacing (sx, sy, sz) in physical units
        """
        self.skeleton_path = skeleton_path
        self.volume_path = volume_path
        self.volume_var = volume_var
        self.volume_spacing = np.array(volume_spacing, dtype=float)

        # Configuration
        self.config = {
            'bifurcation_sphere_scale': 2.5,  # Sphere radius = scale * local_radius
            'local_radius_neighbors': 3,      # Vertices to average for local radius
            'skeleton_tube_radius': 2.0,      # Skeleton rendering radius
            'surface_opacity': 0.25,          # Semi-translucent surface (gray regions)
            'colored_opacity': 0.65,          # Opacity for colored bifurcation regions
            'nan_color': 'lightgray',         # Color for non-bifurcation regions
            'nan_opacity': 0.25,              # Opacity for gray regions
            'skeleton_color': 'black',        # Skeleton color
            'skeleton_opacity': 1.0,          # Opaque skeleton
            'default_metric': 'murray_phi',   # Default: Murray ratio
            'default_gamma': 3.0,             # Default: cube law
            'colormap': 'plasma',              # Color map for metrics
            'background_color': 'white',      # Plotter background
            'window_size': (1400, 900),       # Window dimensions
            'kdtree_k_candidates': 16,        # KD-tree query size for surface transfer
        }

        # Data containers
        self.skeleton = None
        self.skeleton_label = None
        self.volume_data = None
        self.surface_mesh = None
        self.skeleton_mesh = None

        # Bifurcation analysis results
        self.bifurcations = []  # List of bifurcation dicts
        self.bifurcation_metrics = {}  # Metric name -> array of values

        # Load skeleton
        self._load_skeleton()

        # Load volume if provided
        if volume_path:
            self._load_volume()

    def _load_skeleton(self):
        """Load skeleton data from pickle file."""
        print(f"Loading skeleton from {self.skeleton_path}...")

        with open(self.skeleton_path, 'rb') as f:
            skeletons = pickle.load(f)

        if len(skeletons) == 0:
            raise ValueError("No skeletons found in file")

        # Get first skeleton
        self.skeleton_label, self.skeleton = next(iter(skeletons.items()))

        print(f"  Loaded skeleton {self.skeleton_label}")
        print(f"  Vertices: {len(self.skeleton.vertices):,}")
        print(f"  Edges: {len(self.skeleton.edges):,}")

        # Validate radius data
        if hasattr(self.skeleton, 'radius') and self.skeleton.radius is not None:
            self.radius_data = self.skeleton.radius
        elif hasattr(self.skeleton, 'radii') and self.skeleton.radii is not None:
            self.radius_data = self.skeleton.radii
        else:
            raise ValueError("No radius/radii data in skeleton")

        print(f"  Radius range: {self.radius_data.min():.2f} to {self.radius_data.max():.2f}")

    def _load_volume(self):
        """Load binary volume from MATLAB file."""
        print(f"Loading volume from {self.volume_path}...")

        try:
            # Try scipy first (v7.2 and earlier)
            mat = loadmat(self.volume_path)
            self.volume_data = np.asarray(mat[self.volume_var])
            print(f"  Loaded via scipy.io.loadmat")
        except NotImplementedError:
            # MATLAB v7.3 (HDF5)
            with h5py.File(self.volume_path, 'r') as f:
                data = np.array(f[self.volume_var])
                # Transpose for MATLAB column-major order
                self.volume_data = data.T
            print(f"  Loaded via h5py (v7.3) with transpose")

        # Convert to boolean
        if self.volume_data.dtype != bool:
            self.volume_data = (self.volume_data != 0)

        print(f"  Volume shape: {self.volume_data.shape}")
        print(f"  Non-zero voxels: {np.sum(self.volume_data):,}")

    def detect_bifurcations(self):
        """
        Detect degree-3 bifurcations and compute local radii.

        For each bifurcation, identifies the parent (thickest) and two daughters,
        computes local radius estimates, and stores bifurcation geometry.
        """
        print("\nDetecting bifurcations...")

        vertices = self.skeleton.vertices
        edges = self.skeleton.edges
        n_vertices = len(vertices)

        # Build adjacency list
        adjacency = defaultdict(list)
        for v1, v2 in edges:
            adjacency[v1].append(v2)
            adjacency[v2].append(v1)

        # Compute vertex degrees
        degrees = np.array([len(adjacency[i]) for i in range(n_vertices)])

        # Find degree-3 vertices (bifurcations)
        bifurcation_indices = np.where(degrees == 3)[0]
        print(f"  Found {len(bifurcation_indices)} degree-3 bifurcations")

        # Analyze each bifurcation
        self.bifurcations = []

        for bif_idx in bifurcation_indices:
            neighbors = adjacency[bif_idx]

            if len(neighbors) != 3:
                continue

            # Get bifurcation position and radius for distance calculations
            bif_pos = vertices[bif_idx]
            bif_radius = self.radius_data[bif_idx]

            # Identify which branch is the parent (thickest)
            # Use initial neighbor radii as estimate
            initial_radii = [self.radius_data[neighbor_idx] for neighbor_idx in neighbors]
            parent_idx = np.argmax(initial_radii)
            daughter_indices = [i for i in range(3) if i != parent_idx]

            # Second pass: compute local radius for each branch with uniform parameters
            # All branches: 2x bifurcation radius, 5 vertices
            branch_radii = []
            branch_vertices = []  # Track which vertices were used

            for branch_idx, neighbor_idx in enumerate(neighbors):
                # Uniform parameters for all branches
                min_distance = 2.0 * bif_radius  # All branches: march 2x radius
                extra_steps = 5  # Collect 5 vertices for all branches

                # Collect vertices along branch with distances
                branch_verts_all = self._collect_branch_vertices_with_distance(
                    bif_idx, neighbor_idx, adjacency, bif_pos,
                    min_distance=min_distance, extra_steps=extra_steps
                )

                if len(branch_verts_all) >= 1:
                    # Use the vertices that are beyond min_distance
                    verts_for_radius = branch_verts_all
                    local_radius = np.mean([self.radius_data[v] for v in verts_for_radius])
                else:
                    # Fallback: use neighbor vertex directly if branch too short
                    verts_for_radius = [neighbor_idx]
                    local_radius = self.radius_data[neighbor_idx]

                branch_radii.append(local_radius)
                branch_vertices.append(verts_for_radius)  # Store only vertices used for calculation

            branch_radii = np.array(branch_radii)

            r_parent = branch_radii[parent_idx]
            r_daughter1 = branch_radii[daughter_indices[0]]
            r_daughter2 = branch_radii[daughter_indices[1]]

            # Store bifurcation info
            bif_info = {
                'vertex_idx': bif_idx,
                'position': vertices[bif_idx].copy(),
                'r_parent': r_parent,
                'r_daughter1': r_daughter1,
                'r_daughter2': r_daughter2,
                'local_radius_estimate': np.mean(branch_radii),  # For sphere sizing
                'neighbor_indices': neighbors,
                'branch_radii': branch_radii,
                'branch_vertices': branch_vertices,  # List of 3 lists (vertices per branch)
                'parent_branch_idx': parent_idx,
                'daughter_branch_indices': daughter_indices,
            }

            self.bifurcations.append(bif_info)

        print(f"  Analyzed {len(self.bifurcations)} bifurcations")

        # Debug: show some statistics
        if len(self.bifurcations) > 0:
            all_rp = np.array([b['r_parent'] for b in self.bifurcations])
            print(f"  Parent radius range: {all_rp.min():.2f} to {all_rp.max():.2f}")

    def _collect_branch_vertices(self, start_idx: int, next_idx: int,
                                  adjacency: Dict, max_steps: int = 3) -> List[int]:
        """
        Collect vertices along a branch starting from bifurcation.

        Parameters
        ----------
        start_idx : int
            Bifurcation vertex index
        next_idx : int
            First neighbor along branch
        adjacency : dict
            Adjacency list
        max_steps : int
            Maximum steps to follow

        Returns
        -------
        list
            Vertex indices along branch
        """
        visited = {start_idx}
        current = next_idx
        branch_verts = [current]

        for _ in range(max_steps - 1):
            if current is None:
                break

            visited.add(current)
            neighbors = adjacency[current]

            # Find unvisited neighbor (follow the branch)
            next_vert = None
            for n in neighbors:
                if n not in visited:
                    next_vert = n
                    break

            if next_vert is None:
                break

            branch_verts.append(next_vert)
            current = next_vert

        return branch_verts

    def _collect_branch_vertices_with_distance(self, start_idx: int, next_idx: int,
                                                adjacency: Dict, bif_pos: np.ndarray,
                                                min_distance: float, extra_steps: int = 3) -> List[int]:
        """
        Collect vertices along a branch, skipping those within min_distance of bifurcation,
        then collecting extra_steps vertices beyond that threshold.

        Parameters
        ----------
        start_idx : int
            Bifurcation vertex index
        next_idx : int
            First neighbor along branch
        adjacency : dict
            Adjacency list
        bif_pos : np.ndarray
            Position of bifurcation point
        min_distance : float
            Minimum distance from bifurcation before collecting vertices
        extra_steps : int
            Number of vertices to collect after passing min_distance

        Returns
        -------
        list
            Vertex indices beyond min_distance (up to extra_steps vertices)
        """
        visited = {start_idx}
        current = next_idx
        vertices_for_averaging = []
        passed_threshold = False

        # Maximum iterations to prevent infinite loops
        max_iterations = 100
        iteration = 0

        while iteration < max_iterations:
            if current is None:
                break

            # Get distance from bifurcation
            current_pos = self.skeleton.vertices[current]
            distance = np.linalg.norm(current_pos - bif_pos)

            # Check if we've passed the threshold distance
            if distance >= min_distance:
                if not passed_threshold:
                    passed_threshold = True

                # Add this vertex to averaging list
                vertices_for_averaging.append(current)

                # Stop if we've collected enough vertices
                if len(vertices_for_averaging) >= extra_steps:
                    break

            visited.add(current)
            neighbors = adjacency[current]

            # Find unvisited neighbor (follow the branch)
            next_vert = None
            for n in neighbors:
                if n not in visited:
                    next_vert = n
                    break

            if next_vert is None:
                break

            current = next_vert
            iteration += 1

        return vertices_for_averaging

    def compute_murray_metrics(self, gamma: float = 3.0):
        """
        Compute Murray's Law metrics for all bifurcations.

        Computes both the ratio form φ and attempts to estimate the
        best-fit exponent γ for each bifurcation.

        Parameters
        ----------
        gamma : float
            Exponent for Murray ratio computation (default: 3.0 for cube law)
        """
        print(f"\nComputing Murray metrics (gamma={gamma})...")

        n_bif = len(self.bifurcations)

        if n_bif == 0:
            print("  No bifurcations to analyze")
            return

        # Arrays for metrics
        murray_phi = np.zeros(n_bif)
        murray_residual = np.zeros(n_bif)
        murray_gamma = np.full(n_bif, np.nan)

        for i, bif in enumerate(self.bifurcations):
            rp = bif['r_parent']
            r1 = bif['r_daughter1']
            r2 = bif['r_daughter2']

            # Murray ratio φ = (r1^γ + r2^γ) / rp^γ
            phi = (r1**gamma + r2**gamma) / (rp**gamma)
            murray_phi[i] = phi

            # Murray residual ε = 1 - φ
            murray_residual[i] = 1.0 - phi

            # Estimate best-fit exponent γ
            # Solve: rp^γ = r1^γ + r2^γ
            # Only valid if rp > max(r1, r2)
            if rp > max(r1, r2):
                try:
                    # Define function to solve
                    def murray_equation(g):
                        return rp**g - (r1**g + r2**g)

                    # Search for root in reasonable range [1.5, 5.0]
                    gamma_est = brentq(murray_equation, 1.5, 5.0)
                    murray_gamma[i] = gamma_est
                except:
                    # Failed to converge or invalid
                    murray_gamma[i] = np.nan
            else:
                murray_gamma[i] = np.nan

        # Store metrics
        self.bifurcation_metrics = {
            'murray_phi': murray_phi,
            'murray_residual': murray_residual,
            'murray_gamma': murray_gamma,
        }

        # Statistics
        print(f"  Murray ratio φ (γ={gamma}):")
        valid_phi = murray_phi[~np.isnan(murray_phi)]
        if len(valid_phi) > 0:
            print(f"    Range: {valid_phi.min():.3f} to {valid_phi.max():.3f}")
            print(f"    Mean: {valid_phi.mean():.3f}, Median: {np.median(valid_phi):.3f}")
            print(f"    Near 1.0 (consistent): {np.sum(np.abs(valid_phi - 1.0) < 0.1)} / {len(valid_phi)}")

        print(f"  Murray exponent γ:")
        valid_gamma = murray_gamma[~np.isnan(murray_gamma)]
        if len(valid_gamma) > 0:
            print(f"    Valid estimates: {len(valid_gamma)} / {n_bif}")
            print(f"    Range: {valid_gamma.min():.3f} to {valid_gamma.max():.3f}")
            print(f"    Mean: {valid_gamma.mean():.3f}, Median: {np.median(valid_gamma):.3f}")
        else:
            print(f"    No valid estimates (requires rp > max(r1, r2))")

    def extract_surface(self, method: str = 'flying_edges'):
        """
        Extract vessel surface from binary volume using marching cubes/flying edges.

        Parameters
        ----------
        method : str
            Contouring method: 'flying_edges' or 'marching_cubes'
        """
        if self.volume_data is None:
            print("No volume data loaded, cannot extract surface")
            return

        print(f"\nExtracting surface using {method}...")

        # Create PyVista ImageData from volume
        grid = pv.ImageData()
        grid.dimensions = np.array(self.volume_data.shape)
        grid.spacing = self.volume_spacing
        grid.origin = (0.0, 0.0, 0.0)

        # Attach volume as point data (Fortran order for correct indexing)
        grid.point_data['mask'] = self.volume_data.astype(np.float32).flatten(order='F')

        # Extract isosurface at 0.5
        print(f"  Contouring at iso-value 0.5...")
        surface = grid.contour(isosurfaces=[0.5], scalars='mask', method=method)

        # Clean up (remove unused points)
        surface = surface.clean()

        print(f"  Extracted surface: {surface.n_points:,} points, {surface.n_cells:,} cells")

        self.surface_mesh = surface

    def assign_metrics_to_surface(self, metric_name: str = 'murray_phi'):
        """
        Assign bifurcation metric values to surface points using sphere transfer.

        For each surface point:
        - Find nearby bifurcation centers using KD-tree
        - Check if point is inside any bifurcation sphere
        - If yes, assign the metric of the closest containing sphere
        - If no, assign NaN (will render as gray)

        Parameters
        ----------
        metric_name : str
            Which metric to transfer ('murray_phi', 'murray_residual', 'murray_gamma')
        """
        if self.surface_mesh is None:
            print("No surface mesh available")
            return

        if metric_name not in self.bifurcation_metrics:
            print(f"Metric '{metric_name}' not computed")
            return

        print(f"\nAssigning metric '{metric_name}' to surface via sphere transfer...")

        surface_points = self.surface_mesh.points
        n_points = len(surface_points)

        # Metric values for each bifurcation
        metric_values = self.bifurcation_metrics[metric_name]

        # Bifurcation centers and sphere radii
        bif_centers = np.array([b['position'] for b in self.bifurcations])
        bif_radii = np.array([
            self.config['bifurcation_sphere_scale'] * b['local_radius_estimate']
            for b in self.bifurcations
        ])

        print(f"  Surface points: {n_points:,}")
        print(f"  Bifurcations: {len(self.bifurcations)}")
        print(f"  Sphere radius range: {bif_radii.min():.2f} to {bif_radii.max():.2f}")

        # Build KD-tree over bifurcation centers
        kdtree = cKDTree(bif_centers)

        # Initialize all surface scalars to NaN (gray by default)
        surface_scalars = np.full(n_points, np.nan, dtype=np.float32)

        # Query k nearest bifurcations for each surface point
        k = min(self.config['kdtree_k_candidates'], len(self.bifurcations))

        print(f"  Querying {k} nearest bifurcations per surface point...")
        distances, indices = kdtree.query(surface_points, k=k)

        # Handle single neighbor case
        if k == 1:
            distances = distances.reshape(-1, 1)
            indices = indices.reshape(-1, 1)

        # Check containment for each point
        n_assigned = 0

        for i in range(n_points):
            # Check each candidate bifurcation
            for j in range(k):
                bif_idx = indices[i, j]
                dist = distances[i, j]
                radius = bif_radii[bif_idx]

                # If inside this sphere
                if dist <= radius:
                    # Assign metric value (use closest containing sphere)
                    surface_scalars[i] = metric_values[bif_idx]
                    n_assigned += 1
                    break  # Only assign to closest containing sphere

        print(f"  Assigned metric to {n_assigned:,} / {n_points:,} surface points ({100*n_assigned/n_points:.1f}%)")

        # Attach to mesh
        self.surface_mesh[metric_name] = surface_scalars

    def create_skeleton_mesh(self):
        """Create skeleton mesh as tubes for visualization."""
        print("\nCreating skeleton mesh...")

        vertices = self.skeleton.vertices
        edges = self.skeleton.edges

        # Build VTK line format: [2, i, j, 2, i, j, ...]
        lines = []
        for edge in edges:
            v1, v2 = int(edge[0]), int(edge[1])
            lines.extend([2, v1, v2])

        lines_array = np.array(lines, dtype=np.int32)

        # Create PolyData
        poly = pv.PolyData(vertices)
        poly.lines = lines_array

        # Convert to tubes
        skeleton_tube = poly.tube(radius=self.config['skeleton_tube_radius'])

        print(f"  Created skeleton with {len(edges):,} edges as tubes")

        self.skeleton_mesh = skeleton_tube

    def create_bifurcation_points_mesh(self):
        """
        Create mesh showing bifurcation center points (cyan spheres).

        Returns
        -------
        pyvista.PolyData
            Point cloud with spheres at bifurcation vertices
        """
        if len(self.bifurcations) == 0:
            print("No bifurcations detected, cannot create bifurcation points mesh")
            return None

        # Collect bifurcation vertices
        bif_vertex_indices = [bif['vertex_idx'] for bif in self.bifurcations]
        bif_positions = self.skeleton.vertices[bif_vertex_indices]

        print(f"\nCreating bifurcation points mesh...")
        print(f"  Bifurcation vertices: {len(bif_vertex_indices):,}")

        # Create point cloud
        point_cloud = pv.PolyData(bif_positions)

        # Create spheres at each point (5.0x skeleton tube radius for better visibility)
        sphere_radius = self.config['skeleton_tube_radius'] * 5.0
        sphere_geom = pv.Sphere(radius=sphere_radius, phi_resolution=8, theta_resolution=8)
        bif_spheres = point_cloud.glyph(geom=sphere_geom)

        return bif_spheres

    def create_parent_points_mesh(self):
        """
        Create mesh showing points used for parent branch radius calculations (yellow spheres).

        Returns
        -------
        pyvista.PolyData
            Point cloud with spheres at parent branch vertices used for radius averaging
        """
        if len(self.bifurcations) == 0:
            print("No bifurcations detected, cannot create parent points mesh")
            return None

        # Collect parent branch vertices used in Murray calculations
        parent_vertex_indices = set()

        for bif in self.bifurcations:
            # Add parent branch vertices (index is stored in parent_branch_idx)
            parent_branch_idx = bif['parent_branch_idx']
            parent_verts = bif['branch_vertices'][parent_branch_idx]
            parent_vertex_indices.update(parent_verts)

        # Get positions of these vertices
        parent_indices = list(parent_vertex_indices)
        parent_positions = self.skeleton.vertices[parent_indices]

        print(f"\nCreating parent calculation points mesh...")
        print(f"  Parent branch vertices: {len(parent_indices):,}")

        # Create point cloud
        point_cloud = pv.PolyData(parent_positions)

        # Create spheres at each point (5.0x skeleton tube radius for better visibility)
        sphere_radius = self.config['skeleton_tube_radius'] * 5.0
        sphere_geom = pv.Sphere(radius=sphere_radius, phi_resolution=8, theta_resolution=8)
        parent_spheres = point_cloud.glyph(geom=sphere_geom)

        return parent_spheres

    def create_daughter_points_mesh(self):
        """
        Create mesh showing points used for daughter branch radius calculations (magenta spheres).

        Returns
        -------
        pyvista.PolyData
            Point cloud with spheres at daughter branch vertices used for radius averaging
        """
        if len(self.bifurcations) == 0:
            print("No bifurcations detected, cannot create daughter points mesh")
            return None

        # Collect daughter branch vertices used in Murray calculations
        daughter_vertex_indices = set()

        for bif in self.bifurcations:
            # Add daughter branch vertices (indices stored in daughter_branch_indices)
            for daughter_idx in bif['daughter_branch_indices']:
                daughter_verts = bif['branch_vertices'][daughter_idx]
                daughter_vertex_indices.update(daughter_verts)

        # Get positions of these vertices
        daughter_indices = list(daughter_vertex_indices)
        daughter_positions = self.skeleton.vertices[daughter_indices]

        print(f"\nCreating daughter calculation points mesh...")
        print(f"  Daughter branch vertices: {len(daughter_indices):,}")

        # Create point cloud
        point_cloud = pv.PolyData(daughter_positions)

        # Create spheres at each point (5.0x skeleton tube radius for better visibility)
        sphere_radius = self.config['skeleton_tube_radius'] * 5.0
        sphere_geom = pv.Sphere(radius=sphere_radius, phi_resolution=8, theta_resolution=8)
        daughter_spheres = point_cloud.glyph(geom=sphere_geom)

        return daughter_spheres

    def visualize(self, metric_name: str = 'murray_phi', gamma: float = None):
        """
        Create and display the interactive visualization.

        Parameters
        ----------
        metric_name : str
            Metric to visualize ('murray_phi', 'murray_residual', 'murray_gamma')
        gamma : float, optional
            If provided, recompute Murray ratio with this exponent
        """
        # Recompute metrics if gamma changed
        if gamma is not None:
            self.compute_murray_metrics(gamma=gamma)

        # Ensure surface has the metric
        if self.surface_mesh is None or metric_name not in self.surface_mesh.array_names:
            self.assign_metrics_to_surface(metric_name)

        # Ensure skeleton mesh exists
        if self.skeleton_mesh is None:
            self.create_skeleton_mesh()

        print(f"\nCreating visualization for metric '{metric_name}'...")

        # Create plotter
        plotter = pv.Plotter(
            title=f"Murray Bifurcation Analysis - {metric_name}",
            window_size=self.config['window_size']
        )

        plotter.set_background(self.config['background_color'])

        # Determine color limits
        metric_values = self.bifurcation_metrics.get(metric_name, None)

        if metric_values is not None:
            valid_values = metric_values[~np.isnan(metric_values)]

            if len(valid_values) > 0:
                # Robust percentile-based limits
                vmin = np.percentile(valid_values, 5)
                vmax = np.percentile(valid_values, 95)

                # For Murray ratio, center around 1.0
                if metric_name == 'murray_phi':
                    # Symmetric around 1.0
                    max_dev = max(abs(vmin - 1.0), abs(vmax - 1.0))
                    vmin = 1.0 - max_dev
                    vmax = 1.0 + max_dev
                # For residual, center around 0.0
                elif metric_name == 'murray_residual':
                    max_dev = max(abs(vmin), abs(vmax))
                    vmin = -max_dev
                    vmax = max_dev

                clim = (vmin, vmax)
            else:
                clim = (0, 1)
        else:
            clim = (0, 1)

        print(f"  Color limits: {clim[0]:.3f} to {clim[1]:.3f}")

        # Split surface into colored (near bifurcations) and gray (rest) regions
        # for separate opacity control
        metric_values = self.surface_mesh[metric_name]
        colored_mask = ~np.isnan(metric_values)
        gray_mask = np.isnan(metric_values)

        n_colored = np.sum(colored_mask)
        n_gray = np.sum(gray_mask)

        print(f"  Colored points: {n_colored:,} ({100*n_colored/len(metric_values):.1f}%)")
        print(f"  Gray points: {n_gray:,} ({100*n_gray/len(metric_values):.1f}%)")

        # Extract colored region (near bifurcations) - higher opacity
        if n_colored > 0:
            colored_surface = self.surface_mesh.extract_points(colored_mask, adjacent_cells=True)
            plotter.add_mesh(
                colored_surface,
                scalars=metric_name,
                cmap=self.config['colormap'],
                clim=clim,
                opacity=self.config['colored_opacity'],
                show_scalar_bar=True,
                scalar_bar_args={
                    'title': self._get_metric_title(metric_name),
                    'n_labels': 5,
                }
            )

        # Extract gray region (away from bifurcations) - lower opacity
        if n_gray > 0:
            gray_surface = self.surface_mesh.extract_points(gray_mask, adjacent_cells=True)
            plotter.add_mesh(
                gray_surface,
                color=self.config['nan_color'],
                opacity=self.config['nan_opacity'],
                show_scalar_bar=False
            )

        # Add skeleton mesh
        plotter.add_mesh(
            self.skeleton_mesh,
            color=self.config['skeleton_color'],
            opacity=self.config['skeleton_opacity'],
        )

        # Add bifurcation center points (cyan spheres)
        bifurcation_points_mesh = self.create_bifurcation_points_mesh()
        if bifurcation_points_mesh is not None:
            plotter.add_mesh(
                bifurcation_points_mesh,
                color='cyan',
                opacity=0.9,
                show_scalar_bar=False,
            )
            print(f"  Added bifurcation center points in cyan")

        # Add parent branch calculation points (yellow spheres)
        parent_points_mesh = self.create_parent_points_mesh()
        if parent_points_mesh is not None:
            plotter.add_mesh(
                parent_points_mesh,
                color='yellow',
                opacity=0.9,
                show_scalar_bar=False,
            )
            print(f"  Added parent branch calculation points in yellow")

        # Add daughter branch calculation points (magenta spheres)
        daughter_points_mesh = self.create_daughter_points_mesh()
        if daughter_points_mesh is not None:
            plotter.add_mesh(
                daughter_points_mesh,
                color='magenta',
                opacity=0.9,
                show_scalar_bar=False,
            )
            print(f"  Added daughter branch calculation points in magenta")

        # Configure camera
        plotter.camera_position = 'iso'
        # Note: depth peeling disabled due to VTK 9.5.x compatibility issues
        # plotter.enable_depth_peeling()
        plotter.show_axes()

        # Add keyboard shortcuts
        def toggle_skeleton():
            """Toggle skeleton visibility."""
            # This would require tracking actor handle, simplified for now
            pass

        plotter.add_key_event('s', toggle_skeleton)

        # Enable point picking for debugging
        def picking_callback(point):
            """Print metric value at clicked point."""
            print(f"Clicked point: {point}")
            # Find nearest surface point and print its scalar
            if self.surface_mesh is not None and metric_name in self.surface_mesh.array_names:
                idx = self.surface_mesh.find_closest_point(point)
                value = self.surface_mesh[metric_name][idx]
                if np.isnan(value):
                    print(f"  {metric_name}: NaN (gray region, not near bifurcation)")
                else:
                    print(f"  {metric_name}: {value:.3f}")

        plotter.enable_point_picking(callback=picking_callback, show_message=True)

        print("\nInteractive controls:")
        print("  Mouse: Left-drag=rotate, Right-drag=pan, Scroll=zoom")
        print("  Click surface point to see metric value")
        print("  Press 'q' to quit")

        # Show
        plotter.show()

    def _get_metric_title(self, metric_name: str) -> str:
        """Get display title for metric."""
        titles = {
            'murray_phi': 'Murray Ratio φ',
            'murray_residual': 'Murray Residual ε',
            'murray_gamma': 'Murray Exponent γ',
        }
        return titles.get(metric_name, metric_name)

    def print_bifurcation_summary(self):
        """Print summary statistics for bifurcations."""
        print("\n" + "="*60)
        print("BIFURCATION ANALYSIS SUMMARY")
        print("="*60)

        n_bif = len(self.bifurcations)
        print(f"\nTotal bifurcations: {n_bif}")

        if n_bif == 0:
            print("No bifurcations detected")
            return

        # Radius statistics
        all_rp = np.array([b['r_parent'] for b in self.bifurcations])
        all_r1 = np.array([b['r_daughter1'] for b in self.bifurcations])
        all_r2 = np.array([b['r_daughter2'] for b in self.bifurcations])

        print(f"\nRadius statistics:")
        print(f"  Parent:     {all_rp.min():.2f} to {all_rp.max():.2f}, mean {all_rp.mean():.2f}")
        print(f"  Daughter 1: {all_r1.min():.2f} to {all_r1.max():.2f}, mean {all_r1.mean():.2f}")
        print(f"  Daughter 2: {all_r2.min():.2f} to {all_r2.max():.2f}, mean {all_r2.mean():.2f}")

        # Metric statistics
        for metric_name, values in self.bifurcation_metrics.items():
            valid = values[~np.isnan(values)]

            if len(valid) > 0:
                print(f"\n{self._get_metric_title(metric_name)}:")
                print(f"  Valid values: {len(valid)} / {len(values)}")
                print(f"  Range: {valid.min():.3f} to {valid.max():.3f}")
                print(f"  Mean: {valid.mean():.3f}, Median: {np.median(valid):.3f}")
                print(f"  Std: {valid.std():.3f}")

                if metric_name == 'murray_phi':
                    # Count near-consistent bifurcations
                    near_one = np.sum(np.abs(valid - 1.0) < 0.1)
                    print(f"  Near 1.0 (±0.1): {near_one} / {len(valid)} ({100*near_one/len(valid):.1f}%)")

        print("="*60)


def main():
    """Main execution function."""
    from .config import get_dataset_config

    print("Murray's Law Bifurcation Visualizer")
    print("="*60)

    # Load configuration from centralized config
    # To use a different dataset, pass dataset_name to get_dataset_config()
    # e.g., config = get_dataset_config('GS40_liver')
    config = get_dataset_config()  # Uses DEFAULT_DATASET

    # Default metric and gamma
    metric = 'murray_phi'
    gamma = 3.0

    try:
        # Create visualizer
        print(f"\nInitializing visualizer...")
        viz = MurrayBifurcationVisualizer(
            skeleton_path=config['skeleton_path'],
            volume_path=config['volume_path'],
            volume_var=config['volume_var'],
            volume_spacing=config['volume_spacing']
        )

        # Detect bifurcations
        viz.detect_bifurcations()

        # Compute Murray metrics
        viz.compute_murray_metrics(gamma=gamma)

        # Print summary
        viz.print_bifurcation_summary()

        # Extract surface
        viz.extract_surface(method='flying_edges')

        # Assign metrics to surface
        viz.assign_metrics_to_surface(metric_name=metric)

        # Create skeleton mesh
        viz.create_skeleton_mesh()

        # Visualize
        print(f"\nLaunching visualization...")
        viz.visualize(metric_name=metric, gamma=gamma)

        print("\nVisualization complete!")

    except KeyboardInterrupt:
        print("\nInterrupted by user")
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
