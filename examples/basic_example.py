#!/usr/bin/env python3
"""
Murray Bifurcation Visualizer - Example Usage

Demonstrates:
1. Basic visualization with Murray ratio
2. Switching between metrics
3. Changing exponents
4. Configuration options
5. Exporting statistics

Uncomment different sections to try various features.

SETUP INSTRUCTIONS:
1. Place your skeleton (.pkl) files in the ../data/ directory
2. Update file paths in src/murray_viz/config.py (centralized configuration)
3. Each example function uses get_dataset_config() to load paths
4. Available datasets: chicken_liver, GS40_liver, GS55_brain, GS55_skeleton_only
"""

import sys
from pathlib import Path

# Add src directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from murray_viz import MurrayBifurcationVisualizer, get_dataset_config
import numpy as np


def basic_example():
    """Basic usage: visualize Murray ratio with default cube law."""
    print("="*60)
    print("BASIC EXAMPLE: Murray Ratio (φ) with Cube Law (γ=3)")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('chicken_liver')

    # Create visualizer
    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path'],
        volume_var=config['volume_var'],
        volume_spacing=config['volume_spacing']
    )

    # Customize visualization settings for better visibility
    viz.config['colormap'] = 'plasma'

    # Analysis pipeline
    viz.detect_bifurcations()
    viz.compute_murray_metrics(gamma=3.0)
    viz.print_bifurcation_summary()

    # Visualization
    viz.extract_surface()
    viz.assign_metrics_to_surface('murray_phi')
    viz.create_skeleton_mesh()
    viz.visualize(metric_name='murray_phi')


def metric_comparison():
    """Compare different Murray metrics."""
    print("="*60)
    print("METRIC COMPARISON: φ, ε, and γ")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('GS40_liver')

    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path'],
        volume_var=config['volume_var'],
        volume_spacing=config['volume_spacing']
    )

    viz.detect_bifurcations()
    viz.compute_murray_metrics()
    viz.extract_surface()
    viz.create_skeleton_mesh()

    # Visualize each metric
    metrics = ['murray_phi', 'murray_residual', 'murray_gamma']

    for metric in metrics:
        print(f"\n{'='*60}")
        print(f"Visualizing: {metric}")
        print(f"{'='*60}")

        viz.assign_metrics_to_surface(metric)
        viz.visualize(metric_name=metric)


def custom_exponent():
    """Test different Murray exponents."""
    print("="*60)
    print("CUSTOM EXPONENT: Compare γ=2.5 vs γ=3.0")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('GS40_liver')

    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path'],
        volume_var=config['volume_var'],
        volume_spacing=config['volume_spacing']
    )

    viz.detect_bifurcations()
    viz.extract_surface()
    viz.create_skeleton_mesh()

    # Try different exponents
    exponents = [2.5, 3.0, 3.5]

    for gamma in exponents:
        print(f"\n{'='*60}")
        print(f"Murray Ratio with γ={gamma}")
        print(f"{'='*60}")

        viz.compute_murray_metrics(gamma=gamma)

        # Print statistics
        phi = viz.bifurcation_metrics['murray_phi']
        valid = phi[~np.isnan(phi)]
        near_one = np.sum(np.abs(valid - 1.0) < 0.1)

        print(f"Bifurcations near φ=1.0: {near_one}/{len(valid)} ({100*near_one/len(valid):.1f}%)")

        viz.visualize(metric_name='murray_phi', gamma=gamma)


def custom_configuration():
    """Demonstrate configuration options."""
    print("="*60)
    print("CUSTOM CONFIGURATION")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('GS40_liver')

    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path'],
        volume_var=config['volume_var'],
        volume_spacing=config['volume_spacing']
    )

    # Customize appearance
    print("\nCustomizing visualization settings...")

    # Larger bifurcation spheres
    viz.config['bifurcation_sphere_scale'] = 2.0
    print(f"  Sphere scale: {viz.config['bifurcation_sphere_scale']}")

    # More transparent surface
    viz.config['surface_opacity'] = 0.2
    print(f"  Surface opacity: {viz.config['surface_opacity']}")

    # Different colormap
    viz.config['colormap'] = 'plasma'
    print(f"  Colormap: {viz.config['colormap']}")

    # Different skeleton color
    viz.config['skeleton_color'] = 'darkred'
    print(f"  Skeleton color: {viz.config['skeleton_color']}")

    # More neighbors for local radius estimation
    viz.config['local_radius_neighbors'] = 5
    print(f"  Local radius neighbors: {viz.config['local_radius_neighbors']}")

    # Analysis
    viz.detect_bifurcations()
    viz.compute_murray_metrics()
    viz.extract_surface()
    viz.create_skeleton_mesh()

    # Visualize
    viz.assign_metrics_to_surface('murray_phi')
    viz.visualize(metric_name='murray_phi')


def export_statistics():
    """Export bifurcation statistics to CSV."""
    print("="*60)
    print("EXPORT STATISTICS")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('GS40_liver')

    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path'],
        volume_var=config['volume_var'],
        volume_spacing=config['volume_spacing']
    )

    viz.detect_bifurcations()
    viz.compute_murray_metrics()

    # Export to CSV
    import pandas as pd

    data = []
    for i, bif in enumerate(viz.bifurcations):
        row = {
            'bifurcation_id': i,
            'x': bif['position'][0],
            'y': bif['position'][1],
            'z': bif['position'][2],
            'r_parent': bif['r_parent'],
            'r_daughter1': bif['r_daughter1'],
            'r_daughter2': bif['r_daughter2'],
            'local_radius': bif['local_radius_estimate'],
        }

        # Add metrics
        for metric_name, values in viz.bifurcation_metrics.items():
            row[metric_name] = values[i]

        data.append(row)

    df = pd.DataFrame(data)

    # Save
    output_path = 'bifurcation_murray_analysis.csv'
    df.to_csv(output_path, index=False)
    print(f"\nExported {len(df)} bifurcations to {output_path}")

    # Print summary statistics
    print("\nSummary Statistics:")
    print(df.describe())

    # Visualize
    viz.extract_surface()
    viz.create_skeleton_mesh()
    viz.assign_metrics_to_surface('murray_phi')
    viz.visualize(metric_name='murray_phi')


def skeleton_only():
    """Visualize without volume (skeleton + bifurcation analysis only)."""
    print("="*60)
    print("SKELETON-ONLY MODE (No Volume)")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('GS55_skeleton_only')

    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path']
    )

    viz.detect_bifurcations()
    viz.compute_murray_metrics()
    viz.print_bifurcation_summary()

    # Create skeleton mesh
    viz.create_skeleton_mesh()

    # Create simple bifurcation markers
    import pyvista as pv

    bif_centers = np.array([b['position'] for b in viz.bifurcations])
    bif_radii = np.array([b['local_radius_estimate'] for b in viz.bifurcations])
    bif_phi = viz.bifurcation_metrics['murray_phi']

    # Create spheres at bifurcations
    bif_cloud = pv.PolyData(bif_centers)
    bif_cloud['murray_phi'] = bif_phi
    bif_cloud['radius'] = bif_radii

    sphere_geom = pv.Sphere(radius=1.0, phi_resolution=16, theta_resolution=16)
    bif_spheres = bif_cloud.glyph(geom=sphere_geom, scale='radius', factor=2.0)

    # Visualize
    plotter = pv.Plotter(window_size=(1400, 900))
    plotter.set_background('white')

    # Add skeleton
    plotter.add_mesh(viz.skeleton_mesh, color='black', opacity=0.5)

    # Add bifurcation spheres colored by Murray ratio
    valid_phi = bif_phi[~np.isnan(bif_phi)]
    vmin = np.percentile(valid_phi, 5)
    vmax = np.percentile(valid_phi, 95)

    plotter.add_mesh(
        bif_spheres,
        scalars='murray_phi',
        cmap='plasma',
        clim=(vmin, vmax),
        opacity=0.8,
        show_scalar_bar=True,
        scalar_bar_args={'title': 'Murray Ratio φ'}
    )

    plotter.camera_position = 'iso'
    plotter.show_axes()
    plotter.show()


def visualize_murray_exponent():
    """Visualize the best-fit Murray exponent (gamma) at each bifurcation."""
    print("="*60)
    print("MURRAY EXPONENT VISUALIZATION: Best-fit γ per bifurcation")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('chicken_liver')

    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path'],
        volume_var=config['volume_var'],
        volume_spacing=config['volume_spacing']
    )

    # Use plasma colormap
    viz.config['colormap'] = 'plasma'

    # Analysis pipeline
    viz.detect_bifurcations()
    viz.compute_murray_metrics(gamma=3.0)
    viz.print_bifurcation_summary()

    # Report exponent statistics
    gamma_vals = viz.bifurcation_metrics['murray_gamma']
    valid_gamma = gamma_vals[~np.isnan(gamma_vals)]
    if len(valid_gamma) > 0:
        cube_law = np.sum(np.abs(valid_gamma - 3.0) < 0.2)
        print(f"\nMurray Exponent (γ) Statistics:")
        print(f"  Valid estimates: {len(valid_gamma)} / {len(gamma_vals)}")
        print(f"  Range: {valid_gamma.min():.3f} to {valid_gamma.max():.3f}")
        print(f"  Mean: {valid_gamma.mean():.3f}, Median: {np.median(valid_gamma):.3f}")
        print(f"  Near cube law γ=3.0 (±0.2): {cube_law} / {len(valid_gamma)} ({100*cube_law/len(valid_gamma):.1f}%)")

    # Visualization
    viz.extract_surface()
    viz.assign_metrics_to_surface('murray_gamma')
    viz.create_skeleton_mesh()
    viz.visualize(metric_name='murray_gamma')


def quick_analysis():
    """Quick analysis without visualization (for batch processing)."""
    print("="*60)
    print("QUICK ANALYSIS (No Visualization)")
    print("="*60)

    # Load configuration from centralized config
    config = get_dataset_config('GS55_skeleton_only')

    viz = MurrayBifurcationVisualizer(
        skeleton_path=config['skeleton_path'],
        volume_path=config['volume_path']
    )

    viz.detect_bifurcations()
    viz.compute_murray_metrics(gamma=3.0)
    viz.print_bifurcation_summary()

    # Additional analysis
    phi = viz.bifurcation_metrics['murray_phi']
    gamma = viz.bifurcation_metrics['murray_gamma']

    # Count consistency categories
    valid_phi = phi[~np.isnan(phi)]

    consistent = np.sum(np.abs(valid_phi - 1.0) < 0.1)
    under = np.sum(valid_phi < 0.9)
    over = np.sum(valid_phi > 1.1)

    print(f"\nMurray Consistency Categories:")
    print(f"  Consistent (0.9 < φ < 1.1): {consistent} / {len(valid_phi)} ({100*consistent/len(valid_phi):.1f}%)")
    print(f"  Under-perfused (φ < 0.9):   {under} / {len(valid_phi)} ({100*under/len(valid_phi):.1f}%)")
    print(f"  Over-perfused (φ > 1.1):    {over} / {len(valid_phi)} ({100*over/len(valid_phi):.1f}%)")

    # Exponent statistics
    valid_gamma = gamma[~np.isnan(gamma)]
    if len(valid_gamma) > 0:
        cube_law = np.sum(np.abs(valid_gamma - 3.0) < 0.2)
        print(f"\nCube Law Consistency:")
        print(f"  Near γ=3.0 (2.8 < γ < 3.2): {cube_law} / {len(valid_gamma)} ({100*cube_law/len(valid_gamma):.1f}%)")


# ============================================================================
# MAIN: Choose which example to run
# ============================================================================

if __name__ == "__main__":
    print("\nMurray Bifurcation Visualizer - Example Scripts")
    print("="*60)
    print("\nAvailable examples:")
    print("  1. basic_example()              - Basic Murray ratio visualization")
    print("  2. metric_comparison()          - Compare φ, ε, and γ metrics")
    print("  3. custom_exponent()            - Test different exponents")
    print("  4. custom_configuration()       - Customize appearance")
    print("  5. export_statistics()          - Export data to CSV")
    print("  6. skeleton_only()              - Visualize without volume")
    print("  7. visualize_murray_exponent()  - Visualize best-fit Murray exponent γ")
    print("  8. quick_analysis()             - Quick stats without visualization")
    print("="*60)

    # Uncomment the example you want to run:

    basic_example()
    # metric_comparison()
    # custom_exponent()
    # custom_configuration()
    # export_statistics()
    # skeleton_only()
    # visualize_murray_exponent()
    # quick_analysis()

    # Or run the default main from the visualizer
    # from murray_bifurcation_visualizer import main
    # main()
