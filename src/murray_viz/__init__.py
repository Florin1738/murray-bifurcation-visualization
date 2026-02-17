"""Murray Bifurcation 3D Visualization

Interactive visualization of Murray's Law compliance at vascular bifurcations.

This package provides tools for analyzing and visualizing Murray's Law at
vascular bifurcations, including:
- 3D vessel surface rendering
- Bifurcation detection and analysis
- Murray metrics computation (ratio, residual, exponent)
- Interactive PyVista-based visualization
- ML-based driver analysis

Example:
    >>> from murray_viz import MurrayBifurcationVisualizer
    >>> viz = MurrayBifurcationVisualizer(
    ...     skeleton_path="data/skeleton.pkl",
    ...     volume_path="data/volume.mat"
    ... )
    >>> viz.detect_bifurcations()
    >>> viz.compute_murray_metrics()
    >>> viz.visualize(metric_name='murray_phi')
"""

from .visualizer import MurrayBifurcationVisualizer
from .config import get_dataset_config, list_datasets, DATASETS, DEFAULT_DATASET

__version__ = "0.1.0"
__author__ = "Johns Hopkins University"
__all__ = ["MurrayBifurcationVisualizer", "get_dataset_config", "list_datasets", "DATASETS", "DEFAULT_DATASET"]
