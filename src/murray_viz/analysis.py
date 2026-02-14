"""
Murray Exponent Driver Analysis - Diagnostic Plots
===================================================

This script analyzes what factors drive variation in the fitted Murray exponent (alpha)
at vascular bifurcations across multiple species samples.

Key Analyses:
1. Alpha vs Parent Radius (2D hexbin) - Tests for scale-dependent effects
2. Alpha vs Daughter Asymmetry (2D hexbin) - Tests for geometry-dependent effects
3. Permutation Importance (bar chart) - Quantifies which features predict alpha
4. Conditional Feature Importance (CFI vs PFI) - Separates unique from redundant
   predictive signal using fippy's conditional permutation (addresses correlated features)

Mathematical Framework:
- Generalized Murray relation: r_parent^α ≈ r_child1^α + r_child2^α
- Classical Murray's law assumes α = 3
- Local alpha is fitted per bifurcation by solving: r_p^α = r_1^α + r_2^α
- Uses log-space root finding for numerical stability

Author: Generated for vascular morphology comparative analysis
Date: 2026-01-22
"""

import os
import sys

# CRITICAL: Set matplotlib backend BEFORE importing matplotlib
os.environ['MPLBACKEND'] = 'Agg'

import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg', force=True)

# Configure matplotlib
import matplotlib as mpl
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 10
mpl.rcParams['text.usetex'] = False

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_agg import FigureCanvasAgg
from collections import defaultdict
from typing import Dict, List, Tuple, Any, Optional
from scipy.optimize import brentq
from scipy.spatial.distance import cosine
import warnings

warnings.filterwarnings('ignore')

print("Matplotlib backend:", mpl.get_backend(), flush=True)
print("Matplotlib version:", mpl.__version__, flush=True)

# Try to import scikit-learn for permutation importance
try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import train_test_split, RepeatedKFold
    from sklearn.inspection import permutation_importance
    from sklearn.pipeline import make_pipeline
    from sklearn.metrics import mean_squared_error
    SKLEARN_AVAILABLE = True
    print("scikit-learn: Available", flush=True)
except ImportError:
    SKLEARN_AVAILABLE = False
    print("scikit-learn: Not available - permutation importance plot will be skipped", flush=True)

# Try to import fippy for conditional feature importance (CFI)
# CFI uses conditional sampling so permuted feature values remain plausible given
# other features, avoiding the correlation-induced misattribution of naive PFI.
try:
    from fippy.explainers import Explainer
    from fippy.samplers import UnivRFSampler, ContUnivRFSampler, SequentialSampler
    import pandas as pd
    FIPPY_AVAILABLE = True
    print("fippy: Available", flush=True)
except ImportError:
    FIPPY_AVAILABLE = False
    print("fippy: Not available - conditional feature importance plots will be skipped", flush=True)


# ============================================================================
# CONFIGURATION
# ============================================================================
#
# NOTE: Update these paths to point to your skeleton files in the data/ directory.
# You can use relative paths from the repository root:
#   Example: '../../data/alligator_liver_vessels.pkl'
# Or absolute paths to your data location.

SKELETON_FILES = {
    'Alligator': '../../data/alligator_liver_vessels.pkl',
    'GS40 Monkey': '../../data/GS40_liver_skeleton.pkl',
    'GS55 Monkey': '../../data/GS55_liver_skeleton.pkl',
    'Chicken': '../../data/chicken_liver_skeleton.pkl',
    'Turtle': '../../data/turtle_liver_vessels.pkl'
}

VOXEL_SPACING = {
    'Alligator': (2.0, 2.0, 1.0),
    'GS40 Monkey': (2.0, 2.0, 1.0),
    'GS55 Monkey': (2.0, 2.0, 1.0),
    'Chicken': (2.0, 2.0, 1.0),
    'Turtle': (2.0, 2.0, 1.0)
}

# Analysis parameters
JUNCTION_RADIUS_POINTS = 5  # Points to average for local radius estimation
MIN_PARENT_RADIUS = 2.0     # Minimum parent radius to include (filter noise)
ALPHA_BRACKET = (0.5, 6.0)  # Bracket for alpha root finding
BRANCH_DIRECTION_POINTS = 3 # Points to use for estimating branch direction
TORTUOSITY_POINTS = 10      # Points to use for tortuosity calculation

# Output settings
OUTPUT_DIR = r'c:\Users\Florin\OneDrive - Johns Hopkins\Documents\kimimaro_skeletonization\plots'
INDIVIDUAL_PLOTS_DIR = os.path.join(OUTPUT_DIR, 'individual_samples')
FIGURE_DPI = 300
FIGURE_FORMAT = 'png'

# ============================================================================
# FEATURE IMPORTANCE SETTINGS
# ============================================================================
#
# Feature importance is computed using cross-validation for both:
#   - PFI (Permutation Feature Importance): Marginal importance via naive shuffling
#   - CFI (Conditional Feature Importance): Conditional importance accounting for correlations
#
# Cross-validation parameters (K-fold, n_repeats) are AUTO-DETERMINED based on sample size:
#
#   Sample Size | PFI: K×Repeats | CFI: K×Repeats | Notes
#   ------------|----------------|----------------|------------------
#   n = 50-79   | 3×5            | Skip           | Minimum viable
#   n = 80-99   | 3×7            | Skip           | Conservative
#   n = 100-119 | 5×10           | 3×7            | Typical range
#   n = 120-149 | 5×10           | 4×7            | Standard
#   n ≥ 150     | 5×10           | 5×10           | Both use same
#
# CFI uses smaller K (larger training sets) because it trains additional Random
# Forests to learn conditional distributions p(x_j | x_{-j}) for each feature.
#

# Default metric for importance plots:
#   'pct_increase' - % increase in MSE (best for comparing across samples)
#   'delta_rmse'   - RMSE increase in alpha units (most interpretable)
#   'delta_mse'    - Raw MSE increase
IMPORTANCE_METRIC = 'pct_increase'

# Minimum sample sizes:
MIN_SAMPLES_PFI = 50   # Minimum for PFI analysis
MIN_SAMPLES_CFI = 100  # Minimum for CFI analysis (needs larger samples for conditional samplers)

# Permutation repeats per fold:
IMPORTANCE_PERM_REPEATS = 20


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def safe_save_figure(fig, filepath: str, dpi: int = 300) -> bool:
    """
    Save figure using direct buffer access to avoid matplotlib hanging issues.
    """
    try:
        from PIL import Image

        fig.set_dpi(dpi)
        canvas = FigureCanvasAgg(fig)
        canvas.draw()

        buf = canvas.buffer_rgba()
        width, height = canvas.get_width_height()
        image = Image.frombuffer('RGBA', (width, height), buf, 'raw', 'RGBA', 0, 1)

        if filepath.lower().endswith('.png'):
            image.save(filepath, 'PNG')
        elif filepath.lower().endswith('.jpg') or filepath.lower().endswith('.jpeg'):
            rgb_image = image.convert('RGB')
            rgb_image.save(filepath, 'JPEG', quality=95)
        else:
            image.save(filepath, 'PNG')

        return True
    except Exception as e:
        print(f"  ERROR in safe_save_figure: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return False


def load_skeleton(filepath: str) -> Any:
    """Load skeleton from pickle file."""
    with open(filepath, 'rb') as f:
        skeletons = pickle.load(f)

    first_label = list(skeletons.keys())[0]
    skeleton = skeletons[first_label]

    print(f"  Loaded skeleton: {len(skeleton.vertices)} vertices, {len(skeleton.edges)} edges")
    return skeleton


def get_radius_data(skeleton) -> np.ndarray:
    """Safely extract radius data from skeleton object."""
    if hasattr(skeleton, 'radius') and skeleton.radius is not None:
        return skeleton.radius
    elif hasattr(skeleton, 'radii') and skeleton.radii is not None:
        return skeleton.radii
    else:
        raise ValueError("No radius data found in skeleton")


def compute_vertex_degrees(edges: np.ndarray, n_vertices: int) -> np.ndarray:
    """Compute degree (number of connections) for each vertex."""
    degrees = np.zeros(n_vertices, dtype=int)
    for edge in edges:
        degrees[edge[0]] += 1
        degrees[edge[1]] += 1
    return degrees


def build_adjacency_list(edges: np.ndarray, n_vertices: int) -> Dict[int, List[int]]:
    """Build adjacency list representation."""
    adjacency = defaultdict(list)
    for v1, v2 in edges:
        adjacency[v1].append(v2)
        adjacency[v2].append(v1)
    return adjacency


# ============================================================================
# BRANCH GEOMETRY CALCULATIONS
# ============================================================================

def estimate_local_radius_at_junction(junction_idx: int, neighbor_idx: int,
                                      adjacency: Dict, radius_data: np.ndarray,
                                      vertices: np.ndarray,
                                      k: int = 5) -> float:
    """
    Estimate local radius for a branch near junction.

    Marches away from bifurcation by 2x junction radius, then averages k points.
    This reduces bias from junction artifacts.
    """
    junction_pos = vertices[junction_idx]
    junction_radius = radius_data[junction_idx]
    min_distance = 2.0 * junction_radius

    visited = {junction_idx}
    current = neighbor_idx
    vertices_for_averaging = []
    passed_threshold = False

    max_iterations = 100
    iteration = 0

    while iteration < max_iterations:
        if current is None:
            break

        current_pos = vertices[current]
        distance = np.linalg.norm(current_pos - junction_pos)

        if distance >= min_distance:
            if not passed_threshold:
                passed_threshold = True

            vertices_for_averaging.append(current)

            if len(vertices_for_averaging) >= k:
                break

        visited.add(current)
        neighbors = adjacency[current]

        next_vert = None
        for n in neighbors:
            if n not in visited:
                next_vert = n
                break

        if next_vert is None:
            break

        current = next_vert
        iteration += 1

    if len(vertices_for_averaging) > 0:
        local_radii = [radius_data[i] for i in vertices_for_averaging]
        return np.mean(local_radii)
    else:
        # Fallback: use neighbor vertex directly
        return radius_data[neighbor_idx]


def estimate_branch_direction(junction_idx: int, neighbor_idx: int,
                              adjacency: Dict, vertices: np.ndarray,
                              k: int = 3) -> np.ndarray:
    """
    Estimate direction vector for a branch emanating from junction.

    Uses first k vertices along the branch to compute mean direction.
    Returns normalized direction vector.
    """
    visited = {junction_idx}
    current = neighbor_idx
    branch_vertices = [vertices[junction_idx], vertices[neighbor_idx]]

    max_iterations = k
    iteration = 0

    while iteration < max_iterations:
        if current is None:
            break

        visited.add(current)
        neighbors = adjacency[current]

        next_vert = None
        for n in neighbors:
            if n not in visited:
                next_vert = n
                break

        if next_vert is None:
            break

        branch_vertices.append(vertices[next_vert])
        current = next_vert
        iteration += 1

    if len(branch_vertices) < 2:
        # Fallback: use just junction to neighbor
        direction = vertices[neighbor_idx] - vertices[junction_idx]
    else:
        # Use overall direction from junction to last point
        direction = branch_vertices[-1] - branch_vertices[0]

    # Normalize
    norm = np.linalg.norm(direction)
    if norm > 0:
        return direction / norm
    else:
        return np.array([0., 0., 0.])


def compute_branch_angle(direction1: np.ndarray, direction2: np.ndarray) -> float:
    """
    Compute angle between two branch direction vectors in degrees.

    Uses dot product: cos(theta) = (v1 · v2) / (|v1| |v2|)
    """
    norm1 = np.linalg.norm(direction1)
    norm2 = np.linalg.norm(direction2)

    if norm1 == 0 or norm2 == 0:
        return np.nan

    cos_angle = np.dot(direction1, direction2) / (norm1 * norm2)
    # Clamp to [-1, 1] to avoid numerical issues with arccos
    cos_angle = np.clip(cos_angle, -1.0, 1.0)

    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)

    return angle_deg


def estimate_local_tortuosity(junction_idx: int, neighbor_idx: int,
                              adjacency: Dict, vertices: np.ndarray,
                              k: int = 10) -> float:
    """
    Estimate local tortuosity near junction.

    Tortuosity = path_length / straight_line_distance

    Computes path length from summed point-to-point distances over k points.
    Computes straight-line distance from first to last point.
    """
    visited = {junction_idx}
    current = neighbor_idx
    path_points = [vertices[junction_idx], vertices[neighbor_idx]]

    max_iterations = k
    iteration = 0

    while iteration < max_iterations:
        if current is None:
            break

        visited.add(current)
        neighbors = adjacency[current]

        next_vert = None
        for n in neighbors:
            if n not in visited:
                next_vert = n
                break

        if next_vert is None:
            break

        path_points.append(vertices[next_vert])
        current = next_vert
        iteration += 1

    if len(path_points) < 2:
        return 1.0  # No tortuosity

    # Compute path length
    path_length = 0.0
    for i in range(len(path_points) - 1):
        segment_length = np.linalg.norm(path_points[i+1] - path_points[i])
        path_length += segment_length

    # Compute straight-line distance
    straight_line_distance = np.linalg.norm(path_points[-1] - path_points[0])

    if straight_line_distance == 0:
        return 1.0

    tortuosity = path_length / straight_line_distance
    return tortuosity


def estimate_radius_taper(junction_idx: int, neighbor_idx: int,
                          adjacency: Dict, radius_data: np.ndarray,
                          vertices: np.ndarray,
                          k: int = 5) -> float:
    """
    Estimate radius taper along a branch.

    Taper = (r_near - r_far) / distance

    Positive taper means radius decreases along branch (typical).
    Negative taper means radius increases (atypical, may indicate venous).
    """
    # Get near radius (close to junction)
    r_near = estimate_local_radius_at_junction(junction_idx, neighbor_idx,
                                                adjacency, radius_data, vertices, k=2)

    # Get far radius (away from junction)
    visited = {junction_idx}
    current = neighbor_idx
    path_distance = 0.0

    max_iterations = k * 2
    iteration = 0

    while iteration < max_iterations:
        if current is None:
            break

        visited.add(current)
        neighbors = adjacency[current]

        next_vert = None
        for n in neighbors:
            if n not in visited:
                next_vert = n
                break

        if next_vert is None:
            break

        # Accumulate distance
        path_distance += np.linalg.norm(vertices[next_vert] - vertices[current])
        current = next_vert
        iteration += 1

    if current is None or path_distance == 0:
        return 0.0  # No taper measurable

    r_far = radius_data[current]

    taper = (r_near - r_far) / path_distance
    return taper


# ============================================================================
# MURRAY EXPONENT FITTING
# ============================================================================

def fit_local_alpha(r_parent: float, r_child1: float, r_child2: float,
                   bracket: Tuple[float, float] = (0.5, 6.0)) -> float:
    """
    Fit local Murray exponent alpha for a single bifurcation.

    Solves: r_p^α = r_1^α + r_2^α

    Uses log-space formulation for numerical stability:
    f(α) = log(r_1^α + r_2^α) - α·log(r_p) = 0

    Parameters
    ----------
    r_parent : float
        Parent branch radius
    r_child1, r_child2 : float
        Child branch radii
    bracket : tuple
        Search bracket for alpha

    Returns
    -------
    alpha : float
        Fitted exponent, or np.nan if no solution found
    """
    if r_parent <= 0 or r_child1 <= 0 or r_child2 <= 0:
        return np.nan

    def objective(alpha):
        """f(α) = log(r_1^α + r_2^α) - α·log(r_p)"""
        try:
            term1 = r_child1**alpha + r_child2**alpha
            if term1 <= 0:
                return np.nan
            return np.log(term1) - alpha * np.log(r_parent)
        except (OverflowError, ValueError):
            return np.nan

    # Check if objective is valid at bracket endpoints
    f_low = objective(bracket[0])
    f_high = objective(bracket[1])

    if np.isnan(f_low) or np.isnan(f_high):
        return np.nan

    # Check if signs are opposite (required for Brent's method)
    if f_low * f_high > 0:
        # No sign change in bracket - no unique solution
        return np.nan

    try:
        alpha_opt = brentq(objective, bracket[0], bracket[1], xtol=1e-6, maxiter=100)
        return alpha_opt
    except (ValueError, RuntimeError):
        return np.nan


# ============================================================================
# BIFURCATION FEATURE EXTRACTION
# ============================================================================

def extract_bifurcation_features(skeleton, voxel_spacing: Tuple[float, float, float],
                                sample_name: str) -> List[Dict]:
    """
    Extract comprehensive features for each degree-3 bifurcation.

    For each bifurcation, computes:
    - Parent and child radii (using offset from junction)
    - Fitted local alpha
    - phi_3 (Murray ratio at alpha=3)
    - Asymmetry index
    - Branch angle between daughters
    - Local tortuosity
    - Radius taper

    Returns
    -------
    bifurcation_features : list of dict
        Each dict contains all computed features for one bifurcation
    """
    vertices = skeleton.vertices
    edges = skeleton.edges
    radius_data = get_radius_data(skeleton)
    n_vertices = len(vertices)

    # Compute degrees
    degrees = compute_vertex_degrees(edges, n_vertices)

    # Build adjacency
    adjacency = build_adjacency_list(edges, n_vertices)

    # Find degree-3 junctions
    bifurcation_nodes = np.where(degrees == 3)[0]

    bifurcation_features = []

    print(f"  Processing {len(bifurcation_nodes)} bifurcations...")

    for junction_idx in bifurcation_nodes:
        neighbors = adjacency[junction_idx]

        if len(neighbors) != 3:
            continue

        # Estimate local radius for each branch
        local_radii = []
        for neighbor in neighbors:
            local_r = estimate_local_radius_at_junction(
                junction_idx, neighbor, adjacency, radius_data, vertices,
                k=JUNCTION_RADIUS_POINTS
            )
            local_radii.append(local_r)

        # Assign parent as largest radius (assumption: arterial tapering)
        sorted_indices = np.argsort(local_radii)[::-1]
        r_parent = local_radii[sorted_indices[0]]
        r_child1 = local_radii[sorted_indices[1]]
        r_child2 = local_radii[sorted_indices[2]]

        # Filter by minimum parent radius
        if r_parent < MIN_PARENT_RADIUS:
            continue

        # Compute fitted local alpha
        alpha = fit_local_alpha(r_parent, r_child1, r_child2, bracket=ALPHA_BRACKET)

        # Compute phi_3 (classical Murray ratio)
        phi_3 = (r_child1**3 + r_child2**3) / r_parent**3

        # Compute asymmetry index: min(r1, r2) / max(r1, r2)
        asymmetry = min(r_child1, r_child2) / max(r_child1, r_child2)

        # Estimate branch directions for daughters
        daughter_neighbor_1 = neighbors[sorted_indices[1]]
        daughter_neighbor_2 = neighbors[sorted_indices[2]]

        direction1 = estimate_branch_direction(
            junction_idx, daughter_neighbor_1, adjacency, vertices,
            k=BRANCH_DIRECTION_POINTS
        )
        direction2 = estimate_branch_direction(
            junction_idx, daughter_neighbor_2, adjacency, vertices,
            k=BRANCH_DIRECTION_POINTS
        )

        # Compute angle between daughters
        branch_angle = compute_branch_angle(direction1, direction2)

        # Estimate local tortuosity for parent branch
        parent_neighbor = neighbors[sorted_indices[0]]
        tortuosity = estimate_local_tortuosity(
            junction_idx, parent_neighbor, adjacency, vertices,
            k=TORTUOSITY_POINTS
        )

        # Estimate radius taper for parent branch
        radius_taper = estimate_radius_taper(
            junction_idx, parent_neighbor, adjacency, radius_data, vertices,
            k=5
        )

        bifurcation_features.append({
            'sample': sample_name,
            'junction_idx': junction_idx,
            'r_parent': r_parent,
            'r_child1': r_child1,
            'r_child2': r_child2,
            'alpha': alpha,
            'phi_3': phi_3,
            'asymmetry': asymmetry,
            'branch_angle': branch_angle,
            'tortuosity': tortuosity,
            'radius_taper': radius_taper,
            'log_r_parent': r_parent  # Store raw parent radius (not log-transformed)
        })

    print(f"  Extracted features for {len(bifurcation_features)} bifurcations")

    return bifurcation_features


# ============================================================================
# CORE DATA EXTRACTION AND COMPUTATION HELPERS
# ============================================================================

def extract_plot_data(features: List[Dict], x_key: str, y_key: str = 'alpha',
                      require_positive: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract and filter data arrays from feature dictionaries.

    Returns filtered x and y arrays where both are finite (and optionally positive).
    """
    x = np.array([f[x_key] for f in features])
    y = np.array([f[y_key] for f in features])

    mask = np.isfinite(x) & np.isfinite(y)
    if require_positive:
        mask &= (x > 0) & (y > 0)

    return x[mask], y[mask]


def compute_binned_medians(x: np.ndarray, y: np.ndarray, n_bins: int = 12,
                           log_scale: bool = False, min_count: int = 5
                           ) -> Tuple[List[float], List[float]]:
    """
    Compute binned median trend line.

    Parameters
    ----------
    x, y : arrays
        Data arrays
    n_bins : int
        Number of bins
    log_scale : bool
        If True, use log-spaced bins (for radius); else linear (for asymmetry)
    min_count : int
        Minimum points required per bin

    Returns
    -------
    centers, medians : lists
        Bin centers and median y values for bins with sufficient data
    """
    if log_scale:
        bins = np.logspace(np.log10(x.min()), np.log10(x.max()), n_bins + 1)
        get_center = lambda lo, hi: np.sqrt(lo * hi)  # geometric mean
    else:
        bins = np.linspace(x.min(), x.max(), n_bins + 1)
        get_center = lambda lo, hi: (lo + hi) / 2

    centers, medians = [], []
    for i in range(n_bins):
        mask = (x >= bins[i]) & (x < bins[i + 1])
        if np.sum(mask) >= min_count:
            centers.append(get_center(bins[i], bins[i + 1]))
            medians.append(np.median(y[mask]))

    return centers, medians


def prepare_regression_data(features: List[Dict],
                            feature_names: List[str] = None
                            ) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Prepare feature matrix X and target vector y for regression.

    Returns X, y, and the feature names used.
    """
    if feature_names is None:
        feature_names = ['log_r_parent', 'asymmetry', 'branch_angle', 'tortuosity', 'radius_taper']

    X_list, y_list = [], []

    for f in features:
        if not np.isfinite(f.get('alpha', np.nan)):
            continue

        vec = [f.get(name, np.nan) for name in feature_names]
        if all(np.isfinite(vec)):
            X_list.append(vec)
            y_list.append(f['alpha'])

    return np.array(X_list), np.array(y_list), feature_names


# ============================================================================
# HELPER FUNCTIONS FOR FEATURE IMPORTANCE
# ============================================================================

def _adaptive_cv_params(n_samples: int, method: str) -> Tuple[Optional[int], Optional[int]]:
    """Determine optimal K-fold and n_repeats based on sample size and method.

    Args:
        n_samples: Number of samples in dataset
        method: 'pfi' or 'cfi'

    Returns:
        (n_splits, n_repeats) or (None, None) if sample too small

    Rationale:
        - CFI requires larger training sets to learn conditional distributions p(x_j | x_{-j})
        - Smaller K → larger training sets per fold → better for CFI
        - PFI is more robust to smaller training sets
    """
    if method == 'cfi':
        # CFI needs larger training sets for conditional samplers
        if n_samples < 100:
            return None, None  # Skip CFI entirely
        elif n_samples < 120:
            return 3, 7  # ~67-80 samples per train fold
        elif n_samples < 150:
            return 4, 7  # ~75-112 samples per train fold
        else:
            return 5, 10  # ~80-120+ samples per train fold
    else:  # pfi
        if n_samples < 80:
            return 3, 5  # Minimum viable
        elif n_samples < 100:
            return 3, 7  # Conservative
        elif n_samples < 130:
            return 5, 10  # Standard for n~100-120
        else:
            return 5, 10  # Standard for larger samples


def _compute_derived_metrics(
    delta_mse_mean: np.ndarray,
    delta_mse_se: np.ndarray,
    mse_baseline: float
) -> Dict[str, np.ndarray]:
    """Compute ΔRMSE and % increase from ΔMSE using error propagation.

    Args:
        delta_mse_mean: Mean MSE increase per feature (shape: n_features)
        delta_mse_se: Standard error of MSE increase (shape: n_features)
        mse_baseline: Baseline MSE before permutation

    Returns:
        Dictionary with delta_rmse_mean, delta_rmse_se, pct_increase_mean, pct_increase_se

    Formulas:
        ΔRMSE = sqrt(MSE₀ + ΔMSE) - sqrt(MSE₀)
        SE(ΔRMSE) ≈ SE(ΔMSE) / (2 × sqrt(MSE₀ + ΔMSE))  [Taylor approximation]

        %Δ = 100 × ΔMSE / MSE₀
        SE(%Δ) = 100 × SE(ΔMSE) / MSE₀
    """
    # ΔRMSE computation
    delta_rmse_mean = np.sqrt(mse_baseline + delta_mse_mean) - np.sqrt(mse_baseline)
    # Error propagation for ΔRMSE (first-order Taylor approximation)
    delta_rmse_se = delta_mse_se / (2 * np.sqrt(mse_baseline + delta_mse_mean))

    # Percent increase computation
    pct_increase_mean = 100 * delta_mse_mean / mse_baseline
    pct_increase_se = 100 * delta_mse_se / mse_baseline

    return {
        'delta_rmse_mean': delta_rmse_mean,
        'delta_rmse_se': delta_rmse_se,
        'pct_increase_mean': pct_increase_mean,
        'pct_increase_se': pct_increase_se
    }


def compute_feature_importance_cv(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: List[str],
    method: str = 'pfi',
    n_splits: int = None,
    n_repeats_cv: int = None,
    n_repeats_perm: int = 20,
    random_state: int = 42
) -> Optional[Dict[str, Any]]:
    """
    Compute feature importance with cross-validation for PFI or CFI.

    Uses RepeatedKFold cross-validation to provide stable importance estimates
    and proper uncertainty quantification. Supports both:
      - PFI (Permutation Feature Importance): Marginal importance via naive shuffling
      - CFI (Conditional Feature Importance): Conditional importance accounting for correlations

    CV parameters (n_splits, n_repeats_cv) are auto-determined based on sample size
    if not provided. CFI uses smaller K (larger training sets) to ensure adequate
    samples for learning conditional distributions p(x_j | x_{-j}).

    Parameters
    ----------
    X : np.ndarray
        Feature matrix (n_samples, n_features)
    y : np.ndarray
        Target array (Murray exponent alpha)
    feature_names : List[str]
        Names of features
    method : str
        'pfi' for Permutation FI or 'cfi' for Conditional FI
    n_splits : int, optional
        Number of CV folds (auto-determined from sample size if None)
    n_repeats_cv : int, optional
        Number of CV repetitions (auto-determined from sample size if None)
    n_repeats_perm : int
        Number of permutation repeats per fold
    random_state : int
        Random seed for reproducibility

    Returns
    -------
    Dict with keys:
        'method': str - 'pfi' or 'cfi'
        'mse_baseline': float - Mean baseline MSE across folds
        'mse_baseline_se': float - Standard error of baseline MSE
        'rmse_baseline': float - sqrt(mse_baseline), in alpha units
        'feature_names': List[str] - Sorted by importance (descending)
        'delta_mse_mean': np.ndarray - Mean MSE increase per feature
        'delta_mse_se': np.ndarray - Standard error of MSE increase
        'delta_rmse_mean': np.ndarray - Mean RMSE increase per feature
        'delta_rmse_se': np.ndarray - Standard error of RMSE increase
        'pct_increase_mean': np.ndarray - Percent MSE increase per feature
        'pct_increase_se': np.ndarray - Standard error of percent increase
        'cv_r2_mean': float - Mean CV R-squared
        'cv_r2_se': float - Standard error of CV R-squared
        'n_samples': int
        'n_folds': int
        'n_repeats': int
        'n_total_splits': int
        'low_confidence': bool - True if n < 100

    Returns None if sample size insufficient or method requirements not met.
    """
    n_samples = len(y)

    # Determine adaptive CV parameters if not provided
    if n_splits is None or n_repeats_cv is None:
        auto_splits, auto_repeats = _adaptive_cv_params(n_samples, method)
        if auto_splits is None:
            # Sample too small for this method
            return None
        if n_splits is None:
            n_splits = auto_splits
        if n_repeats_cv is None:
            n_repeats_cv = auto_repeats

    # Setup cross-validation
    cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats_cv,
                       random_state=random_state)

    # Flag low confidence for small samples
    low_confidence = n_samples < 100

    # Storage for all folds
    mse_baselines = []
    r2_scores = []
    all_delta_mse = []  # Shape: (n_folds, n_features) for CFI or (n_folds, n_features, n_repeats_perm) for PFI

    # Branch based on method
    if method == 'pfi':
        # =====================================================================
        # PFI: Permutation Feature Importance (marginal)
        # =====================================================================
        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Train model
            rf = RandomForestRegressor(
                n_estimators=100, max_depth=10, min_samples_split=20,
                min_samples_leaf=10, random_state=random_state, n_jobs=-1
            )
            rf.fit(X_train, y_train)

            # Baseline metrics
            y_pred = rf.predict(X_test)
            mse_baselines.append(mean_squared_error(y_test, y_pred))
            r2_scores.append(rf.score(X_test, y_test))

            # Permutation importance (naive shuffling)
            perm_result = permutation_importance(
                rf, X_test, y_test,
                scoring='neg_mean_squared_error',
                n_repeats=n_repeats_perm,
                random_state=random_state + fold_idx,
                n_jobs=-1
            )

            # perm_result.importances shape: (n_features, n_repeats_perm)
            all_delta_mse.append(perm_result.importances)

    elif method == 'cfi':
        # =====================================================================
        # CFI: Conditional Feature Importance (respects correlations)
        # =====================================================================
        if not FIPPY_AVAILABLE:
            warnings.warn("CFI requires fippy library. Returning None.", UserWarning)
            return None

        try:
            import pandas as pd
            from fippy.samplers import SequentialSampler, ContUnivRFSampler
            from fippy.worker import FIComputeWorker
        except ImportError:
            warnings.warn("CFI requires fippy.samplers and fippy.worker. Returning None.", UserWarning)
            return None

        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Train model
            rf = RandomForestRegressor(
                n_estimators=100, max_depth=10, min_samples_split=20,
                min_samples_leaf=10, random_state=random_state, n_jobs=-1
            )
            rf.fit(X_train, y_train)

            # Baseline metrics
            y_pred = rf.predict(X_test)
            mse_baseline = mean_squared_error(y_test, y_pred)
            mse_baselines.append(mse_baseline)
            r2_scores.append(rf.score(X_test, y_test))

            # Build conditional samplers
            X_train_df = pd.DataFrame(X_train, columns=feature_names)
            y_train_series = pd.Series(y_train)

            sampler = SequentialSampler(X_train_df, ContUnivRFSampler)
            wrk = FIComputeWorker(rf, sampler, loss_fn=mean_squared_error)

            # Compute CFI
            X_test_df = pd.DataFrame(X_test, columns=feature_names)
            y_test_series = pd.Series(y_test)
            cfi_result = wrk.cfi(X_test_df, y_test_series)

            # Extract CFI means and stds
            cfi_means, cfi_stds = _extract_fi_arrays(cfi_result.fi_means_stds(), feature_names)

            # CFI returns importance as delta_MSE directly
            # Shape: (n_features,) - single value per feature per fold
            all_delta_mse.append(cfi_means)

    else:
        raise ValueError(f"method must be 'pfi' or 'cfi', got '{method}'")

    # =========================================================================
    # AGGREGATE RESULTS ACROSS FOLDS
    # =========================================================================

    # Baseline MSE and R²
    mse_baseline_mean = np.mean(mse_baselines)
    mse_baseline_se = np.std(mse_baselines, ddof=1) / np.sqrt(len(mse_baselines))
    rmse_baseline = np.sqrt(mse_baseline_mean)

    r2_mean = np.mean(r2_scores)
    r2_se = np.std(r2_scores, ddof=1) / np.sqrt(len(r2_scores))

    # Pool importance values across folds (and repeats for PFI)
    if method == 'pfi':
        # all_delta_mse is list of arrays, each shape (n_features, n_repeats_perm)
        # Stack to (n_folds, n_features, n_repeats_perm)
        all_importances = np.stack(all_delta_mse, axis=0)
        n_cv_folds = all_importances.shape[0]
        n_features = all_importances.shape[1]

        # Reshape to (n_features, n_folds * n_repeats_perm) for aggregation
        all_importances_flat = all_importances.transpose(1, 0, 2).reshape(n_features, -1)
        n_total_repeats = all_importances_flat.shape[1]

        # Aggregate: mean and SE across all folds and permutation repeats
        delta_mse_mean = np.mean(all_importances_flat, axis=1)
        delta_mse_se = np.std(all_importances_flat, axis=1, ddof=1) / np.sqrt(n_total_repeats)

    elif method == 'cfi':
        # all_delta_mse is list of arrays, each shape (n_features,)
        # Stack to (n_folds, n_features)
        all_importances = np.stack(all_delta_mse, axis=0)
        n_cv_folds = all_importances.shape[0]
        n_features = all_importances.shape[1]
        n_total_repeats = n_cv_folds

        # Aggregate: mean and SE across folds
        delta_mse_mean = np.mean(all_importances, axis=0)
        delta_mse_se = np.std(all_importances, axis=0, ddof=1) / np.sqrt(n_cv_folds)

    # Sort features by importance (descending)
    sorted_idx = np.argsort(delta_mse_mean)[::-1]
    sorted_names = [feature_names[i] for i in sorted_idx]
    delta_mse_mean = delta_mse_mean[sorted_idx]
    delta_mse_se = delta_mse_se[sorted_idx]

    # Compute derived metrics using helper function
    derived = _compute_derived_metrics(delta_mse_mean, delta_mse_se, mse_baseline_mean)

    return {
        'method': method,
        'mse_baseline': mse_baseline_mean,
        'mse_baseline_se': mse_baseline_se,
        'rmse_baseline': rmse_baseline,
        'feature_names': sorted_names,
        'delta_mse_mean': delta_mse_mean,
        'delta_mse_se': delta_mse_se,
        'delta_rmse_mean': derived['delta_rmse_mean'],
        'delta_rmse_se': derived['delta_rmse_se'],
        'pct_increase_mean': derived['pct_increase_mean'],
        'pct_increase_se': derived['pct_increase_se'],
        'cv_r2_mean': r2_mean,
        'cv_r2_se': r2_se,
        'n_samples': n_samples,
        'n_folds': n_cv_folds,
        'n_repeats': n_repeats_cv,
        'n_total_splits': n_total_repeats,
        'low_confidence': low_confidence
    }


# ============================================================================
# REUSABLE PLOTTING COMPONENTS
# ============================================================================

def create_hexbin_plot(ax, x: np.ndarray, y: np.ndarray, gridsize: int = 60,
                       cmap: str = 'plasma', xscale: str = None,
                       show_trend: bool = True, n_bins: int = 12,
                       trend_style: str = 'white'):
    """
    Create a hexbin plot with optional binned median trend line.

    Parameters
    ----------
    ax : matplotlib axis
    x, y : arrays
    gridsize : int
    cmap : str
    xscale : str or None
        'log' for log x-axis, None for linear
    show_trend : bool
        Whether to overlay binned median trend
    n_bins : int
        Bins for trend line
    trend_style : str
        'white' for standalone plots, 'cyan' for panel plots
    """
    hexbin = ax.hexbin(x, y, gridsize=gridsize, cmap=cmap, mincnt=1,
                       xscale=xscale if xscale else 'linear', edgecolors='none')

    if show_trend and len(x) > 0:
        centers, medians = compute_binned_medians(x, y, n_bins=n_bins,
                                                   log_scale=(xscale == 'log'), min_count=3)
        if centers:
            if trend_style == 'white':
                ax.plot(centers, medians, 'w-', lw=2.5, alpha=0.9, label='Binned median')
                ax.plot(centers, medians, 'k--', lw=1.5, alpha=0.7)
            else:  # cyan style for panels
                ax.plot(centers, medians, 'k-', lw=3, alpha=0.9, zorder=10)
                ax.plot(centers, medians, 'cyan', lw=1.5, alpha=0.9, zorder=11)

    ax.axhline(3.0, color='red', linestyle='--', lw=2 if trend_style == 'white' else 1.5,
               alpha=0.7, label='Classical Murray (α=3)' if trend_style == 'white' else None)

    return hexbin


def add_importance_bars(ax, sorted_names: List[str], sorted_means: np.ndarray,
                        sorted_stds: np.ndarray, test_score: float = None,
                        fontsize: int = 12, capsize: int = 5):
    """Add horizontal bar chart for feature importance."""
    y_pos = np.arange(len(sorted_names))
    ax.barh(y_pos, sorted_means, xerr=sorted_stds, color='steelblue', alpha=0.7,
            edgecolor='black', linewidth=1.2, capsize=capsize)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_names, fontsize=fontsize, weight='bold')
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')


# ============================================================================
# CONDITIONAL FEATURE IMPORTANCE (CFI) HELPERS
# ============================================================================
#
# CFI uses a conditional sampler (fippy library) so that when a feature is
# permuted, the replacement values are drawn conditional on the remaining
# features.  This keeps the joint feature distribution realistic and avoids
# the correlation-induced importance misattribution that naive (marginal)
# permutation importance (PFI) suffers from.
#
# Interpretation key:
#   - PFI high, CFI drops a lot  -> feature was predictive mainly through
#     redundancy with correlated features
#   - CFI stays high             -> feature carries unique predictive
#     information beyond the correlated set
#   - CFI rises relative to PFI  -> naive shuffling created unrealistic
#     feature combos that masked the feature's real contribution
# ============================================================================

def _extract_fi_arrays(fi_result, feature_names: List[str]
                       ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Robustly extract mean and std arrays from a fippy fi_means_stds() result.

    fippy's fi_means_stds() returns a pandas object indexed by feature name.
    The internal layout can vary (DataFrame with 'std' column, or other
    structures depending on version).  This helper tries multiple access
    patterns so the code doesn't break if fippy changes its return format.

    Parameters
    ----------
    fi_result : pandas DataFrame or Series
        Output of Explanation.fi_means_stds().
    feature_names : list of str
        Feature names in the desired output order.

    Returns
    -------
    means, stds : np.ndarray, np.ndarray
        Arrays of length len(feature_names), aligned to feature_names order.
    """
    means = np.zeros(len(feature_names))
    stds = np.zeros(len(feature_names))

    for i, fn in enumerate(feature_names):
        row = fi_result.loc[fn]

        # Case 1: row is a Series (feature has multiple columns like 'mean', 'std')
        if hasattr(row, 'index') and 'std' in row.index:
            # The non-'std' entries are the mean(s); take the first one
            mean_vals = row.drop('std')
            means[i] = float(mean_vals.iloc[0]) if len(mean_vals) > 0 else float(mean_vals)
            stds[i] = float(row['std'])
        # Case 2: row is a scalar (only mean, no std column)
        elif np.isscalar(row):
            means[i] = float(row)
            stds[i] = 0.0
        # Case 3: row is a Series but without 'std' label
        else:
            means[i] = float(row.iloc[0]) if hasattr(row, 'iloc') else float(row)
            stds[i] = 0.0

    return means, stds


def plot_cfi_pfi_heatmap(pfi_sample_results: Dict[str, Dict],
                         cfi_sample_results: Dict[str, Dict],
                         output_path: str,
                         feature_order: List[str] = None,
                         metric: str = 'delta_mse'):
    """
    Plot side-by-side PFI and CFI heatmaps across all samples.

    Left panel:  PFI importance (features x samples)
    Right panel: CFI importance (features x samples)

    The comparison reveals where correlation distorts naive importance:
      - Feature bright in PFI but dim in CFI -> importance was from redundancy
      - Feature bright in both -> genuinely unique predictive signal
      - Feature brighter in CFI than PFI -> naive shuffling masked its role

    Parameters
    ----------
    pfi_sample_results : dict
        Maps sample_name -> pfi_result dict from compute_feature_importance_cv(..., method='pfi').
    cfi_sample_results : dict
        Maps sample_name -> cfi_result dict from compute_feature_importance_cv(..., method='cfi').
    output_path : str
        Where to save the figure.
    feature_order : list of str, optional
        Fixed ordering of features on the y-axis. If None, uses mean CFI
        across samples (descending) so the most important features are on top.
    metric : str
        Which metric to plot: 'delta_mse', 'delta_rmse', or 'pct_increase' (default: 'delta_mse')
    """
    # Filter to samples that have valid results in BOTH PFI and CFI
    valid_samples = set(k for k, v in pfi_sample_results.items() if v is not None) & \
                    set(k for k, v in cfi_sample_results.items() if v is not None)

    if len(valid_samples) == 0:
        print("    No valid PFI+CFI results to plot heatmap")
        return

    sample_names = sorted(list(valid_samples))

    # Use the feature names from the first valid result (they are identical)
    all_feature_names = list(cfi_sample_results[sample_names[0]]['feature_names'])

    # Determine metric key
    if metric == 'pct_increase':
        metric_key = 'pct_increase_mean'
    elif metric == 'delta_rmse':
        metric_key = 'delta_rmse_mean'
    else:  # delta_mse
        metric_key = 'delta_mse_mean'

    # Determine feature ordering
    # Default: sort by mean CFI across samples, descending
    if feature_order is None:
        mean_cfi_per_feature = {}
        for fn in all_feature_names:
            vals = []
            for sname in sample_names:
                cfi_res = cfi_sample_results[sname]
                fnames = list(cfi_res['feature_names'])
                idx = fnames.index(fn)
                vals.append(cfi_res[metric_key][idx])
            mean_cfi_per_feature[fn] = np.mean(vals)
        feature_order = sorted(all_feature_names,
                               key=lambda fn: mean_cfi_per_feature[fn],
                               reverse=True)

    n_features = len(feature_order)
    n_samples = len(sample_names)

    # Build heatmap matrices (features x samples)
    pfi_matrix = np.zeros((n_features, n_samples))
    cfi_matrix = np.zeros((n_features, n_samples))

    for j, sname in enumerate(sample_names):
        pfi_res = pfi_sample_results[sname]
        cfi_res = cfi_sample_results[sname]
        pfi_fnames = list(pfi_res['feature_names'])
        cfi_fnames = list(cfi_res['feature_names'])

        for i, fn in enumerate(feature_order):
            pfi_idx = pfi_fnames.index(fn)
            cfi_idx = cfi_fnames.index(fn)
            pfi_matrix[i, j] = pfi_res[metric_key][pfi_idx]
            cfi_matrix[i, j] = cfi_res[metric_key][cfi_idx]

    # --- Determine shared color scale so PFI and CFI panels are comparable ---
    vmin = min(pfi_matrix.min(), cfi_matrix.min())
    vmax = max(pfi_matrix.max(), cfi_matrix.max())
    # Ensure non-negative floor (importance should be >= 0 in expectation)
    vmin = min(vmin, 0)

    # --- Create figure with two side-by-side heatmaps ---
    fig, (ax_pfi, ax_cfi) = plt.subplots(1, 2, figsize=(14, 6),
                                          sharey=True)

    # PFI heatmap (left panel)
    im_pfi = ax_pfi.imshow(pfi_matrix, aspect='auto', cmap='YlOrRd',
                            vmin=vmin, vmax=vmax, interpolation='nearest')
    ax_pfi.set_xticks(np.arange(n_samples))
    ax_pfi.set_xticklabels(sample_names, rotation=45, ha='right', fontsize=10)
    ax_pfi.set_yticks(np.arange(n_features))
    ax_pfi.set_yticklabels(feature_order, fontsize=11, weight='bold')
    ax_pfi.set_title('PFI (Marginal Permutation)', fontsize=13, weight='bold',
                     pad=10)

    # Annotate PFI cells with numeric values for readability
    for i in range(n_features):
        for j in range(n_samples):
            val = pfi_matrix[i, j]
            # Use white text on dark cells, black on light cells
            text_color = 'white' if val > (vmin + vmax) / 2 else 'black'
            ax_pfi.text(j, i, f'{val:.3f}', ha='center', va='center',
                        fontsize=8, color=text_color)

    # CFI heatmap (right panel)
    im_cfi = ax_cfi.imshow(cfi_matrix, aspect='auto', cmap='YlOrRd',
                            vmin=vmin, vmax=vmax, interpolation='nearest')
    ax_cfi.set_xticks(np.arange(n_samples))
    ax_cfi.set_xticklabels(sample_names, rotation=45, ha='right', fontsize=10)
    ax_cfi.set_title('CFI (Conditional Permutation)', fontsize=13,
                     weight='bold', pad=10)

    # Annotate CFI cells with numeric values
    for i in range(n_features):
        for j in range(n_samples):
            val = cfi_matrix[i, j]
            text_color = 'white' if val > (vmin + vmax) / 2 else 'black'
            ax_cfi.text(j, i, f'{val:.3f}', ha='center', va='center',
                        fontsize=8, color=text_color)

    # Shared colorbar on the right edge
    fig.colorbar(im_cfi, ax=[ax_pfi, ax_cfi], label='Mean increase in MSE',
                 shrink=0.8, pad=0.02)

    fig.suptitle(
        'Conditional vs Marginal Feature Importance Across Samples\n'
        'PFI high but CFI low → importance from correlated redundancy; '
        'CFI high → unique predictive signal',
        fontsize=12, weight='bold', y=1.02
    )

    fig.tight_layout(rect=[0, 0, 0.92, 0.95])
    if safe_save_figure(fig, output_path, dpi=FIGURE_DPI):
        print(f"  ✓ Saved: {output_path}")
    plt.close()


def plot_cfi_pfi_bars(pfi_result: Dict, cfi_result: Dict, output_path: str,
                      sample_name: str = None, metric: str = 'delta_mse'):
    """
    Plot side-by-side horizontal bar chart comparing PFI and CFI for one dataset.

    This is the single-sample (or pooled) analog of the heatmap. Each feature
    gets two bars: one for PFI (marginal) and one for CFI (conditional), so
    you can directly see which features lose or gain importance when
    correlations are properly accounted for.

    Parameters
    ----------
    pfi_result : dict
        Output from compute_feature_importance_cv(..., method='pfi').
    cfi_result : dict
        Output from compute_feature_importance_cv(..., method='cfi').
    output_path : str
        Where to save the figure.
    sample_name : str, optional
        Name for the title (None = pooled across all samples).
    metric : str
        Which metric to plot: 'delta_mse', 'delta_rmse', or 'pct_increase' (default: 'delta_mse')
    """
    if pfi_result is None or cfi_result is None:
        print("    Skipping PFI vs CFI comparison (one or both results unavailable)")
        return

    title_prefix = f"{sample_name}: " if sample_name else ""

    # Extract importance values based on metric
    if metric == 'pct_increase':
        pfi_values = pfi_result['pct_increase_mean']
        pfi_errors = pfi_result['pct_increase_se'] * 2
        cfi_values = cfi_result['pct_increase_mean']
        cfi_errors = cfi_result['pct_increase_se'] * 2
        xlabel = '% Increase in MSE when feature shuffled'
    elif metric == 'delta_rmse':
        pfi_values = pfi_result['delta_rmse_mean']
        pfi_errors = pfi_result['delta_rmse_se'] * 2
        cfi_values = cfi_result['delta_rmse_mean']
        cfi_errors = cfi_result['delta_rmse_se'] * 2
        xlabel = 'Increase in RMSE (α units) when feature shuffled'
    else:  # delta_mse
        pfi_values = pfi_result['delta_mse_mean']
        pfi_errors = pfi_result['delta_mse_se'] * 2
        cfi_values = cfi_result['delta_mse_mean']
        cfi_errors = cfi_result['delta_mse_se'] * 2
        xlabel = 'Increase in MSE when feature shuffled'

    feature_names = cfi_result['feature_names']
    n_samples = cfi_result['n_samples']
    cv_r2_mean = cfi_result['cv_r2_mean']
    cv_r2_se = cfi_result['cv_r2_se']
    rmse_baseline = cfi_result['rmse_baseline']

    # Sort features by CFI importance (descending) so the most
    # conditionally important feature is at the top of the chart
    sorted_idx = np.argsort(cfi_values)  # ascending; plot bottom-to-top
    sorted_names = [feature_names[i] for i in sorted_idx]
    sorted_pfi_values = pfi_values[sorted_idx]
    sorted_pfi_errors = pfi_errors[sorted_idx]
    sorted_cfi_values = cfi_values[sorted_idx]
    sorted_cfi_errors = cfi_errors[sorted_idx]

    n_features = len(sorted_names)
    y_pos = np.arange(n_features)
    bar_height = 0.35  # half-height so PFI and CFI bars don't overlap

    fig, ax = plt.subplots(figsize=(10, 6))

    # PFI bars (offset up)
    ax.barh(y_pos + bar_height / 2, sorted_pfi_values, height=bar_height,
            xerr=sorted_pfi_errors, color='steelblue', alpha=0.7,
            edgecolor='black', linewidth=0.8, capsize=4, label='PFI (Permutation)')

    # CFI bars (offset down)
    ax.barh(y_pos - bar_height / 2, sorted_cfi_values, height=bar_height,
            xerr=sorted_cfi_errors, color='darkorange', alpha=0.7,
            edgecolor='black', linewidth=0.8, capsize=4, label='CFI (Conditional)')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_names, fontsize=11, weight='bold')
    ax.set_xlabel(xlabel, fontsize=11, weight='bold')
    ax.set_title(
        f'{title_prefix}PFI vs CFI - Feature Importance for Murray Exponent α\n'
        f'Random Forest with CV (n={n_samples:,}, CV R²={cv_r2_mean:.3f})',
        fontsize=13, weight='bold', pad=15
    )
    ax.legend(loc='lower right', fontsize=10, framealpha=0.9)
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')

    # Annotation box with model performance (removed n_folds per user request)
    ax.text(0.98, 0.02,
            f'Baseline RMSE: {rmse_baseline:.3f}\n'
            f'CV R²: {cv_r2_mean:.3f} ± {cv_r2_se:.3f}\n'
            f'n = {n_samples:,}',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig.tight_layout()
    if safe_save_figure(fig, output_path, dpi=FIGURE_DPI):
        print(f"{'    ' if sample_name else ''}✓ Saved: "
              f"{os.path.basename(output_path) if sample_name else output_path}")
    plt.close()


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_alpha_vs_parent_radius_hexbin(all_features: List[Dict], output_path: str,
                                       sample_name: str = None):
    """Plot 2D hexbin: Alpha vs Parent Radius (works for pooled or single sample)."""
    title_prefix = f"{sample_name}: " if sample_name else ""
    print(f"\n{'='*70}")
    print(f"PLOT: {title_prefix}Alpha vs Parent Radius (Hexbin)")
    print("="*70)

    r_parent, alpha = extract_plot_data(all_features, 'r_parent')

    if len(r_parent) == 0:
        print(f"    No valid data" + (f" for {sample_name}" if sample_name else ""))
        return

    print(f"Valid bifurcations: {len(r_parent)}")

    fig = plt.figure(figsize=(12, 8))
    ax = plt.subplot(111)

    hexbin = create_hexbin_plot(ax, r_parent, alpha, gridsize=60 if sample_name else 80,
                                cmap='plasma', xscale='log')

    ax.set_xlabel('Parent radius', fontsize=13, weight='bold')
    ax.set_ylabel('Fitted Murray exponent α', fontsize=13, weight='bold')
    ax.set_title(f'{title_prefix}Murray Exponent vs Vessel Scale\n'
                 f'Tests for scale-dependent exponent variation',
                 fontsize=14, weight='bold', pad=15)
    ax.set_ylim(0.5, 5.5)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9)

    plt.colorbar(hexbin, ax=ax, label='Count').ax.tick_params(labelsize=10)

    ax.text(0.02, 0.98,
            f'n = {len(r_parent):,} bifurcations\n'
            f'Median α = {np.median(alpha):.2f}\n'
            f'Parent r range: {r_parent.min():.1f}–{r_parent.max():.1f}',
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    fig.tight_layout()
    if safe_save_figure(fig, output_path, dpi=FIGURE_DPI):
        print(f"{'    ' if sample_name else ''}✓ Saved: {os.path.basename(output_path) if sample_name else output_path}")
    plt.close()


def plot_alpha_vs_asymmetry_hexbin(all_features: List[Dict], output_path: str,
                                   sample_name: str = None):
    """Plot 2D hexbin: Alpha vs Daughter Asymmetry (works for pooled or single sample)."""
    title_prefix = f"{sample_name}: " if sample_name else ""
    print(f"\n{'='*70}")
    print(f"PLOT: {title_prefix}Alpha vs Daughter Asymmetry (Hexbin)")
    print("="*70)

    asymmetry, alpha = extract_plot_data(all_features, 'asymmetry')

    if len(asymmetry) == 0:
        print(f"    No valid data" + (f" for {sample_name}" if sample_name else ""))
        return

    print(f"Valid bifurcations: {len(asymmetry)}")

    fig = plt.figure(figsize=(12, 8))
    ax = plt.subplot(111)

    hexbin = create_hexbin_plot(ax, asymmetry, alpha, gridsize=60 if sample_name else 80,
                                cmap='viridis', xscale=None)

    ax.set_xlabel('Daughter asymmetry (min/max radius)', fontsize=13, weight='bold')
    ax.set_ylabel('Fitted Murray exponent α', fontsize=13, weight='bold')
    ax.set_title(f'{title_prefix}Murray Exponent vs Branching Geometry\n'
                 f'Tests for asymmetry-dependent exponent variation',
                 fontsize=14, weight='bold', pad=15)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 5.5)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9)

    plt.colorbar(hexbin, ax=ax, label='Count').ax.tick_params(labelsize=10)

    info_text = (f'n = {len(asymmetry):,} bifurcations\n'
                 f'Median α = {np.median(alpha):.2f}\n'
                 f'Median asymmetry = {np.median(asymmetry):.3f}')
    if not sample_name:
        info_text += f'\nSymmetric (>0.8): {np.sum(asymmetry > 0.8)/len(asymmetry):.1%}'

    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    fig.tight_layout()
    if safe_save_figure(fig, output_path, dpi=FIGURE_DPI):
        print(f"{'    ' if sample_name else ''}✓ Saved: {os.path.basename(output_path) if sample_name else output_path}")
    plt.close()


def plot_feature_importance(
    result: Dict[str, Any],
    output_path: str,
    sample_name: str = None,
    metric: str = 'pct_increase'
):
    """
    Plot feature importance (PFI or CFI) with cross-validation and interpretable metrics.

    Parameters
    ----------
    result : Dict[str, Any]
        Result dict from compute_feature_importance_cv() containing importance values,
        baseline metrics, and metadata
    output_path : str
        Where to save figure
    sample_name : str, optional
        Sample name for title
    metric : str
        Which importance metric to display:
        - 'pct_increase': Percent MSE increase (default, best for cross-sample comparison)
        - 'delta_rmse': RMSE increase in alpha units
        - 'delta_mse': Raw MSE increase
    """
    if result is None:
        print("    No result to plot (sample too small or method unavailable)")
        return

    method = result.get('method', 'pfi')
    method_label = 'PFI (Permutation)' if method == 'pfi' else 'CFI (Conditional)'

    title_prefix = f"{sample_name}: " if sample_name else ""
    print(f"\n{'='*70}")
    print(f"PLOT: {title_prefix}{method_label} - {metric}")
    print("="*70)

    # Select metric to plot
    if metric == 'pct_increase':
        values = result['pct_increase_mean']
        errors = result['pct_increase_se'] * 2  # 2×SE for ~95% CI
        xlabel = '% Increase in MSE when feature shuffled'
    elif metric == 'delta_rmse':
        values = result['delta_rmse_mean']
        errors = result['delta_rmse_se'] * 2
        xlabel = 'Increase in RMSE (α units) when feature shuffled'
    else:  # delta_mse
        values = result['delta_mse_mean']
        errors = result['delta_mse_se'] * 2
        xlabel = 'Increase in MSE when feature shuffled'

    print(f"CV R²: {result['cv_r2_mean']:.3f} ± {result['cv_r2_se']:.3f}")
    print(f"Baseline RMSE: {result['rmse_baseline']:.3f} (MSE: {result['mse_baseline']:.4f})")
    if result['low_confidence']:
        print(f"    ⚠ LOW CONFIDENCE: Small sample size")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Horizontal bar chart
    y_pos = np.arange(len(result['feature_names']))
    ax.barh(y_pos, values, xerr=errors, color='steelblue', alpha=0.7,
            edgecolor='black', linewidth=1.2, capsize=5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(result['feature_names'], fontsize=12, weight='bold')
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')

    ax.set_xlabel(xlabel, fontsize=12, weight='bold')

    # Title with sample info
    confidence_note = " [LOW CONFIDENCE]" if result['low_confidence'] else ""
    ax.set_title(
        f'{title_prefix}{method_label} - Feature Importance for Murray Exponent α{confidence_note}\n'
        f'Random Forest with CV (n={result["n_samples"]:,}, CV R²={result["cv_r2_mean"]:.3f})',
        fontsize=13, weight='bold', pad=15
    )

    # Annotation box with baseline metrics (removed n_folds info per user request)
    annotation_text = (
        f'Baseline RMSE: {result["rmse_baseline"]:.3f}\n'
        f'Baseline MSE: {result["mse_baseline"]:.4f}\n'
        f'CV R²: {result["cv_r2_mean"]:.3f} ± {result["cv_r2_se"]:.3f}\n'
        f'n = {result["n_samples"]:,}'
    )

    bbox_color = 'mistyrose' if result['low_confidence'] else 'lightyellow'
    ax.text(0.98, 0.02, annotation_text,
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            family='monospace',
            bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.9))

    # Add warning banner if low confidence
    if result['low_confidence']:
        ax.text(0.5, 0.98, '⚠ LOW CONFIDENCE: Small sample size',
                transform=ax.transAxes, fontsize=10, ha='center', va='top',
                color='darkred', weight='bold',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    fig.tight_layout()
    if safe_save_figure(fig, output_path, dpi=FIGURE_DPI):
        print(f"{'    ' if sample_name else ''}✓ Saved: {os.path.basename(output_path) if sample_name else output_path}")
    plt.close()

    return result  # Return result for further analysis if needed


def plot_three_panel_summary(all_features: List[Dict], output_path: str,
                             sample_name: str = None):
    """Create three-panel summary figure (works for pooled or single sample)."""
    title_prefix = f"{sample_name}: " if sample_name else ""
    print(f"\n{'='*70}")
    print(f"CREATING {title_prefix.upper()}THREE-PANEL SUMMARY FIGURE")
    print("="*70)

    r_parent, alpha_A = extract_plot_data(all_features, 'r_parent')
    asymmetry, alpha_B = extract_plot_data(all_features, 'asymmetry')

    fig = plt.figure(figsize=(18, 6))

    # Panel A: Alpha vs Parent Radius
    ax_A = plt.subplot(1, 3, 1)
    if len(r_parent) > 0:
        hexbin_A = create_hexbin_plot(ax_A, r_parent, alpha_A, gridsize=50,
                                      cmap='plasma', xscale='log', trend_style='cyan')
        plt.colorbar(hexbin_A, ax=ax_A).set_label('Count', fontsize=9)
    ax_A.set_xlabel('Parent radius', fontsize=11, weight='bold')
    ax_A.set_ylabel('Murray exponent α', fontsize=11, weight='bold')
    ax_A.set_title('A. Scale Dependence', fontsize=12, weight='bold', pad=10)
    ax_A.set_ylim(0.5, 5.5)
    ax_A.grid(True, alpha=0.25)

    # Panel B: Alpha vs Asymmetry
    ax_B = plt.subplot(1, 3, 2)
    if len(asymmetry) > 0:
        hexbin_B = create_hexbin_plot(ax_B, asymmetry, alpha_B, gridsize=50,
                                      cmap='viridis', xscale=None, trend_style='cyan')
        plt.colorbar(hexbin_B, ax=ax_B).set_label('Count', fontsize=9)
    ax_B.set_xlabel('Daughter asymmetry', fontsize=11, weight='bold')
    ax_B.set_ylabel('Murray exponent α', fontsize=11, weight='bold')
    ax_B.set_title('B. Geometry Dependence', fontsize=12, weight='bold', pad=10)
    ax_B.set_xlim(0, 1)
    ax_B.set_ylim(0.5, 5.5)
    ax_B.grid(True, alpha=0.25)

    # Panel C: Permutation Importance
    ax_C = plt.subplot(1, 3, 3)
    short_names = ['r_parent', 'asymmetry', 'angle', 'tortuosity', 'taper']

    if SKLEARN_AVAILABLE:
        X, y, _ = prepare_regression_data(all_features)
        if len(y) >= 100:
            sorted_names, sorted_means, sorted_stds, _, test_score = \
                train_and_get_importance(X, y, short_names, n_repeats=15)
            add_importance_bars(ax_C, sorted_names, sorted_means, sorted_stds, fontsize=10, capsize=4)
            ax_C.set_xlabel('Importance', fontsize=11, weight='bold')
            ax_C.set_title(f'C. Feature Importance\n(test R²={test_score:.3f})',
                          fontsize=12, weight='bold', pad=10)
        else:
            ax_C.text(0.5, 0.5, f'n={len(y)} < 100\nInsufficient data',
                     ha='center', va='center', transform=ax_C.transAxes, fontsize=11)
            ax_C.set_title('C. Feature Importance', fontsize=12, weight='bold', pad=10)
    else:
        ax_C.text(0.5, 0.5, 'scikit-learn\nnot available',
                 ha='center', va='center', transform=ax_C.transAxes, fontsize=12, style='italic')
        ax_C.set_title('C. Feature Importance', fontsize=12, weight='bold', pad=10)

    # Overall title
    n_bif = len(all_features)
    subtitle = f'{n_bif:,} bifurcations' + ('' if sample_name else ' across multiple species')
    fig.suptitle(f'{title_prefix}What Drives Murray Exponent Variation?\n{subtitle}',
                 fontsize=15, weight='bold', y=0.98)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    if safe_save_figure(fig, output_path, dpi=FIGURE_DPI):
        print(f"{'    ' if sample_name else ''}✓ Saved: {os.path.basename(output_path) if sample_name else output_path}")
    plt.close()


# ============================================================================
# INDIVIDUAL SAMPLE ANALYSIS
# ============================================================================

def generate_individual_sample_plots(sample_features: List[Dict], sample_name: str,
                                     output_dir: str):
    """
    Generate diagnostic plots for a single sample.

    Returns (pfi_result, cfi_result) tuple so that main() can reuse them
    for the cross-sample heatmap without recomputing.
    """
    print(f"\n{'='*70}")
    print(f"Generating individual plots for: {sample_name}")
    print(f"{'='*70}")
    print(f"Number of bifurcations: {len(sample_features)}")

    sample_dir = os.path.join(output_dir, sample_name.replace(' ', '_'))
    os.makedirs(sample_dir, exist_ok=True)
    print(f"Output directory: {sample_dir}")

    if len(sample_features) < 10:
        print(f"WARNING: Only {len(sample_features)} bifurcations - skipping plots")
        return None

    safe_name = sample_name.replace(' ', '_').replace('/', '_')

    # All plots use the unified functions with sample_name parameter
    plot_alpha_vs_parent_radius_hexbin(
        sample_features,
        os.path.join(sample_dir, f"{safe_name}_1_alpha_vs_radius.{FIGURE_FORMAT}"),
        sample_name=sample_name
    )

    plot_alpha_vs_asymmetry_hexbin(
        sample_features,
        os.path.join(sample_dir, f"{safe_name}_2_alpha_vs_asymmetry.{FIGURE_FORMAT}"),
        sample_name=sample_name
    )

    # Feature importance plots - compute PFI and CFI with cross-validation
    pfi_result = None
    cfi_result = None

    # Prepare regression data
    if SKLEARN_AVAILABLE and len(sample_features) >= MIN_SAMPLES_PFI:
        X, y, feature_names = prepare_regression_data(sample_features)

        # Compute PFI (Permutation Feature Importance)
        print(f"\n  Computing PFI (Permutation Feature Importance) for {sample_name}...")
        try:
            pfi_result = compute_feature_importance_cv(
                X, y, feature_names,
                method='pfi',
                n_repeats_perm=IMPORTANCE_PERM_REPEATS
            )
            if pfi_result is not None:
                plot_feature_importance(
                    pfi_result,
                    os.path.join(sample_dir, f"{safe_name}_3_pfi.{FIGURE_FORMAT}"),
                    sample_name=sample_name,
                    metric=IMPORTANCE_METRIC
                )
        except Exception as e:
            print(f"  ERROR computing PFI for {sample_name}: {str(e)}")
            import traceback
            traceback.print_exc()

        # Compute CFI (Conditional Feature Importance)
        if FIPPY_AVAILABLE and len(sample_features) >= MIN_SAMPLES_CFI:
            print(f"\n  Computing CFI (Conditional Feature Importance) for {sample_name}...")
            try:
                cfi_result = compute_feature_importance_cv(
                    X, y, feature_names,
                    method='cfi',
                    n_repeats_perm=IMPORTANCE_PERM_REPEATS
                )
                if cfi_result is not None:
                    # Standalone CFI plot
                    plot_feature_importance(
                        cfi_result,
                        os.path.join(sample_dir, f"{safe_name}_4_cfi.{FIGURE_FORMAT}"),
                        sample_name=sample_name,
                        metric=IMPORTANCE_METRIC
                    )

                    # PFI vs CFI comparison plot
                    if pfi_result is not None:
                        plot_cfi_pfi_bars(
                            pfi_result, cfi_result,
                            os.path.join(sample_dir, f"{safe_name}_5_cfi_vs_pfi.{FIGURE_FORMAT}"),
                            sample_name=sample_name,
                            metric=IMPORTANCE_METRIC
                        )
            except Exception as e:
                print(f"  ERROR computing CFI for {sample_name}: {str(e)}")
                import traceback
                traceback.print_exc()
        elif not FIPPY_AVAILABLE:
            print(f"  Skipping CFI (fippy not available)")
        else:
            print(f"  Skipping CFI (need ≥{MIN_SAMPLES_CFI} bifurcations, have {len(sample_features)})")
    else:
        print(f"  Skipping feature importance (need ≥{MIN_SAMPLES_PFI} bifurcations, have {len(sample_features)})")

    plot_three_panel_summary(
        sample_features,
        os.path.join(sample_dir, f"{safe_name}_summary_3panel.{FIGURE_FORMAT}"),
        sample_name=sample_name
    )

    print(f"✓ Completed plots for {sample_name}")
    return pfi_result, cfi_result


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def main():
    """Main execution function."""
    print("\n" + "="*100)
    print("MURRAY EXPONENT DRIVER ANALYSIS")
    print("="*100)
    print(f"\nAnalyzing {len(SKELETON_FILES)} samples:")
    for name in SKELETON_FILES.keys():
        print(f"  - {name}")

    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_DIR}")

    # Collect bifurcation features from all samples
    all_features = []
    sample_features_dict = {}  # Store features per sample for individual plots

    for sample_name, skeleton_path in SKELETON_FILES.items():
        print(f"\n{'='*70}")
        print(f"Processing: {sample_name}")
        print(f"{'='*70}")

        try:
            skeleton = load_skeleton(skeleton_path)
            voxel_spacing = VOXEL_SPACING[sample_name]

            features = extract_bifurcation_features(skeleton, voxel_spacing, sample_name)
            all_features.extend(features)
            sample_features_dict[sample_name] = features  # Store for individual plots

            print(f"  Collected {len(features)} bifurcations from {sample_name}")

        except Exception as e:
            print(f"  ERROR processing {sample_name}: {str(e)}")
            import traceback
            traceback.print_exc()
            sample_features_dict[sample_name] = []  # Empty list on error

    print(f"\n{'='*70}")
    print(f"TOTAL BIFURCATIONS COLLECTED: {len(all_features):,}")
    print(f"{'='*70}")

    if len(all_features) == 0:
        print("ERROR: No bifurcations collected. Exiting.")
        return

    # Summary statistics
    valid_alpha = [f['alpha'] for f in all_features if np.isfinite(f['alpha'])]
    print(f"\nValid alpha values: {len(valid_alpha):,}")
    if len(valid_alpha) > 0:
        print(f"Alpha range: {np.min(valid_alpha):.2f} - {np.max(valid_alpha):.2f}")
        print(f"Alpha median: {np.median(valid_alpha):.2f}")
        print(f"Alpha mean: {np.mean(valid_alpha):.2f} ± {np.std(valid_alpha):.2f}")

    # Generate plots
    print(f"\n{'='*70}")
    print("GENERATING DIAGNOSTIC PLOTS")
    print(f"{'='*70}")

    # Plot 1: Alpha vs Parent Radius
    plot1_path = os.path.join(OUTPUT_DIR, f"murray_drivers_1_alpha_vs_radius.{FIGURE_FORMAT}")
    try:
        plot_alpha_vs_parent_radius_hexbin(all_features, plot1_path)
    except Exception as e:
        print(f"ERROR generating Plot 1: {str(e)}")
        import traceback
        traceback.print_exc()

    # Plot 2: Alpha vs Asymmetry
    plot2_path = os.path.join(OUTPUT_DIR, f"murray_drivers_2_alpha_vs_asymmetry.{FIGURE_FORMAT}")
    try:
        plot_alpha_vs_asymmetry_hexbin(all_features, plot2_path)
    except Exception as e:
        print(f"ERROR generating Plot 2: {str(e)}")
        import traceback
        traceback.print_exc()

    # Three-panel summary
    summary_path = os.path.join(OUTPUT_DIR, f"murray_drivers_summary_3panel.{FIGURE_FORMAT}")
    try:
        plot_three_panel_summary(all_features, summary_path)
    except Exception as e:
        print(f"ERROR generating summary figure: {str(e)}")
        import traceback
        traceback.print_exc()

    # ========================================================================
    # FEATURE IMPORTANCE - POOLED ACROSS ALL SAMPLES
    # ========================================================================
    print(f"\n{'='*70}")
    print("FEATURE IMPORTANCE - ALL SAMPLES POOLED")
    print(f"{'='*70}")

    pooled_pfi_result = None
    pooled_cfi_result = None

    if SKLEARN_AVAILABLE and len(all_features) >= MIN_SAMPLES_PFI:
        X, y, feature_names = prepare_regression_data(all_features)

        # Plot 3: PFI (Permutation Feature Importance)
        plot3_path = os.path.join(OUTPUT_DIR, f"murray_drivers_3_pfi.{FIGURE_FORMAT}")
        try:
            pooled_pfi_result = compute_feature_importance_cv(
                X, y, feature_names,
                method='pfi',
                n_repeats_perm=IMPORTANCE_PERM_REPEATS
            )
            if pooled_pfi_result is not None:
                plot_feature_importance(
                    pooled_pfi_result, plot3_path,
                    metric=IMPORTANCE_METRIC
                )
        except Exception as e:
            print(f"ERROR generating pooled PFI plot: {str(e)}")
            import traceback
            traceback.print_exc()

        # Plot 4: CFI (Conditional Feature Importance)
        if FIPPY_AVAILABLE and len(all_features) >= MIN_SAMPLES_CFI:
            plot4_path = os.path.join(OUTPUT_DIR, f"murray_drivers_4_cfi.{FIGURE_FORMAT}")
            try:
                pooled_cfi_result = compute_feature_importance_cv(
                    X, y, feature_names,
                    method='cfi',
                    n_repeats_perm=IMPORTANCE_PERM_REPEATS
                )
                if pooled_cfi_result is not None:
                    plot_feature_importance(
                        pooled_cfi_result, plot4_path,
                        metric=IMPORTANCE_METRIC
                    )

                    # Plot 5: PFI vs CFI comparison
                    if pooled_pfi_result is not None:
                        plot5_path = os.path.join(OUTPUT_DIR, f"murray_drivers_5_cfi_vs_pfi.{FIGURE_FORMAT}")
                        plot_cfi_pfi_bars(
                            pooled_pfi_result, pooled_cfi_result,
                            plot5_path,
                            metric=IMPORTANCE_METRIC
                        )
            except Exception as e:
                print(f"ERROR generating pooled CFI plot: {str(e)}")
                import traceback
                traceback.print_exc()
        elif not FIPPY_AVAILABLE:
            print("  Skipping pooled CFI (fippy not available)")
        else:
            print(f"  Skipping pooled CFI (need ≥{MIN_SAMPLES_CFI} bifurcations)")
    else:
        print(f"  Skipping pooled feature importance (need ≥{MIN_SAMPLES_PFI} bifurcations)")

    # ========================================================================
    # INDIVIDUAL SAMPLE PLOTS
    # ========================================================================
    print(f"\n{'='*100}")
    print("GENERATING INDIVIDUAL SAMPLE PLOTS")
    print(f"{'='*100}")

    # Create individual plots directory
    os.makedirs(INDIVIDUAL_PLOTS_DIR, exist_ok=True)
    print(f"Individual plots directory: {INDIVIDUAL_PLOTS_DIR}\n")

    # generate_individual_sample_plots returns (pfi_result, cfi_result) tuple
    # so we can reuse them for the cross-sample heatmap without recomputing.
    # CFI is expensive (trains per-feature RF conditional samplers), so
    # caching here avoids doing it twice per sample.
    per_sample_pfi_results = {}
    per_sample_cfi_results = {}

    for sample_name, features in sample_features_dict.items():
        if len(features) > 0:
            try:
                pfi_result, cfi_result = generate_individual_sample_plots(
                    features, sample_name, INDIVIDUAL_PLOTS_DIR)
                if pfi_result is not None:
                    per_sample_pfi_results[sample_name] = pfi_result
                if cfi_result is not None:
                    per_sample_cfi_results[sample_name] = cfi_result
            except Exception as e:
                print(f"ERROR generating individual plots for {sample_name}: {str(e)}")
                import traceback
                traceback.print_exc()
        else:
            print(f"\nSkipping {sample_name} (no features extracted)")

    # ========================================================================
    # CROSS-SAMPLE CFI vs PFI HEATMAP
    # ========================================================================
    # Uses the per-sample PFI and CFI results already computed above during
    # individual sample plotting. The side-by-side heatmap shows how
    # conditional vs marginal importance varies across species.
    print(f"\n{'='*100}")
    print("GENERATING CROSS-SAMPLE CFI vs PFI HEATMAP")
    print(f"{'='*100}")

    if FIPPY_AVAILABLE and len(per_sample_pfi_results) >= 2 and len(per_sample_cfi_results) >= 2:
        heatmap_path = os.path.join(
            OUTPUT_DIR,
            f"murray_drivers_6_cfi_pfi_heatmap.{FIGURE_FORMAT}")
        try:
            plot_cfi_pfi_heatmap(
                per_sample_pfi_results,
                per_sample_cfi_results,
                heatmap_path,
                metric=IMPORTANCE_METRIC
            )
        except Exception as e:
            print(f"ERROR generating CFI/PFI heatmap: {str(e)}")
            import traceback
            traceback.print_exc()
    elif not FIPPY_AVAILABLE:
        print("  fippy not available - skipping cross-sample heatmap")
    else:
        print("  Need >=2 samples with sufficient data for heatmap, skipping")

    # ========================================================================
    # FEATURE CORRELATION TABLES
    # ========================================================================
    print(f"\n{'='*100}")
    print("FEATURE CORRELATION TABLES (Per Sample)")
    print(f"{'='*100}")

    feature_keys = ['alpha', 'r_parent', 'asymmetry', 'branch_angle', 'tortuosity', 'radius_taper']
    feature_labels = ['alpha', 'r_parent', 'asymmetry', 'angle', 'tortuosity', 'taper']

    for sample_name, features in sample_features_dict.items():
        if len(features) < 10:
            print(f"\n{sample_name}: Insufficient data (n={len(features)})")
            continue

        # Extract feature arrays
        feature_arrays = {}
        for key in feature_keys:
            feature_arrays[key] = np.array([f.get(key, np.nan) for f in features])

        # Create mask for valid data (all features finite)
        valid_mask = np.ones(len(features), dtype=bool)
        for key in feature_keys:
            valid_mask &= np.isfinite(feature_arrays[key])

        n_valid = np.sum(valid_mask)
        if n_valid < 10:
            print(f"\n{sample_name}: Insufficient valid data (n={n_valid})")
            continue

        # Filter to valid data
        filtered_arrays = {key: feature_arrays[key][valid_mask] for key in feature_keys}

        # Compute correlation matrix
        n_features = len(feature_keys)
        corr_matrix = np.zeros((n_features, n_features))

        for i, key_i in enumerate(feature_keys):
            for j, key_j in enumerate(feature_keys):
                corr_matrix[i, j] = np.corrcoef(filtered_arrays[key_i], filtered_arrays[key_j])[0, 1]

        # Print table header
        print(f"\n{sample_name} (n={n_valid} bifurcations)")
        print("-" * (12 + 10 * n_features))

        # Column headers
        header = f"{'Feature':<12}"
        for label in feature_labels:
            header += f"{label:>10}"
        print(header)
        print("-" * (12 + 10 * n_features))

        # Print rows
        for i, label in enumerate(feature_labels):
            row = f"{label:<12}"
            for j in range(n_features):
                row += f"{corr_matrix[i, j]:>10.3f}"
            print(row)

        print("-" * (12 + 10 * n_features))

    print(f"\n{'='*100}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*100}")
    print(f"\nGenerated combined plots:")
    print(f"  - {plot1_path}")
    print(f"  - {plot2_path}")
    print(f"  - {plot3_path}")
    print(f"  - {summary_path}")
    if FIPPY_AVAILABLE:
        print(f"  - {pooled_cfi_path}  (CFI vs PFI bar chart, pooled)")
        if len(per_sample_cfi_results) >= 2:
            print(f"  - {os.path.join(OUTPUT_DIR, f'murray_drivers_5_cfi_pfi_heatmap.{FIGURE_FORMAT}')}  (CFI vs PFI heatmap across samples)")
    print(f"\nGenerated individual sample plots in:")
    print(f"  - {INDIVIDUAL_PLOTS_DIR}")
    print(f"    (up to 5 plots per sample: alpha_vs_radius, alpha_vs_asymmetry,")
    print(f"     permutation_importance, summary, cfi_vs_pfi)")
    print()


if __name__ == "__main__":
    main()
