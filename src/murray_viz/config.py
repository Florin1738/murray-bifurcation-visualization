#!/usr/bin/env python3
"""
Centralized configuration for Murray Bifurcation Visualization.

This module provides default file paths and configurations for the visualization.
Update these paths to match your data location.
"""

from pathlib import Path

# Data directory (relative to project root)
DATA_DIR = Path(__file__).parent.parent.parent / "data"

# Default dataset configurations
DATASETS = {
    'chicken_liver': {
        'skeleton_path': str(DATA_DIR / "chicken_liver_skeleton.pkl"),
        'volume_path': None,  # Skeleton-only mode
        'volume_var': 'rescaled_vol',
        'volume_spacing': (4, 4, 8)
    },

    'GS40_liver': {
        'skeleton_path': str(DATA_DIR / "GS40_liver_skeleton.pkl"),
        'volume_path': r"\\10.162.80.16\Andre_expansion\students\Marco\data\FINAL_VOLS\GS40S1\whole_species\liver_blood_vessels\downsampled\GS40_liver_vesels_4xy.mat",
        'volume_var': 'rescaled_vol',
        'volume_spacing': (16, 16, 8)
    },

    'GS55_brain': {
        'skeleton_path': str(DATA_DIR / "GS55_skeleton2.pkl"),
        'volume_path': r"\\10.162.80.11\Andre_kit\data\monkey_fetus\bissected_monkey_GS55\10x_python\registeredE\int_2\5x_downsampled_images\classification_MODEL1_5x_GS40_GS55_06_10_2025_45_big_tiles_inceptionresnetv2\volume\Processed_Volumes\02_brain_smooth_4xy.mat",
        'volume_var': 'rescaled_vol',
        'volume_spacing': (16, 16, 16)
    },

    'GS55_skeleton_only': {
        'skeleton_path': str(DATA_DIR / "GS55_skeleton2.pkl"),
        'volume_path': None,  # Skeleton-only mode
        'volume_var': 'rescaled_vol',
        'volume_spacing': (16, 16, 16)
    }
}

# Default dataset to use
DEFAULT_DATASET = 'chicken_liver'


def get_dataset_config(dataset_name=None):
    """
    Get configuration for a specific dataset.

    Parameters
    ----------
    dataset_name : str, optional
        Name of the dataset. If None, returns DEFAULT_DATASET config.

    Returns
    -------
    dict
        Dataset configuration with keys: skeleton_path, volume_path, volume_var, volume_spacing

    Raises
    ------
    ValueError
        If dataset_name is not found in DATASETS
    """
    if dataset_name is None:
        dataset_name = DEFAULT_DATASET

    if dataset_name not in DATASETS:
        available = ', '.join(DATASETS.keys())
        raise ValueError(f"Dataset '{dataset_name}' not found. Available datasets: {available}")

    return DATASETS[dataset_name].copy()


def list_datasets():
    """
    List available dataset names.

    Returns
    -------
    list
        List of available dataset names
    """
    return list(DATASETS.keys())
