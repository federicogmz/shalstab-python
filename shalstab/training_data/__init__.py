"""
SHALSTAB Training Data

This module provides easy access to training data for the SHALSTAB package.
The training data includes sample DEM and geology files for learning and testing.

Files available:
- dem.tif: Digital Elevation Model
- geology.geojson: Geological units with geotechnical parameters
"""

from pathlib import Path

# Define paths to training data files
_data_dir = Path(__file__).parent
training_dem = str(_data_dir / "dem.tif")
training_geology = str(_data_dir / "geology.geojson")

__all__ = ["training_dem", "training_geology"]
