"""
SHALSTAB - Slope Stability Analysis

Simple Python package for slope stability analysis.

Usage
-----
    import shalstab

    analyzer = shalstab.Analyzer("dem.tif", "geology.geojson")
    results = analyzer.calculate_critical_rainfall()

"""

from .slope_analyzer import Analyzer

# Training data
try:
    from .training_data import training_dem, training_geology

    __all__ = ["Analyzer", "training_dem", "training_geology"]
except ImportError:
    __all__ = ["Analyzer"]

__version__ = "1.0.0"
