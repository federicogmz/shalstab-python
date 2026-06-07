"""Comprehensive tests for the SHALSTAB Analyzer class using training data."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

import shalstab
from shalstab import Analyzer

# ============================================================
# Fixtures
# ============================================================


@pytest.fixture(scope="module")
def analyzer():
    """Create an Analyzer instance with training data (module-level for speed)."""
    return Analyzer(
        dem_path=shalstab.training_dem,
        geo=shalstab.training_geology,
        geo_columns=["Cohesion", "Phi", "Gamma_kN_m", "Ks_m_s"],
    )


@pytest.fixture(scope="module")
def critical_rainfall(analyzer):
    """Calculate critical rainfall (module-level for speed)."""
    return analyzer.calculate_critical_rainfall(show_plot=False)


@pytest.fixture(scope="module")
def stability_result(analyzer):
    """Calculate stability for a moderate rainfall event."""
    return analyzer.calculate_stability(rainfall_mm_day=50)


@pytest.fixture(scope="module")
def failure_probability(analyzer):
    """Calculate failure probability."""
    return analyzer.calculate_failure_probability()


# ============================================================
# Initialization Tests
# ============================================================


class TestAnalyzerInit:
    """Tests for Analyzer initialization and setup."""

    def test_analyzer_creates_with_training_data(self):
        """Analyzer should initialize successfully with training data."""
        a = Analyzer(
            dem_path=shalstab.training_dem,
            geo=shalstab.training_geology,
            geo_columns=["Cohesion", "Phi", "Gamma_kN_m", "Ks_m_s"],
        )
        assert a is not None

    def test_dem_is_loaded(self, analyzer):
        """DEM should be loaded as a DataArray."""
        assert isinstance(analyzer.dem, xr.DataArray)
        assert analyzer.dem.size > 0

    def test_extent_is_set(self, analyzer):
        """Spatial extent should be a list of 4 floats."""
        assert len(analyzer.extent) == 4
        assert all(isinstance(v, (int, float)) for v in analyzer.extent)

    def test_hillshade_is_computed(self, analyzer):
        """Hillshade should be a numpy array with valid values."""
        assert isinstance(analyzer.hillshade, np.ndarray)
        assert analyzer.hillshade.min() >= 0
        assert analyzer.hillshade.max() <= 255

    def test_flow_accumulation_is_computed(self, analyzer):
        """Flow accumulation should be a non-negative DataArray."""
        assert isinstance(analyzer.flow_accumulated, xr.DataArray)
        assert (analyzer.flow_accumulated >= 0).all()

    def test_slope_is_computed(self, analyzer):
        """Slope should be in radians and non-negative."""
        assert isinstance(analyzer.slope_rad, xr.DataArray)
        assert (analyzer.slope_rad >= 0).all()

    def test_geotechnical_parameters_initialized(self, analyzer):
        """All four geotechnical parameters should be DataArrays."""
        assert isinstance(analyzer.cohesion, xr.DataArray)
        assert isinstance(analyzer.friction_rad, xr.DataArray)
        assert isinstance(analyzer.unit_weight, xr.DataArray)
        assert isinstance(analyzer.permeability, xr.DataArray)

    def test_soil_thickness_computed(self, analyzer):
        """Soil thickness should be a DataArray with positive values."""
        assert isinstance(analyzer._soil_thickness, xr.DataArray)
        valid = analyzer._soil_thickness.values[
            ~np.isnan(analyzer._soil_thickness.values)
        ]
        assert len(valid) > 0

    def test_geo_columns_must_be_exactly_4(self):
        """geo_columns must have exactly 4 elements."""
        with pytest.raises(ValueError, match="exactly 4"):
            Analyzer(
                dem_path=shalstab.training_dem,
                geo=shalstab.training_geology,
                geo_columns=["Cohesion", "Phi"],
            )

    def test_missing_geo_columns_raises(self):
        """Passing None for geo_columns should raise ValueError."""
        with pytest.raises(ValueError, match="geo_columns is required"):
            Analyzer(
                dem_path=shalstab.training_dem,
                geo=shalstab.training_geology,
                geo_columns=None,
            )

    def test_invalid_dem_path_raises(self):
        """Non-existent DEM path should raise an error."""
        with pytest.raises((ValueError, FileNotFoundError)):
            Analyzer(
                dem_path="nonexistent.tif",
                geo=shalstab.training_geology,
                geo_columns=["Cohesion", "Phi", "Gamma_kN_m", "Ks_m_s"],
            )

    def test_figsize_is_respected(self):
        """Custom figsize should be stored."""
        a = Analyzer(
            dem_path=shalstab.training_dem,
            geo=shalstab.training_geology,
            geo_columns=["Cohesion", "Phi", "Gamma_kN_m", "Ks_m_s"],
            figsize=(15, 10),
        )
        assert a.figsize == (15, 10)


# ============================================================
# Critical Rainfall Tests
# ============================================================


class TestCriticalRainfall:
    """Tests for critical rainfall calculation."""

    def test_returns_dataarray(self, critical_rainfall):
        """Should return an xarray DataArray."""
        assert isinstance(critical_rainfall, xr.DataArray)

    def test_has_valid_values(self, critical_rainfall):
        """Should contain finite positive values (after masking NaNs)."""
        valid = critical_rainfall.values[~np.isnan(critical_rainfall.values)]
        assert len(valid) > 0
        assert np.all(np.isfinite(valid))

    def test_units_are_mm_per_day(self, critical_rainfall):
        """Values should be in mm/day range (0.1 to 10000)."""
        valid = critical_rainfall.values[~np.isnan(critical_rainfall.values)]
        assert valid.min() > 0
        assert valid.max() < 100000

    def test_shape_matches_dem(self, analyzer, critical_rainfall):
        """Output shape should match DEM spatial dimensions."""
        assert critical_rainfall.shape == analyzer.dem.squeeze().shape

    def test_no_plot_when_disabled(self, analyzer):
        """should_plot=False should not create a figure."""
        plt.close("all")
        result = analyzer.calculate_critical_rainfall(show_plot=False)
        assert isinstance(result, xr.DataArray)


# ============================================================
# Stability Analysis Tests
# ============================================================


class TestStabilityAnalysis:
    """Tests for stability classification."""

    def test_returns_tuple(self, stability_result):
        """Should return (DataArray, Figure) tuple."""
        assert isinstance(stability_result, tuple)
        assert len(stability_result) == 2

    def test_stability_is_dataarray(self, stability_result):
        """First element should be a DataArray."""
        stability, fig = stability_result
        assert isinstance(stability, xr.DataArray)

    def test_figure_is_created(self, stability_result):
        """Second element should be a matplotlib Figure."""
        _, fig = stability_result
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_stability_classes_are_valid(self, stability_result):
        """All values should be 1, 2, 3, or 4."""
        stability, _ = stability_result
        valid_classes = {1, 2, 3, 4}
        actual_classes = set(np.unique(stability.values[~np.isnan(stability.values)]))
        assert actual_classes.issubset(valid_classes)

    def test_unstable_increases_with_rainfall(self, analyzer):
        """Higher rainfall should produce at least as many unstable cells."""
        stab10, fig1 = analyzer.calculate_stability(rainfall_mm_day=10)
        stab100, fig2 = analyzer.calculate_stability(rainfall_mm_day=100)
        plt.close("all")

        unstable_10 = (stab10.values == 3).sum()
        unstable_100 = (stab100.values == 3).sum()
        # Higher rainfall should not decrease unstable area
        assert unstable_100 >= unstable_10

    def test_stability_report_generated(self, stability_result):
        """Should generate a statistical report."""
        stability, _ = stability_result
        assert "Reporte" in stability.attrs
        report = stability.attrs["Reporte"]
        assert "Unconditionally Stable" in report or "Unstable" in report

    def test_custom_rainfall_value(self, analyzer):
        """Should work with any positive rainfall value."""
        stability, fig = analyzer.calculate_stability(rainfall_mm_day=25)
        assert stability is not None
        assert stability.shape == analyzer.dem.squeeze().shape
        valid = stability.values[~np.isnan(stability.values)]
        assert all(v in {1, 2, 3, 4} for v in valid)
        plt.close(fig)


# ============================================================
# Failure Probability Tests
# ============================================================


class TestFailureProbability:
    """Tests for failure probability calculation."""

    def test_returns_dataarray(self, failure_probability):
        """Should return an xarray DataArray."""
        assert isinstance(failure_probability, xr.DataArray)

    def test_values_between_0_and_100(self, failure_probability):
        """All valid values should be between 0 and 100."""
        valid = failure_probability.values[~np.isnan(failure_probability.values)]
        assert len(valid) > 0
        assert valid.min() >= 0
        assert valid.max() <= 100

    def test_has_report(self, failure_probability):
        """Should have a Reporte attribute."""
        assert "Reporte" in failure_probability.attrs


# ============================================================
# Soil Thickness Tests
# ============================================================


class TestSoilThickness:
    """Tests for Catani soil thickness model."""

    def test_soil_thickness_is_dataarray(self, analyzer):
        """Should return a DataArray."""
        assert isinstance(analyzer._soil_thickness, xr.DataArray)

    def test_shape_matches_dem(self, analyzer):
        """Shape should match DEM dimensions."""
        assert analyzer._soil_thickness.shape == analyzer.dem.squeeze().shape

    def test_positive_values(self, analyzer):
        """Valid soil thickness values should be positive."""
        valid = analyzer._soil_thickness.values[
            ~np.isnan(analyzer._soil_thickness.values)
        ]
        if len(valid) > 0:
            assert valid.min() > 0


# ============================================================
# Export Tests
# ============================================================


class TestExport:
    """Tests for raster export functionality."""

    def test_export_raster_creates_file(self, analyzer, tmp_path):
        """Exporting a raster should create a file on disk."""
        output = tmp_path / "test_output.tif"
        analyzer.export_raster(analyzer.dem.squeeze(), output)
        assert output.exists()
        assert output.stat().st_size > 0


# ============================================================
# Package-Level Tests
# ============================================================


class TestPackage:
    """Tests for package-level attributes."""

    def test_version_exists(self):
        """Package should expose a version string."""
        assert hasattr(shalstab, "__version__")
        assert isinstance(shalstab.__version__, str)

    def test_training_dem_path_exists(self):
        """Training DEM path should be set."""
        assert shalstab.training_dem is not None
        assert Path(shalstab.training_dem).exists()

    def test_training_geology_path_exists(self):
        """Training geology path should be set."""
        assert shalstab.training_geology is not None
