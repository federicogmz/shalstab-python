"""Basic smoke tests for shalstab package."""
import pytest


def test_package_imports():
    """Test that the shalstab package can be imported."""
    import shalstab
    assert hasattr(shalstab, "Analyzer")


def test_analyzer_class_exists():
    """Test that the Analyzer class is accessible."""
    from shalstab import Analyzer
    assert callable(Analyzer)


def test_training_data_available():
    """Test that training data paths are provided."""
    import shalstab
    assert shalstab.training_dem is not None
    assert shalstab.training_geology is not None


def test_package_has_version():
    """Test that the package exposes a version."""
    import shalstab
    assert hasattr(shalstab, "__version__")
