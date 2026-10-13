"""
SHALSTAB - Slope Stability Analysis

Simple Python package for slope stability analysis.

Usage
-----
    import shalstab

    analyzer = shalstab.Analyzer("dem.tif", "geology.geojson")
    results = analyzer.calculate_critical_rainfall()

"""

import importlib.util
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .slope_analyzer import Analyzer

__version__ = "1.1.0.dev0"
__all__ = ["Analyzer", "__version__"]

try:
    if importlib.util.find_spec("shalstab.training_data") is not None:
        __all__.extend(["training_dem", "training_geology"])
except Exception:
    pass


def __getattr__(name: str) -> Any:
    """Lazy import heavy modules only when accessed to keep startup fast and lean."""
    if name == "Analyzer":
        from .slope_analyzer import Analyzer

        return Analyzer
    elif name in ("training_dem", "training_geology"):
        try:
            from . import training_data

            return getattr(training_data, name)
        except ImportError:
            raise AttributeError(
                f"module {__name__!r} has no attribute {name!r}"
            ) from None
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
