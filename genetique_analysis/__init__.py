"""
Genetique Analysis - A comprehensive genetics analysis package.

This package provides tools for genetic data analysis including:
- Frequency and heterozygosity calculations
- Pairwise difference analysis
- Recapture analysis
- ML-related analysis and reliability testing
- Statistical testing and relationship analysis

The package is designed to be configuration-driven, allowing users to
control analysis steps and parameters through JSON configuration files.
"""

__version__ = "1.0.0"
__author__ = "Genetics Analysis Team"
__email__ = "genetics@example.com"

# Import core modules lazily to avoid dependency issues
__all__ = [
    "FileConfiguration",
    "AGGREGATION_TYPES",
    "ANALYSIS_STEPS",
    "CSV_SEPARATOR",
    "DEFAULT_FIGURE_SIZE",
    "DEFAULT_NB_FAMILIES",
    "DEFAULT_NB_INDIVIDUALS",
    "FILE_ENCODING",
    "LARGE_FIGURE_SIZE",
    "PLOT_DPI",
]


# Lazy imports to avoid dependency issues at package level
def __getattr__(name):
    if name == "FileConfiguration":
        from .core.config import FileConfiguration

        return FileConfiguration
    elif name == "AGGREGATION_TYPES":
        from .core.constants import AGGREGATION_TYPES

        return AGGREGATION_TYPES
    elif name == "ANALYSIS_STEPS":
        from .core.constants import ANALYSIS_STEPS

        return ANALYSIS_STEPS
    elif name == "CSV_SEPARATOR":
        from .core.constants import CSV_SEPARATOR

        return CSV_SEPARATOR
    elif name == "DEFAULT_FIGURE_SIZE":
        from .core.constants import DEFAULT_FIGURE_SIZE

        return DEFAULT_FIGURE_SIZE
    elif name == "DEFAULT_NB_FAMILIES":
        from .core.constants import DEFAULT_NB_FAMILIES

        return DEFAULT_NB_FAMILIES
    elif name == "DEFAULT_NB_INDIVIDUALS":
        from .core.constants import DEFAULT_NB_INDIVIDUALS

        return DEFAULT_NB_INDIVIDUALS
    elif name == "FILE_ENCODING":
        from .core.constants import FILE_ENCODING

        return FILE_ENCODING
    elif name == "LARGE_FIGURE_SIZE":
        from .core.constants import LARGE_FIGURE_SIZE

        return LARGE_FIGURE_SIZE
    elif name == "PLOT_DPI":
        from .core.constants import PLOT_DPI

        return PLOT_DPI
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
