"""
Core module for genetique analysis.

This module contains the core configuration and constants used throughout
the genetics analysis pipeline.
"""

# Import modules lazily to avoid dependency issues
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


# Lazy imports to avoid dependency issues at module level
def __getattr__(name):
    if name == "FileConfiguration":
        from .config import FileConfiguration

        return FileConfiguration
    elif name == "AGGREGATION_TYPES":
        from .constants import AGGREGATION_TYPES

        return AGGREGATION_TYPES
    elif name == "ANALYSIS_STEPS":
        from .constants import ANALYSIS_STEPS

        return ANALYSIS_STEPS
    elif name == "CSV_SEPARATOR":
        from .constants import CSV_SEPARATOR

        return CSV_SEPARATOR
    elif name == "DEFAULT_FIGURE_SIZE":
        from .constants import DEFAULT_FIGURE_SIZE

        return DEFAULT_FIGURE_SIZE
    elif name == "DEFAULT_NB_FAMILIES":
        from .constants import DEFAULT_NB_FAMILIES

        return DEFAULT_NB_FAMILIES
    elif name == "DEFAULT_NB_INDIVIDUALS":
        from .constants import DEFAULT_NB_INDIVIDUALS

        return DEFAULT_NB_INDIVIDUALS
    elif name == "FILE_ENCODING":
        from .constants import FILE_ENCODING

        return FILE_ENCODING
    elif name == "LARGE_FIGURE_SIZE":
        from .constants import LARGE_FIGURE_SIZE

        return LARGE_FIGURE_SIZE
    elif name == "PLOT_DPI":
        from .constants import PLOT_DPI

        return PLOT_DPI
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
