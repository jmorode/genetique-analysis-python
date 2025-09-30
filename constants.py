"""
Constants module for genetics analysis.

This module contains all constant values used throughout the genetics analysis pipeline,
including ML-related constants, family relationship mappings, and other fixed values.
"""

from typing import Any, Dict

import numpy as np

# ML Relate Constants
# Reference array for ML-related relationship comparisons
REFERENCE_ARRAY = np.array(
    [
        ["-", np.nan, np.nan, np.nan, np.nan, np.nan],
        ["U", "-", np.nan, np.nan, np.nan, np.nan],
        ["U", "U", "-", np.nan, np.nan, np.nan],
        ["PO", "PO", "U", "-", np.nan, np.nan],
        ["PO", "PO", "U", "FS", "-", np.nan],
        ["U", "PO", "PO", "HS", "HS", "-"],
    ],
    dtype=object,
)

# Family member ordering for ML-related analysis
ORDER_INDS_FAMILY: Dict[str, int] = {
    "P1": 1,
    "M": 2,
    "P2": 3,
    "F1": 4,
    "F2": 5,
    "DF": 6,
}

# Number of individuals per relationship type
NB_IND_RELATION: Dict[str, int] = {"PO": 6, "U": 6, "FS": 1, "HS": 2}

# Generation Constants
# Mapping for renaming parent individuals in family generation
RENAME_PARENTS: Dict[str, str] = {"Ind1": "P1", "Ind2": "M", "Ind3": "P2"}

# Analysis Constants
# Supported aggregation types
AGGREGATION_TYPES = ["all", "pops", "pop_years", "subcat"]

# Analysis steps
ANALYSIS_STEPS = [
    "step_1_frequency_heterozygosity",
    "step_2_pairwise_differences",
    "step_3_plot_pairwise_distances",
    "step_4_recaptures",
    "step_5_ml_relate_input",
    "step_6_ml_relate_output",
    "step_7_relationship_analysis",
]

# File extensions and separators
CSV_SEPARATOR = ";"
FILE_ENCODING = "utf-8"

# Plotting constants
DEFAULT_FIGURE_SIZE = (8, 6)
LARGE_FIGURE_SIZE = (30, 10)
PLOT_DPI = 300

# Simulation constants
DEFAULT_NB_FAMILIES = 1000
DEFAULT_NB_INDIVIDUALS = 100

# Error simulation constants
ERROR_RATES = [0, 0.01, 0.1]
ERROR_TYPES = ["homozygotes", "heterozygotes", "random"]
HETEROZYGOTE_ERROR_TYPES = ["hetero-random", "hetero-biased"]

# Statistical constants
KS_TEST_ALPHA = 0.05
Q95_THRESHOLD = 0.95

# Backward compatibility - keeping original names for existing code
reference_array = REFERENCE_ARRAY
order_inds_family = ORDER_INDS_FAMILY
nb_ind_relation = NB_IND_RELATION
rename_parents = RENAME_PARENTS
