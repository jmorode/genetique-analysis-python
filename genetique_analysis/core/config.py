"""
Configuration module for genetics analysis.

This module handles project configuration, file paths, data loading, and project setup.
It provides the FileConfiguration class that manages all aspects of project configuration.
"""

import json
import os
import sys
from typing import Any, Dict, List, Optional

import pandas as pd
from loguru import logger

from ..utils.utils import (
    conversion_two_lines_to_one_lines_genotypes,
    create_folder_if_necessary,
)
from .constants import (
    AGGREGATION_TYPES,
    ANALYSIS_STEPS,
    CSV_SEPARATOR,
    DEFAULT_FIGURE_SIZE,
    DEFAULT_NB_FAMILIES,
    DEFAULT_NB_INDIVIDUALS,
    FILE_ENCODING,
    LARGE_FIGURE_SIZE,
    PLOT_DPI,
)


class FileConfiguration:
    """
    Configuration class for genetics analysis projects.

    This class manages all configuration aspects including file paths, data loading,
    validation, and project setup for genetics analysis workflows.

    Args:
        project_name: Name of the project directory
        agg_type: Aggregation type ("all", "pops", "pop_years", "subcat")
        selection_name: Name of the selection file to use
    """

    def __init__(self, project_name: str, agg_type: str, selection_name: str):
        """
        Initialize FileConfiguration with project parameters.

        Args:
            project_name: Name of the project directory
            agg_type: Aggregation type for analysis
            selection_name: Name of the selection file
        """
        if agg_type not in AGGREGATION_TYPES:
            raise ValueError(
                f"Invalid aggregation type: {agg_type}. Must be one of {AGGREGATION_TYPES}"
            )

        self.project_name = project_name
        self.agg_type = agg_type
        self.selection_name = selection_name

        # Configuration data
        self.config: Dict[str, Any] = {}

        # File paths
        self.input_path = f"./projects/{self.project_name}/inputs"
        self.path_to_config_file = f"{self.input_path}/config.json"
        self.output_path = f"./projects/{self.project_name}/outputs"

        # Data containers
        self.genotypes_two_lines: Optional[pd.DataFrame] = None
        self.genotypes_data_pairwise: Optional[pd.DataFrame] = None
        self.pops_to_select: Optional[pd.DataFrame] = None
        self.dict_pop_samples: Optional[Dict[str, str]] = None
        self.loci_list: Optional[pd.Index] = None

        # Configuration flags
        self.bool_config = True

        # Analysis configuration (will be loaded from config file)
        self.population_order: List[str] = []
        self.analysis_steps: Dict[str, bool] = {}
        self.simulation_parameters: Dict[str, int] = {}
        self.plotting_parameters: Dict[str, Any] = {}
        self.file_parameters: Dict[str, str] = {}

    def _no_config(self) -> None:
        """Set configuration flag to False when config file doesn't exist."""
        self.bool_config = False
        logger.info("Config file does not exist!")

    def _check_config_file(self) -> None:
        """Check if configuration file exists."""
        if not os.path.exists(self.path_to_config_file):
            self._no_config()

    def _retrieve_config(self) -> None:
        """Load configuration from JSON file."""
        try:
            with open(self.path_to_config_file, encoding=FILE_ENCODING) as f:
                data = json.load(f)

            # Load basic configuration
            for key in data.keys():
                self.config[key] = data[key]

            # Load analysis configuration with defaults
            self._load_analysis_config(data)

        except (json.JSONDecodeError, IOError) as e:
            logger.error(f"Error loading config file: {e}")
            self._no_config()

    def _load_analysis_config(self, data: Dict[str, Any]) -> None:
        """Load analysis configuration parameters with defaults."""
        # Population order
        self.population_order = data.get("population_order", [])

        # Analysis steps with defaults (all enabled by default)
        default_steps = {step: True for step in ANALYSIS_STEPS}
        self.analysis_steps = data.get("analysis_steps", default_steps)

        # Simulation parameters with defaults
        self.simulation_parameters = {
            "nb_families": data.get("simulation_parameters", {}).get(
                "nb_families", DEFAULT_NB_FAMILIES
            ),
            "nb_individuals": data.get("simulation_parameters", {}).get(
                "nb_individuals", DEFAULT_NB_INDIVIDUALS
            ),
        }

        # Plotting parameters with defaults
        self.plotting_parameters = {
            "figure_size": data.get("plotting_parameters", {}).get(
                "figure_size", list(DEFAULT_FIGURE_SIZE)
            ),
            "large_figure_size": data.get("plotting_parameters", {}).get(
                "large_figure_size", list(LARGE_FIGURE_SIZE)
            ),
            "dpi": data.get("plotting_parameters", {}).get("dpi", PLOT_DPI),
        }

        # File parameters with defaults
        self.file_parameters = {
            "csv_separator": data.get("file_parameters", {}).get(
                "csv_separator", CSV_SEPARATOR
            ),
            "file_encoding": data.get("file_parameters", {}).get(
                "file_encoding", FILE_ENCODING
            ),
        }

        logger.info("Analysis configuration loaded successfully")

    def should_execute_step(self, step_name: str) -> bool:
        """
        Check if a specific analysis step should be executed.

        Args:
            step_name: Name of the analysis step

        Returns:
            True if the step should be executed, False otherwise
        """
        if not self.bool_config:
            # If no config file, execute all steps by default
            return True
        return self.analysis_steps.get(step_name, True)

    def get_simulation_parameter(self, param_name: str) -> int:
        """
        Get a simulation parameter value.

        Args:
            param_name: Name of the parameter

        Returns:
            Parameter value
        """
        return self.simulation_parameters.get(param_name, DEFAULT_NB_FAMILIES)

    def get_plotting_parameter(self, param_name: str) -> Any:
        """
        Get a plotting parameter value.

        Args:
            param_name: Name of the parameter

        Returns:
            Parameter value
        """
        return self.plotting_parameters.get(param_name, DEFAULT_FIGURE_SIZE)

    def get_file_parameter(self, param_name: str) -> str:
        """
        Get a file parameter value.

        Args:
            param_name: Name of the parameter

        Returns:
            Parameter value
        """
        return self.file_parameters.get(param_name, CSV_SEPARATOR)

    def initialize_folders(self) -> None:
        """Create necessary output directories for the analysis."""
        output_subdirs = [
            "raw_data",
            "heterozygosity",
            "pairwise_differences",
            "pairwise_differences/plots/",
            "recaptures",
            "ml_relate",
            "test_stats",
            "error_recapture_tests",
        ]

        for subdir in output_subdirs:
            create_folder_if_necessary(f"{self.output_path}/{subdir}/")

    def check_population_year(self, selection: pd.DataFrame) -> None:
        """
        Validate that all population-year combinations in selection exist in genotypes.

        Args:
            selection: DataFrame containing population and year selections

        Raises:
            AssertionError: If any population-year combination doesn't exist in genotypes
        """
        if self.genotypes_two_lines is None:
            raise ValueError(
                "Genotypes data not loaded. Call retrieve_and_prepare_inputs first."
            )

        pop_year_geno = self.genotypes_two_lines[
            ["Population", "Year"]
        ].drop_duplicates()

        for pop, year in zip(selection.Population, selection.Year):
            exists = (
                len(
                    pop_year_geno[
                        (pop_year_geno.Population == pop) & (pop_year_geno.Year == year)
                    ]
                )
                != 0
            )

            if not exists:
                raise ValueError(
                    f"Population '{pop}' / Year '{year}' in selection.csv doesn't exist in genotypes.csv"
                )

    def check_population_year_subcat(self, selection: pd.DataFrame) -> None:
        """
        Validate that all population-year-subcategory combinations in selection exist in genotypes.

        Args:
            selection: DataFrame containing population, year, and subcategory selections

        Raises:
            ValueError: If any population-year-subcategory combination doesn't exist in genotypes
        """
        if self.genotypes_two_lines is None:
            raise ValueError(
                "Genotypes data not loaded. Call retrieve_and_prepare_inputs first."
            )

        pop_year_geno = self.genotypes_two_lines[
            ["Population", "Year", "Subcategory"]
        ].drop_duplicates()

        for pop, year, sub in zip(
            selection.Population, selection.Year, selection.Subcategory
        ):
            if pd.notna(sub):
                exists = (
                    len(
                        pop_year_geno[
                            (pop_year_geno.Population == pop)
                            & (pop_year_geno.Year == year)
                            & (pop_year_geno.Subcategory == sub)
                        ]
                    )
                    != 0
                )

                if not exists:
                    raise ValueError(
                        f"Population '{pop}' / Year '{year}' / Subcategory '{sub}' in selection.csv "
                        f"doesn't exist in genotypes.csv"
                    )

    def check_loci_names(self, conversion: pd.DataFrame) -> None:
        """
        Validate that loci names in conversion file match those in genotypes file.

        Args:
            conversion: DataFrame containing locus conversion information

        Raises:
            ValueError: If loci names don't match between conversion and genotypes files
        """
        if self.genotypes_two_lines is None:
            raise ValueError(
                "Genotypes data not loaded. Call retrieve_and_prepare_inputs first."
            )

        metadata_cols = ["Population", "Sample", "Year", "Subcategory"]
        loci_genotypes = set(self.genotypes_two_lines.columns.drop(metadata_cols))
        loci_conversion = set(conversion["Loci"].unique())

        if loci_conversion != loci_genotypes:
            raise ValueError(
                f"Loci mismatch: conversion.csv has {loci_conversion} "
                f"but genotypes.csv has {loci_genotypes}"
            )

    def convert_raw_genotypes(self, conversion_table: pd.DataFrame) -> None:
        """
        Convert raw allele values to corrected values using conversion table.

        Args:
            conversion_table: DataFrame with columns 'Loci', 'allele_raw', 'allele_corrected'
        """
        if self.genotypes_two_lines is None:
            raise ValueError(
                "Genotypes data not loaded. Call retrieve_and_prepare_inputs first."
            )

        for loci in conversion_table.Loci.unique():
            subconv_table = conversion_table[conversion_table.Loci == loci].copy()
            subconv_table = subconv_table.drop("Loci", axis=1)
            subconv_dict = subconv_table.set_index("allele_raw").to_dict("dict")

            self.genotypes_two_lines[loci] = self.genotypes_two_lines[loci].apply(
                lambda x: (
                    int(subconv_dict["allele_corrected"][x]) if pd.notnull(x) else None
                )
            )

        # Save converted genotypes
        output_file = f"{self.output_path}/raw_data/genotypes.csv"
        self.genotypes_two_lines.to_csv(output_file, sep=CSV_SEPARATOR, index=False)
        logger.info(f"Converted genotypes saved to {output_file}")

    def convert_genotypes_two_lines_to_one_line(self) -> None:
        """
        Convert two-line genotype format to one-line format for pairwise analysis.
        """
        if self.genotypes_two_lines is None:
            raise ValueError(
                "Genotypes data not loaded. Call retrieve_and_prepare_inputs first."
            )

        input_data = self.genotypes_two_lines.copy()
        loci_cols = input_data.columns
        meta_cols = ["Sample", "Population", "Year", "Subcategory"]
        loci_cols = loci_cols.drop(meta_cols)

        data = conversion_two_lines_to_one_lines_genotypes(
            input_data, meta_cols, loci_cols
        )

        # Save pairwise format
        output_file = f"{self.output_path}/raw_data/genotypes_pairwise_diff.csv"
        data.to_csv(output_file, sep=CSV_SEPARATOR, index=False, header=False)
        logger.info(f"Pairwise format genotypes saved to {output_file}")

        self.genotypes_data_pairwise = data.copy()

    def retrieve_and_prepare_inputs(self, selection_name: str) -> None:
        """
        Load and prepare all input data for analysis.

        Args:
            selection_name: Name of the selection file to load
        """
        logger.info(f"Loading data from {self.input_path}")

        # Load genotypes file
        genotypes_file = f"{self.input_path}/genotypes.csv"
        if not os.path.exists(genotypes_file):
            logger.error(f"Genotypes file not found: {genotypes_file}")
            sys.exit(1)

        self.genotypes_two_lines = pd.read_csv(genotypes_file, sep=CSV_SEPARATOR)
        logger.info(f"Loaded genotypes: {self.genotypes_two_lines.shape}")

        # Load and apply conversion if available
        conversion_file = f"{self.input_path}/conversion.csv"
        if os.path.exists(conversion_file):
            conversion = pd.read_csv(conversion_file, sep=CSV_SEPARATOR)
            self.check_loci_names(conversion)
            self.convert_raw_genotypes(conversion)
        else:
            logger.info("No conversion file found, skipping allele conversion")

        # Convert to pairwise format
        self.convert_genotypes_two_lines_to_one_line()

        # Save all population/year/subcategory combinations
        self.save_all_pop_year_subcategory_combinations()

        # Load selection data
        selection_file = f"{self.input_path}/selection_{selection_name}.csv"
        if not os.path.exists(selection_file):
            logger.error(f"Selection file not found: {selection_file}")
            sys.exit(1)

        self.pops_to_select = pd.read_csv(selection_file, sep=CSV_SEPARATOR)
        logger.info(f"Loaded selection: {self.pops_to_select.shape}")

        # Validate selection data
        self.check_population_year(self.pops_to_select)
        if self.agg_type == "subcat":
            self.check_population_year_subcat(self.pops_to_select)

    def initialize_global(self) -> None:
        """Initialize the complete configuration and load all data."""
        self._check_config_file()
        if self.bool_config:
            self._retrieve_config()
        self.initialize_folders()
        self.retrieve_and_prepare_inputs(self.selection_name)
        self.init_dict_pop_samples()
        self.initialize_loci_list()

    def save_all_pop_year_subcategory_combinations(self) -> None:
        """Save all unique population/year/subcategory combinations to file."""
        if self.genotypes_two_lines is None:
            raise ValueError("Genotypes data not loaded.")

        combinations = self.genotypes_two_lines[
            ["Population", "Year", "Subcategory"]
        ].drop_duplicates()

        output_file = f"{self.input_path}/select_everything.csv"
        combinations.to_csv(output_file, sep=CSV_SEPARATOR, index=False)
        logger.info(f"Saved all combinations to {output_file}")

    def init_dict_pop_samples(self) -> None:
        """Initialize dictionary mapping samples to population legends based on aggregation type."""
        if self.genotypes_two_lines is None:
            raise ValueError("Genotypes data not loaded.")

        pop_samples = self.genotypes_two_lines.copy()

        if self.agg_type == "pops":
            pop_samples["legend"] = pop_samples["Population"]
        elif self.agg_type in ["all", "pop_years"]:
            pop_samples["legend"] = (
                pop_samples["Population"] + " - " + pop_samples["Year"].astype(str)
            )
        else:  # subcat
            pop_samples["legend"] = (
                pop_samples["Population"]
                + " - "
                + pop_samples["Year"].astype(str)
                + " - "
                + pop_samples["Subcategory"].astype(str)
            )

        pop_samples = pop_samples[["Sample", "legend"]].drop_duplicates()
        self.dict_pop_samples = dict(zip(pop_samples["Sample"], pop_samples["legend"]))
        logger.info(
            f"Initialized population sample mapping for {len(self.dict_pop_samples)} samples"
        )

    def initialize_loci_list(self) -> None:
        """Initialize list of loci from genotype data."""
        if self.genotypes_two_lines is None:
            raise ValueError("Genotypes data not loaded.")

        metadata_cols = ["Population", "Sample", "Year", "Subcategory"]
        self.loci_list = self.genotypes_two_lines.columns.drop(metadata_cols)
        logger.info(f"Initialized loci list with {len(self.loci_list)} loci")
