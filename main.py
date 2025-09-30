"""
Main entry point for genetics analysis pipeline.

This module orchestrates the complete genetics analysis workflow including:
- Frequency and heterozygosity calculations
- Pairwise difference analysis
- Recapture analysis
- ML-related analysis
- Statistical testing
"""

import json
import os
import sys
import warnings

# Local imports
from config import FileConfiguration
from description import compute_frequency_and_heterozygosity, save_proba_same_ind
from generation import (
    plot_cumulated_and_not_frequencies_for_simu_relationships_and_sample_by_pop,
    plot_pairwise_distances_per_relationships,
)
from loguru import logger
from ml_relate import (
    count_nbr_relationships_from_ml_relate_output,
    reliability_ml_relate_based_simulated_families,
    write_ml_relate_input_file_genotypes,
    write_ml_relate_input_file_simu_families,
)
from pairwise_differences import (
    plot_pairwise_distances_all_intra_populations,
    plot_pairwise_distances_one_by_one,
    plot_pairwise_distances_two_by_two,
    recover_pairwise_difference,
)
from recaptures import get_recaptures
from test_stats import compute_ks_test_between_pairwise_distance_and_relationships_simu

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")


def load_configuration(project_name: str) -> FileConfiguration:
    """
    Load configuration from JSON file and create FileConfiguration object.

    Args:
        project_name: Name of the project directory

    Returns:
        FileConfiguration object with loaded parameters

    Raises:
        FileNotFoundError: If configuration file doesn't exist
        ValueError: If configuration file has invalid parameters
    """
    config_path = f"./projects/{project_name}/inputs/config.json"

    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"Configuration file not found: {config_path}\n"
            f"Please create a config.json file in projects/{project_name}/inputs/"
        )

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            config_data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in configuration file: {e}")
    except Exception as e:
        raise ValueError(f"Error reading configuration file: {e}")

    # Extract required parameters
    selection_name = config_data.get("selection_name")
    aggregation_type = config_data.get("aggregation_type")

    if not selection_name:
        raise ValueError(
            "Missing required parameter 'selection_name' in configuration file"
        )
    if not aggregation_type:
        raise ValueError(
            "Missing required parameter 'aggregation_type' in configuration file"
        )

    # Create and initialize configuration
    config = FileConfiguration(project_name, aggregation_type, selection_name)
    config.initialize_global()

    return config


def main(project_name: str) -> None:
    """
    Main analysis pipeline for genetics data.

    This function orchestrates the complete genetics analysis workflow:
    1. Configuration and data loading
    2. Frequency and heterozygosity analysis
    3. Pairwise difference calculations and plotting
    4. Recapture analysis
    5. & 6. ML-related analysis and reliability testing
    7. Statistical testing and relationship analysis

    Args:
        project_name: Name of the project directory
    """
    logger.info("Starting genetics analysis pipeline")

    # Load configuration from file
    config = load_configuration(project_name)
    logger.info("Configuration and data loaded successfully")

    # Get parameters from configuration
    population_order = config.population_order
    nb_families = config.get_simulation_parameter("nb_families")

    # Step 1: Frequency and heterozygosity analysis
    if config.should_execute_step("step_1_frequency_heterozygosity"):
        logger.info("Step 1: Computing frequencies and heterozygosity")
        save_proba_same_ind(config)
        compute_frequency_and_heterozygosity(config, population_order)
    else:
        logger.info("Step 1: Skipped (disabled in configuration)")

    # Step 2: Pairwise difference analysis
    if config.should_execute_step("step_2_pairwise_differences"):
        logger.info("Step 2: Computing pairwise differences")
        recover_pairwise_difference(config)
    else:
        logger.info("Step 2: Skipped (disabled in configuration)")

    # Step 3: Plotting pairwise distances
    if config.should_execute_step("step_3_plot_pairwise_distances"):
        logger.info("Step 3: Plotting pairwise distances")
        plot_pairwise_distances_one_by_one(config)
        if config.agg_type != "all":
            plot_pairwise_distances_two_by_two(config)
        if config.agg_type == "pops":
            plot_pairwise_distances_all_intra_populations(config)
    else:
        logger.info("Step 3: Skipped (disabled in configuration)")

    # Step 4: Recapture analysis
    if config.should_execute_step("step_4_recaptures"):
        logger.info("Step 4: Analyzing recaptures")
        get_recaptures(config)
    else:
        logger.info("Step 4: Skipped (disabled in configuration)")

    # Step 5: ML-related analysis INPUT FILE
    if config.should_execute_step("step_5_ml_relate_input"):
        logger.info("Step 5: ML-related analysis INPUT FILE")
        write_ml_relate_input_file_genotypes(config)
        write_ml_relate_input_file_simu_families(config, nb_families=nb_families)
    else:
        logger.info("Step 5: Skipped (disabled in configuration)")

    # Step 6: ML-related analysis OUTPUT FILE
    if config.should_execute_step("step_6_ml_relate_output"):
        logger.info("Step 6: ML-related analysis OUTPUT FILE")
        # Note: ML relate output file should be added to input folder as ml_relate_output_{selection_name}.csv
        count_nbr_relationships_from_ml_relate_output(config)
        reliability_ml_relate_based_simulated_families(config, nb_families=nb_families)
    else:
        logger.info("Step 6: Skipped (disabled in configuration)")

    # Step 7: Relationship analysis and statistical testing
    if config.should_execute_step("step_7_relationship_analysis"):
        logger.info("Step 7: Relationship analysis and statistical testing")
        plot_pairwise_distances_per_relationships(config, nb_families=nb_families)
        compute_ks_test_between_pairwise_distance_and_relationships_simu(
            config, nb_families=nb_families
        )

        if config.agg_type == "pops":
            plot_cumulated_and_not_frequencies_for_simu_relationships_and_sample_by_pop(
                config
            )
    else:
        logger.info("Step 7: Skipped (disabled in configuration)")

    logger.info("Genetics analysis pipeline completed successfully")


if __name__ == "__main__":
    # Get project name from command line argument or use default
    if len(sys.argv) > 1:
        project_name = sys.argv[1]
    else:
        project_name = "analysis_jan_2025"  # Default project name

    logger.info(f"Starting analysis for project: {project_name}")

    try:
        # Run the analysis pipeline
        # All parameters are loaded from projects/{project_name}/inputs/config.json
        main(project_name)

    except FileNotFoundError as e:
        logger.error(f"Configuration file not found: {e}")
        logger.info(
            "Please create a config.json file in your project's inputs directory."
        )
        logger.info("See CONFIG_README.md for configuration file format and examples.")
        sys.exit(1)

    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        logger.info("Please check your config.json file format.")
        logger.info("See CONFIG_README.md for configuration file format and examples.")
        sys.exit(1)

    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)
