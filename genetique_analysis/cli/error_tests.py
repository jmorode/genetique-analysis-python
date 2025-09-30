"""
Error testing module for genetique analysis.

This module provides functionality for testing the impact of genotyping errors
on recapture rates and other genetic analysis metrics.
"""

from loguru import logger

from ..analysis.description import save_proba_same_ind
from ..analysis.error_generation import test_error_impact_scenarios_on_given_pop
from ..core.config import FileConfiguration


def main_error_recaptures(
    project_name: str,
    selection_name: str,
    pop_name: str,
    nb_inds_simu: int,
    nb_iterations_simu: int,
) -> None:
    """
    Test the impact of genotyping errors on recapture rates.

    Args:
        project_name: Name of the project directory
        selection_name: Name of the selection file
        pop_name: Name of the population to test
        nb_inds_simu: Number of individuals for simulation
        nb_iterations_simu: Number of iterations for simulation
    """
    # Initialize configuration
    config = FileConfiguration(project_name, "pops", selection_name)
    config.initialize_global()
    logger.info("Configuration and data retrieved.")

    # Test impact of errors of genotyping in recapture rates
    save_proba_same_ind(config)
    test_error_impact_scenarios_on_given_pop(
        config, pop_name, nb_inds_simu, nb_iterations_simu
    )

    logger.info("Error impact testing completed.")


def cli_error_tests():
    """Command-line interface for error testing."""
    import sys

    # Default parameters
    name_project = "analysis_nov_2024_more_loci"
    name_selection = "pachon_sabinos"
    pop_name = "Pachon"
    nb_inds_simu = 100
    nb_iterations_simu = 1000

    # Parse command line arguments if provided
    if len(sys.argv) > 1:
        name_project = sys.argv[1]
    if len(sys.argv) > 2:
        name_selection = sys.argv[2]
    if len(sys.argv) > 3:
        pop_name = sys.argv[3]
    if len(sys.argv) > 4:
        nb_inds_simu = int(sys.argv[4])
    if len(sys.argv) > 5:
        nb_iterations_simu = int(sys.argv[5])

    logger.info(f"Starting error impact testing for project: {name_project}")
    logger.info(f"Selection: {name_selection}, Population: {pop_name}")
    logger.info(
        f"Simulation parameters: {nb_inds_simu} individuals, {nb_iterations_simu} iterations"
    )

    try:
        main_error_recaptures(
            name_project, name_selection, pop_name, nb_inds_simu, nb_iterations_simu
        )
    except Exception as e:
        logger.error(f"Error during testing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    cli_error_tests()
