from error_generation import test_error_impact_scenarios_on_given_pop
from config import FileConfiguration
from loguru import logger
from description import save_proba_same_ind


def main_error_recaptures(
    project_name: str, selection_name: str, pop_name: str, nb_inds_simu: int, nb_iterations_simu:int
) -> None:
    # INIT
    config = FileConfiguration(project_name, "pops", selection_name)
    config.initialize_global()
    logger.info("Configuration and data retrieved.")

    # Test impact of errors of genotyping in recapture rates
    save_proba_same_ind(config)
    test_error_impact_scenarios_on_given_pop(config, pop_name, nb_inds_simu, nb_iterations_simu)

    logger.info("The End.")


if __name__ == "__main__":
    name_project = "analysis_nov_2024_more_loci"
    name_selection = "pachon_sabinos"
    name_population = "CF_Pachon"
    nb_generated_inds = 202
    nb_evaluation_recap = 500

    main_error_recaptures(name_project, name_selection, name_population, nb_generated_inds, nb_evaluation_recap)

