from config import FileConfiguration
from description import compute_frequency_and_heterozygosity, save_proba_same_ind
from loguru import logger
from ml_relate import (
    count_nbr_relationships_from_ml_relate_output,
    reliability_ml_relate_based_simulated_families,
    write_ml_relate_input_file_genotypes,
    write_ml_relate_input_file_simu_families,
)
from pairwise_differences import recover_pairwise_difference, plot_pairwise_distances_one_by_one, plot_pairwise_distances_two_by_two, plot_pairwise_distances_all_intra_populations
from recaptures import get_recaptures
from generation import plot_cumulated_and_not_frequencies_for_simu_relationships_and_sample_by_pop, plot_pairwise_distances_per_relationships
from test_stats import compute_ks_test_between_pairwise_distance_and_relationships_simu


def main(
    project_name: str, selection_name: str, agg_type: str, pop_order_list: list[str]
) -> None:
    # INIT
    config = FileConfiguration(project_name, agg_type, selection_name)
    config.initialize_global()
    logger.info("Configuration and data retrieved.")

    # Description
    save_proba_same_ind(config)
    compute_frequency_and_heterozygosity(config, pop_order_list)

    # Pairwise
    recover_pairwise_difference(config)

    plot_pairwise_distances_one_by_one(config)
    if config.agg_type != "all":
        plot_pairwise_distances_two_by_two(config)
    if config.agg_type == "pops":
        plot_pairwise_distances_all_intra_populations(config)

    # Recaptures (need pairwise)
    get_recaptures(config)

    # ML Relate (test output in the software / and then the reliability function)
    # (needs frequencies)
    write_ml_relate_input_file_genotypes(config)
    write_ml_relate_input_file_simu_families(config, nb_families=1000)
    # add ml relate output to input folder with the name ml_relate_output_{selection_name}.csv
    count_nbr_relationships_from_ml_relate_output(config)
    reliability_ml_relate_based_simulated_families(config, nb_families=1000)

    # Pairwise Relationships
    plot_pairwise_distances_per_relationships(config, nb_families=1000)
    compute_ks_test_between_pairwise_distance_and_relationships_simu(config, nb_families=1000)
    if config.agg_type == "pops":
        plot_cumulated_and_not_frequencies_for_simu_relationships_and_sample_by_pop(config)

    logger.info("The End.")


if __name__ == "__main__":
    name_project = "test_sept_2025" #"analysis_nov_2024_more_loci"
    name_selection = "everything"
    # if aggregation_type == all: all config together
    # if aggregation_type == pops: all config together for each population (merge the years)
    # if aggregation_type == pop_years: each population / year
    # if aggregation_type == subcat: account for subcategories
    aggregation_type = "all"  # "all", "pops", "pop_years", "subcat"
    freq_plot_pop_order_list = (
        []
    )
    # ['SF Arroyo-Tampemole', 'SF Pozo Pachon Praxe', 'CF Pachon', 'CF Refugio', 'CF Sabinos', 'CF Tinaja','CF Piedras', 'CF Toro', 'CF Subterraneo ']
    main(name_project, name_selection, aggregation_type, freq_plot_pop_order_list)


    # add an option so that we can directly select all the genotype database (option in config file ?)
    # add analysis options with config file for analysis - group by relevant action
    # Exhaustive documentation and tests
