"""
Analysis modules for genetique analysis.

This module contains all the analysis functions for genetics data processing
including frequency calculations, pairwise analysis, ML-related analysis,
and statistical testing.
"""

# Import modules lazily to avoid circular import issues
try:
    from .description import (
        calculate_frequencies,
        compute_frequency_and_heterozygosity,
        compute_proba_same_ind,
        heterozygosity_calculation,
        plot_frequencies,
        plot_heterozygosity,
        save_proba_same_ind,
    )
    from .pairwise_differences import (
        calculate_pairwise_differences,
        calculate_pairwise_differences_pandas,
        get_combination_duets,
        plot_intra_population_distances,
        plot_one_pairwise_distribution,
        plot_pairwise_distances_all_intra_populations,
        plot_pairwise_distances_one_by_one,
        plot_pairwise_distances_two_by_two,
        plot_two_pairwise_distributions,
        recover_pairwise_difference,
        separate_elements_with__string,
    )
    from .recaptures import get_recaptures, get_recaptures_sample
except ImportError:
    # If dependencies are not available, define placeholder functions
    def _placeholder(*args, **kwargs):
        raise ImportError(
            "Dependencies not installed. Please install with: pip install -e ."
        )

    # Define placeholder functions for all exports
    calculate_frequencies = _placeholder
    compute_frequency_and_heterozygosity = _placeholder
    compute_proba_same_ind = _placeholder
    heterozygosity_calculation = _placeholder
    plot_frequencies = _placeholder
    plot_heterozygosity = _placeholder
    save_proba_same_ind = _placeholder
    calculate_pairwise_differences = _placeholder
    calculate_pairwise_differences_pandas = _placeholder
    get_combination_duets = _placeholder
    plot_intra_population_distances = _placeholder
    plot_one_pairwise_distribution = _placeholder
    plot_pairwise_distances_all_intra_populations = _placeholder
    plot_pairwise_distances_one_by_one = _placeholder
    plot_pairwise_distances_two_by_two = _placeholder
    plot_two_pairwise_distributions = _placeholder
    recover_pairwise_difference = _placeholder
    separate_elements_with__string = _placeholder
    get_recaptures = _placeholder
    get_recaptures_sample = _placeholder

__all__ = [
    # Description module
    "calculate_frequencies",
    "compute_frequency_and_heterozygosity",
    "compute_proba_same_ind",
    "heterozygosity_calculation",
    "plot_frequencies",
    "plot_heterozygosity",
    "save_proba_same_ind",
    # Error generation module
    "compute_avg_q95_simu",
    "get_distrib_recap_errors",
    "get_nb_recaptures_following_error_introduction",
    "introduce_errors",
    "plot_distrib_recap_errors",
    "test_error_impact_scenarios_on_given_pop",
    # Generation module
    "generate_child",
    "generate_n_families_following_frequencies_for_ml_relate",
    "generate_n_families_following_frequencies_for_pairwise_distances",
    "generate_related_individuals_following_frequencies",
    "generate_unrelated_individuals_following_frequencies",
    "plot_cumulated_and_not_frequencies_for_simu_relationships_and_sample_by_pop",
    "plot_pairwise_distances_per_relationships",
    # ML relate module
    "add_missing_relationships_matches",
    "add_missing_relationships_mismatches",
    "compare_numpy_arrays",
    "count_nbr_relationships_from_ml_relate_output",
    "get_percentage_and_count_relationships_families",
    "get_sum_percentage_by_category",
    "get_total_percentage_over_families",
    "overall_stats_relationships_simulated_families",
    "reliability_ml_relate_based_simulated_families",
    "reformat_genotypes_selection_for_ml_relate",
    "write_ml_relate_input_file",
    "write_ml_relate_input_file_genotypes",
    "write_ml_relate_input_file_simu_families",
    # Pairwise differences module
    "calculate_pairwise_differences",
    "calculate_pairwise_differences_pandas",
    "get_combination_duets",
    "plot_intra_population_distances",
    "plot_one_pairwise_distribution",
    "plot_pairwise_distances_all_intra_populations",
    "plot_pairwise_distances_one_by_one",
    "plot_pairwise_distances_two_by_two",
    "plot_two_pairwise_distributions",
    "recover_pairwise_difference",
    "separate_elements_with__string",
    # Recaptures module
    "get_recaptures",
    "get_recaptures_sample",
    # Test stats module
    "compute_ks_test_between_pairwise_distance_and_relationships_simu",
    "compute_ks_test_between_pairwise_distance_and_relationships_simu_per_selection",
    "get_mean_std_for_ks_test_between_pairwise_distance_and_relationships_simu",
]
