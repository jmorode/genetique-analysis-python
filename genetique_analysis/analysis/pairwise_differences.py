import itertools
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..core.config import FileConfiguration
from ..utils.utils import (
    get_extension_if_subcat,
    select_a_given_pop,
    select_a_given_pop_year,
    select_a_given_pop_year_subcat,
    select_and_concat_all_genotypes,
)


def recover_pairwise_difference(config: FileConfiguration) -> None:
    data_pairwise = config.genotypes_data_pairwise.copy()
    selection = config.pops_to_select.copy()

    if config.agg_type == "all":
        df_selection = select_and_concat_all_genotypes(selection, data_pairwise)
        calculate_pairwise_differences(
            config, config.selection_name + "_all", df_selection
        )

    elif config.agg_type == "pops":
        for _pop in selection.Population.unique():
            df_selection = select_a_given_pop(data_pairwise, _pop)
            calculate_pairwise_differences(
                config, config.selection_name + f"_{_pop}", df_selection
            )

    elif config.agg_type == "pop_years":
        selection = selection[["Population", "Year"]].drop_duplicates()
        for _pop, _year in zip(selection.Population, selection.Year):
            df_selection = select_a_given_pop_year(data_pairwise, _pop, _year)
            calculate_pairwise_differences(
                config, config.selection_name + f"_{_pop}_{_year}", df_selection
            )

    else:  # subcat
        for _pop, _year, _sub in zip(
            selection.Population, selection.Year, selection.Subcategory
        ):
            df_selection = select_a_given_pop_year_subcat(
                data_pairwise, _pop, _year, _sub
            )
            ext = get_extension_if_subcat(_sub)
            calculate_pairwise_differences(
                config, config.selection_name + f"_{_pop}_{_year}{ext}", df_selection
            )


def calculate_pairwise_differences_pandas(
    config: FileConfiguration, pop: str, df_data: pd.DataFrame
) -> None:
    df_data = df_data.drop(columns=["Population", "Year", "Subcategory"])
    # data.index = data.iloc[:, 0]  # Set row names
    # data = data.iloc[:, 1:]  # Remove samples names once in index`

    # Convert to NumPy array
    data = df_data.iloc[:, 1:].to_numpy()
    # Number of individuals and loci
    num_individuals, num_loci = data.shape

    # List to store all pairwise differences
    all_differences = []

    # Create a difference matrix
    diff_matrix = np.zeros((num_individuals, num_individuals))

    # Vectorized calculation of differences
    for i in range(num_individuals):
        for j in range(i):
            nb_diff_pair = []
            for k in range(0, num_loci, 2):
                individual1_lock = data[i, k : k + 2].tolist()
                individual2_lock = data[j, k : k + 2].tolist()

                # Check if either individual has missing data at this locus
                # Missing data is represented as (0, 0) or NaN values
                ind1_missing = (
                    (individual1_lock[0] == 0 and individual1_lock[1] == 0)
                    or pd.isna(individual1_lock[0])
                    or pd.isna(individual1_lock[1])
                )
                ind2_missing = (
                    (individual2_lock[0] == 0 and individual2_lock[1] == 0)
                    or pd.isna(individual2_lock[0])
                    or pd.isna(individual2_lock[1])
                )

                # Skip this locus if either individual has missing data
                if ind1_missing or ind2_missing:
                    continue

                # Calculate differences only for loci with known alleles for both individuals
                if len(np.unique(individual1_lock + individual2_lock)) == 1:
                    nb_diff_loc = 0
                else:
                    nb_diff_loc = 2 - len(
                        list(set(individual1_lock) & set(individual2_lock))
                    )

                nb_diff_pair.append(nb_diff_loc)

            # Update difference matrix and all differences list
            # Sum of differences across all valid (non-missing) loci
            nb_diff_pair_sum = np.sum(nb_diff_pair)
            diff_matrix[i, j] = nb_diff_pair_sum
            diff_matrix[j, i] = nb_diff_pair_sum  # Symmetric for upper triangle
            all_differences.append(nb_diff_pair_sum)

    # Convert back to DataFrame and replace diagonal with NaN
    diff_matrix = pd.DataFrame(
        diff_matrix, index=df_data.iloc[:, 0], columns=df_data.iloc[:, 0]
    )
    np.fill_diagonal(diff_matrix.values, np.nan)

    # Save pairwise matrix
    diff_matrix.to_csv(
        f"{config.output_path}/pairwise_differences/pairwise_differences_{pop}.csv",
        sep=";",
        index=True,
    )

    # Calculate frequency of each number of differences
    freq_table = pd.Series(all_differences).value_counts(normalize=True).reset_index()
    freq_table = freq_table.rename(columns={"index": "nb_diff"})
    freq_table.to_csv(
        f"{config.output_path}/pairwise_differences/frequency_pairwise_distances_{pop}.csv",
        sep=";",
        index=False,
    )


def calculate_pairwise_differences(
    config: FileConfiguration, pop: str, df_data: pd.DataFrame
) -> pd.DataFrame:
    data = df_data.iloc[:, 4:].to_numpy()  # Convert to NumPy array
    num_individuals, num_loci = data.shape

    # Create a difference matrix
    diff_matrix = np.zeros((num_individuals, num_individuals))

    # List to store all pairwise differences
    all_differences = []

    # Vectorized calculation of differences
    for i in range(num_individuals):
        for j in range(i):
            nb_diff_pair = []
            for k in range(0, num_loci, 2):
                individual1_lock = data[i, k : k + 2].tolist()
                individual2_lock = data[j, k : k + 2].tolist()

                # Check if either individual has missing data at this locus
                # Missing data is represented as (0, 0) or NaN values
                ind1_missing = (
                    (individual1_lock[0] == 0 and individual1_lock[1] == 0)
                    or pd.isna(individual1_lock[0])
                    or pd.isna(individual1_lock[1])
                )
                ind2_missing = (
                    (individual2_lock[0] == 0 and individual2_lock[1] == 0)
                    or pd.isna(individual2_lock[0])
                    or pd.isna(individual2_lock[1])
                )

                # Skip this locus if either individual has missing data
                if ind1_missing or ind2_missing:
                    continue

                # Calculate differences only for loci with known alleles for both individuals
                if len(np.unique(individual1_lock + individual2_lock)) == 1:
                    nb_diff_loc = 0
                else:
                    nb_diff_loc = 2 - len(
                        list(set(individual1_lock) & set(individual2_lock))
                    )

                nb_diff_pair.append(nb_diff_loc)

            # Update difference matrix and all differences list
            # Sum of differences across all valid (non-missing) loci
            nb_diff_pair_sum = np.sum(nb_diff_pair)
            diff_matrix[i, j] = nb_diff_pair_sum
            diff_matrix[j, i] = nb_diff_pair_sum
            all_differences.append(nb_diff_pair_sum)

    # Convert back to DataFrame and replace diagonal with NaN
    diff_matrix = pd.DataFrame(
        diff_matrix, index=df_data.iloc[:, 0], columns=df_data.iloc[:, 0]
    )
    np.fill_diagonal(diff_matrix.values, np.nan)

    # Save pairwise matrix
    diff_matrix["Population"] = diff_matrix.index.map(config.dict_pop_samples)
    diff_matrix.to_csv(
        f"{config.output_path}/pairwise_differences/pairwise_differences_{pop}.csv",
        sep=";",
        index=True,
    )

    # Calculate frequency of each number of differences
    freq_table = pd.Series(all_differences).value_counts(normalize=True).reset_index()
    freq_table = freq_table.rename(columns={"index": "nb_diff"})
    freq_table.to_csv(
        f"{config.output_path}/pairwise_differences/frequency_pairwise_distances_{pop}.csv",
        sep=";",
        index=False,
    )
    return diff_matrix.drop(columns=["Population"])


def get_combination_duets(list_pops: list[str]) -> list:
    return list(itertools.combinations(list_pops, 2))


def separate_elements_with__string(x: str) -> Tuple[str, str]:
    return x.split("__")[0], x.split("__")[1]


def plot_one_pairwise_distribution(config: FileConfiguration, pop: str) -> None:
    freq_table = pd.read_csv(
        f"{config.output_path}/pairwise_differences/frequency_pairwise_distances_{pop}.csv",
        sep=";",
    )
    plt.figure(figsize=(10, 6))
    plt.bar(freq_table.nb_diff, freq_table.proportion, color="skyblue", ec="k")
    plt.xlabel("Number of differences")
    plt.ylabel("Frequency")
    plt.title(pop)
    plt.xlim(0, len(config.loci_list))
    plt.xticks(
        list(np.linspace(0, 2 * len(config.loci_list), 2 * len(config.loci_list) + 1))
    )
    plt.tight_layout()
    plt.savefig(
        f"{config.output_path}/pairwise_differences/plots/plot_pairwise_frequencies_{pop}.pdf"
    )
    plt.close()


def plot_two_pairwise_distributions(
    config: FileConfiguration, pop1: str, pop2: str
) -> None:
    # get distances
    freq_table1 = pd.read_csv(
        f"{config.output_path}/pairwise_differences/frequency_pairwise_distances_{pop1}.csv",
        sep=";",
    )
    freq_table2 = pd.read_csv(
        f"{config.output_path}/pairwise_differences/frequency_pairwise_distances_{pop2}.csv",
        sep=";",
    )

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.bar(
        freq_table2.nb_diff - 0.2,
        freq_table2.proportion,
        width=0.4,
        color="darkgray",
        alpha=0.4,
        label="Theoretical",
    )
    ax.bar(
        freq_table1.nb_diff + 0.2,
        freq_table1.proportion,
        width=0.4,
        color="blue",
        alpha=0.4,
        label="Measured",
    )
    ax.set_xlabel("Number of differences")
    ax.set_ylabel("Frequency")
    ax.set_xlim(0, len(config.loci_list))
    ax.set_xticks(
        list(np.linspace(0, 2 * len(config.loci_list), 2 * len(config.loci_list) + 1))
    )
    ax.legend()
    ax.set_title(pop1 + " - " + pop2)
    plt.tight_layout()
    plt.savefig(
        f"{config.output_path}/pairwise_differences/plots/plot_pairwise_frequencies_{pop1}_{pop2}.pdf"
    )
    plt.close()


def plot_intra_population_distances(config: FileConfiguration, pop1__pop2: str) -> None:
    pop1, pop2 = separate_elements_with__string(pop1__pop2)
    print(pop1, pop2)

    # calculate pairwise between the pops
    df_selection = config.genotypes_data_pairwise[
        config.genotypes_data_pairwise.Population.isin([pop1, pop2])
    ]
    calculate_pairwise_differences(
        config, config.selection_name + f"_{pop1__pop2}", df_selection
    )
    diff_matrix = pd.read_csv(
        f"{config.output_path}/pairwise_differences/pairwise_differences_{config.selection_name}_{pop1__pop2}.csv",
        sep=";",
    )
    diff_matrix["pop_cat"] = diff_matrix.Sample.apply(
        lambda x: config.dict_pop_samples[x]
    )
    list_ind_pop1 = diff_matrix[diff_matrix["pop_cat"] == pop1]["Sample"].to_list()
    list_ind_pop2 = diff_matrix[diff_matrix["pop_cat"] == pop2]["Sample"].to_list()
    diff_matrix = diff_matrix.set_index("Sample")

    # pop_cat -> Name of the population
    cross_part1 = diff_matrix[list_ind_pop1][diff_matrix["pop_cat"] == pop2].copy()
    cross_part2 = diff_matrix[list_ind_pop2][diff_matrix["pop_cat"] == pop1].copy()
    list_cross_part1 = cross_part1.values.flatten().tolist()
    list_cross_part2 = cross_part2.values.flatten().tolist()
    list_complete = list_cross_part1 + list_cross_part2

    freq_table = pd.Series(list_complete).value_counts(normalize=True).reset_index()
    freq_table = freq_table.rename(columns={"index": "nb_diff"})
    freq_table.to_csv(
        f"{config.output_path}/pairwise_differences/frequency_pairwise_distances_{pop1__pop2}.csv",
        sep=";",
        index=False,
    )
    plot_one_pairwise_distribution(config, pop1__pop2)


def plot_pairwise_distances_one_by_one(config: FileConfiguration) -> None:
    if config.agg_type == "all":
        plot_one_pairwise_distribution(config, config.selection_name + "_all")
    elif config.agg_type == "pops":
        for _pop in config.pops_to_select.Population.unique():
            plot_one_pairwise_distribution(config, config.selection_name + f"_{_pop}")
    elif config.agg_type == "pop_years":
        selection = config.pops_to_select[["Population", "Year"]].drop_duplicates()
        for _pop, _year in zip(selection.Population, selection.Year):
            plot_one_pairwise_distribution(
                config, config.selection_name + f"_{_pop}_{_year}"
            )
    else:  # subcat
        for _pop, _year, _sub in zip(
            config.pops_to_select.Population,
            config.pops_to_select.Year,
            config.pops_to_select.Subcategory,
        ):
            ext = get_extension_if_subcat(_sub)
            plot_one_pairwise_distribution(
                config, config.selection_name + f"_{_pop}_{_year}{ext}"
            )


def plot_pairwise_distances_two_by_two(config: FileConfiguration) -> None:
    if config.agg_type == "all":
        print(
            "'all' aggregation type selected - only got one distribution, so cannot plot two."
        )

    elif config.agg_type == "pops":
        possible_combinations = get_combination_duets(
            config.pops_to_select.Population.unique()
        )
        for _pop1, _pop2 in possible_combinations:
            plot_two_pairwise_distributions(
                config,
                config.selection_name + f"_{_pop1}",
                config.selection_name + f"_{_pop2}",
            )
    elif config.agg_type == "pop_years":
        selection = config.pops_to_select[["Population", "Year"]].drop_duplicates()
        selection["pop_year"] = selection.apply(
            lambda x: x.Population + "__" + str(x.Year), axis=1
        )
        possible_combinations = get_combination_duets(selection.pop_year.unique())
        for _ele1, _ele2 in possible_combinations:
            _pop1, _year1 = separate_elements_with__string(_ele1)
            _pop2, _year2 = separate_elements_with__string(_ele2)
            plot_two_pairwise_distributions(
                config,
                config.selection_name + f"_{_pop1}_{_year1}",
                config.selection_name + f"_{_pop2}_{_year2}",
            )

    else:  # subcat
        selection = config.pops_to_select.copy()
        selection["pop_year_sub"] = selection.apply(
            lambda x: x.Population
            + "__"
            + str(x.Year)
            + get_extension_if_subcat(x.Subcategory),
            axis=1,
        )
        possible_combinations = get_combination_duets(selection.pop_year_sub.unique())
        for _ele1, _ele2 in possible_combinations:
            _pop1, _year_ext1 = separate_elements_with__string(_ele1)
            _pop2, _year_ext2 = separate_elements_with__string(_ele2)
            plot_two_pairwise_distributions(
                config,
                config.selection_name + f"_{_pop1}_{_year_ext1}",
                config.selection_name + f"_{_pop2}_{_year_ext2}",
            )


def plot_pairwise_distances_all_intra_populations(config: FileConfiguration) -> None:
    if config.agg_type == "pops":
        possible_combinations = get_combination_duets(
            config.pops_to_select.Population.unique()
        )
        for _pop1, _pop2 in possible_combinations:
            plot_intra_population_distances(config, f"{_pop1}__{_pop2}")
    else:
        print(
            "plot_pairwise_distances_all_intra_populations - available only for 'pops' aggregation type"
        )
