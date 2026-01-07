"""
Description module for genetics analysis.

This module handles frequency and heterozygosity calculations, including:
- Probability calculations for same individuals
- Allele frequency computations
- Heterozygosity calculations
- Frequency plotting and visualization
"""

from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..core.config import FileConfiguration
from ..core.constants import CSV_SEPARATOR, DEFAULT_FIGURE_SIZE
from ..utils.utils import (
    get_extension_if_subcat,
    select_a_given_pop,
    select_a_given_pop_year,
    select_a_given_pop_year_subcat,
    select_and_concat_all_genotypes,
)


def compute_proba_same_ind(
    config: FileConfiguration, df_allele: pd.DataFrame, name: str
) -> float:
    """
    Compute probability that two randomly chosen individuals are the same.

    This function calculates the probability that two randomly chosen individuals
    from a population have identical genotypes across all loci.

    Args:
        config: FileConfiguration object
        df_allele: DataFrame containing allele data
        name: Name identifier for the population/group

    Returns:
        Probability that two randomly chosen individuals are identical
    """
    # Remove metadata columns
    metadata_cols = ["Population", "Sample", "Year", "Subcategory"]
    df_allele = df_allele.drop(columns=metadata_cols)

    # Initialize probability
    proba = 1.0
    df_by_locus = pd.DataFrame()

    # Loop over loci
    for locus in df_allele.columns:
        # Count alleles (exclude missing values: 0 and NaN)
        # Filter out missing data (0 represents missing alleles)
        locus_data = df_allele[locus].dropna()
        locus_data = locus_data[locus_data != 0]  # Exclude 0 values (missing data)
        loci_alleles = locus_data.value_counts()

        if len(loci_alleles) == 0:
            continue  # Skip loci with no data

        # Calculate term 1: Sum of allele frequencies to the power of 4
        allele_freqs = loci_alleles / loci_alleles.sum()
        term_1 = sum(allele_freqs**4)

        # Calculate term 2: Sum of squared products of allele frequencies
        term_2 = 0
        if len(loci_alleles) > 1:
            for i in range(len(allele_freqs) - 1):
                for j in range(i + 1, len(allele_freqs)):
                    term_2 += (allele_freqs.iloc[i] * allele_freqs.iloc[j]) ** 2

        # Update probability for this locus
        locus_proba = term_1 + 4 * term_2
        proba *= locus_proba

        # Keep info by locus
        df_by_locus = pd.concat(
            [
                df_by_locus,
                pd.DataFrame(
                    {"pop": name, "locus": locus, "proba": locus_proba}, index=[0]
                ),
            ],
            ignore_index=True,
        )

    # Save locus-specific probabilities
    output_file = (
        f"{config.output_path}/heterozygosity/probability_same_inds_by_loci_{name}.csv"
    )
    df_by_locus.to_csv(output_file, index=False, sep=CSV_SEPARATOR)

    return proba


def compute_proba_same_ind_siblings(
    config: FileConfiguration, df_allele: pd.DataFrame, name: str
) -> float:
    """
    Compute probability that two randomly chosen individuals with the same genotype are siblings.

    This function calculates the probability that two randomly chosen individuals
    from a population with identical genotypes across all loci are siblings.

    Args:
        config: FileConfiguration object
        df_allele: DataFrame containing allele data
        name: Name identifier for the population/group

    Returns:
        Probability that two randomly chosen individuals with the same genotypes are siblings
    """
    # Remove metadata columns
    metadata_cols = ["Population", "Sample", "Year", "Subcategory"]
    df_allele = df_allele.drop(columns=metadata_cols)

    # Initialize probability
    proba = 1.0
    df_by_locus = pd.DataFrame()

    # Loop over loci
    for locus in df_allele.columns:
        # Count alleles (exclude missing values: 0 and NaN)
        # Filter out missing data (0 represents missing alleles)
        locus_data = df_allele[locus].dropna()
        locus_data = locus_data[locus_data != 0]  # Exclude 0 values (missing data)
        loci_alleles = locus_data.value_counts()

        if len(loci_alleles) == 0:
            continue  # Skip loci with no data

        # Calculate sum of allele frequencies to the power of 2 and to the power of 4
        allele_freqs = loci_alleles / loci_alleles.sum()
        term_squared = sum(allele_freqs**2)
        term_power4 = sum(allele_freqs**4)

        # Update probability for this locus
        locus_proba = (
            0.25
            + (0.5 * term_squared)
            + (0.5 * np.power(term_squared, 2))
            - (0.25 * term_power4)
        )
        proba *= locus_proba

        # Keep info by locus
        df_by_locus = pd.concat(
            [
                df_by_locus,
                pd.DataFrame(
                    {"pop": name, "locus": locus, "proba": locus_proba}, index=[0]
                ),
            ],
            ignore_index=True,
        )

    # Save locus-specific probabilities
    output_file = f"{config.output_path}/heterozygosity/probability_same_inds_by_loci_sibilings_{name}.csv"
    df_by_locus.to_csv(output_file, index=False, sep=CSV_SEPARATOR)

    return proba


def _compute_proba_all_samples(
    config: FileConfiguration, selection: pd.DataFrame, selection_name: str
) -> pd.DataFrame:
    """Compute probability for all samples combined."""
    df_selection = select_and_concat_all_genotypes(
        selection, config.genotypes_two_lines
    )
    return pd.DataFrame(
        {
            "Group": selection_name,
            "ProbabilitySameInd": compute_proba_same_ind(config, df_selection, "all"),
        },
        index=[0],
    )


def _compute_proba_by_populations(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    """Compute probability for each population separately."""
    df_probas = pd.DataFrame()
    for pop in selection.Population.unique():
        selec = select_a_given_pop(config.genotypes_two_lines, pop)
        df = pd.DataFrame(
            {
                "Population": pop,
                "ProbabilitySameInd": compute_proba_same_ind(config, selec, pop),
            },
            index=[0],
        )
        df_probas = pd.concat([df, df_probas], ignore_index=True)
    return df_probas


def _compute_proba_by_pop_years(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    """Compute probability for each population-year combination."""
    df_probas = pd.DataFrame()
    selection = selection[["Population", "Year"]].drop_duplicates()
    for pop, year in zip(selection.Population, selection.Year):
        selec = select_a_given_pop_year(config.genotypes_two_lines, pop, year)
        df = pd.DataFrame(
            {
                "Population": pop,
                "Year": year,
                "ProbabilitySameInd": compute_proba_same_ind(
                    config, selec, f"{pop}_{year}"
                ),
            },
            index=[0],
        )
        df_probas = pd.concat([df, df_probas], ignore_index=True)
    return df_probas


def _compute_proba_by_subcategories(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    """Compute probability for each population-year-subcategory combination."""
    df_probas = pd.DataFrame()
    for pop, year, sub in zip(
        selection.Population, selection.Year, selection.Subcategory
    ):
        selec = select_a_given_pop_year_subcat(
            config.genotypes_two_lines, pop, year, sub
        )
        df = pd.DataFrame(
            {
                "Population": pop,
                "Year": year,
                "Subcategory": sub,
                "ProbabilitySameInd": compute_proba_same_ind(
                    config, selec, f"{pop}_{year}_{sub}"
                ),
            },
            index=[0],
        )
        df_probas = pd.concat([df, df_probas], ignore_index=True)
    return df_probas


def _compute_proba_siblings_all_samples(
    config: FileConfiguration, selection: pd.DataFrame, selection_name: str
) -> pd.DataFrame:
    """Compute probability for all samples combined (siblings)."""
    df_selection = select_and_concat_all_genotypes(
        selection, config.genotypes_two_lines
    )
    return pd.DataFrame(
        {
            "Group": selection_name,
            "ProbabilitySameIndSiblings": compute_proba_same_ind_siblings(
                config, df_selection, "all"
            ),
        },
        index=[0],
    )


def _compute_proba_siblings_by_populations(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    """Compute probability for each population separately (siblings)."""
    df_probas = pd.DataFrame()
    for pop in selection.Population.unique():
        selec = select_a_given_pop(config.genotypes_two_lines, pop)
        df = pd.DataFrame(
            {
                "Population": pop,
                "ProbabilitySameIndSiblings": compute_proba_same_ind_siblings(
                    config, selec, pop
                ),
            },
            index=[0],
        )
        df_probas = pd.concat([df, df_probas], ignore_index=True)
    return df_probas


def _compute_proba_siblings_by_pop_years(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    """Compute probability for each population-year combination (siblings)."""
    df_probas = pd.DataFrame()
    selection = selection[["Population", "Year"]].drop_duplicates()
    for pop, year in zip(selection.Population, selection.Year):
        selec = select_a_given_pop_year(config.genotypes_two_lines, pop, year)
        df = pd.DataFrame(
            {
                "Population": pop,
                "Year": year,
                "ProbabilitySameIndSiblings": compute_proba_same_ind_siblings(
                    config, selec, f"{pop}_{year}"
                ),
            },
            index=[0],
        )
        df_probas = pd.concat([df, df_probas], ignore_index=True)
    return df_probas


def _compute_proba_siblings_by_subcategories(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    """Compute probability for each population-year-subcategory combination (siblings)."""
    df_probas = pd.DataFrame()
    for pop, year, sub in zip(
        selection.Population, selection.Year, selection.Subcategory
    ):
        selec = select_a_given_pop_year_subcat(
            config.genotypes_two_lines, pop, year, sub
        )
        df = pd.DataFrame(
            {
                "Population": pop,
                "Year": year,
                "Subcategory": sub,
                "ProbabilitySameIndSiblings": compute_proba_same_ind_siblings(
                    config, selec, f"{pop}_{year}_{sub}"
                ),
            },
            index=[0],
        )
        df_probas = pd.concat([df, df_probas], ignore_index=True)
    return df_probas


def save_proba_same_ind_siblings(config: FileConfiguration) -> None:
    """
    Save probability calculations for same individuals (siblings) based on aggregation type.

    Args:
        config: FileConfiguration object containing analysis parameters
    """
    selection = config.pops_to_select.copy()
    selection_name = config.selection_name

    # Compute probabilities based on aggregation type
    if config.agg_type == "all":
        df_probas = _compute_proba_siblings_all_samples(
            config, selection, selection_name
        )
        selection_name += "_all"
    elif config.agg_type == "pops":
        df_probas = _compute_proba_siblings_by_populations(config, selection)
        selection_name += "_pops"
    elif config.agg_type == "pop_years":
        df_probas = _compute_proba_siblings_by_pop_years(config, selection)
        selection_name += "_pops_years"
    else:  # subcat
        df_probas = _compute_proba_siblings_by_subcategories(config, selection)
        selection_name += "_subcat"

    # Save results
    output_file = f"{config.output_path}/heterozygosity/probability_same_inds_siblings_{selection_name}.csv"
    df_probas.to_csv(output_file, index=False, sep=CSV_SEPARATOR)


def save_proba_same_ind(config: FileConfiguration) -> None:
    """
    Save probability calculations for same individuals based on aggregation type.

    Args:
        config: FileConfiguration object containing analysis parameters
    """
    selection = config.pops_to_select.copy()
    selection_name = config.selection_name

    # Compute probabilities based on aggregation type
    if config.agg_type == "all":
        df_probas = _compute_proba_all_samples(config, selection, selection_name)
        selection_name += "_all"
    elif config.agg_type == "pops":
        df_probas = _compute_proba_by_populations(config, selection)
        selection_name += "_pops"
    elif config.agg_type == "pop_years":
        df_probas = _compute_proba_by_pop_years(config, selection)
        selection_name += "_pops_years"
    else:  # subcat
        df_probas = _compute_proba_by_subcategories(config, selection)
        selection_name += "_subcat"

    # Save results
    output_file = f"{config.output_path}/heterozygosity/probability_same_inds_{selection_name}.csv"
    df_probas.to_csv(output_file, index=False, sep=CSV_SEPARATOR)


def calculate_frequencies(
    config: FileConfiguration,
    selection_genotypes_two_lines: pd.DataFrame,
    selection_name: str,
) -> pd.DataFrame:
    """
    Calculate allele frequencies for each locus in the selection.

    Args:
        config: FileConfiguration object
        selection_genotypes_two_lines: DataFrame containing genotype data
        selection_name: Name identifier for the selection

    Returns:
        DataFrame with allele frequencies for each locus
    """
    df = pd.DataFrame()

    for locus in config.loci_list:
        # Count alleles and calculate frequencies (exclude missing values: 0 and NaN)
        # Filter out missing data (0 represents missing alleles)
        locus_data = selection_genotypes_two_lines[locus].dropna()
        locus_data = locus_data[locus_data != 0]  # Exclude 0 values (missing data)
        allele_counts = locus_data.value_counts()
        if len(allele_counts) == 0:
            continue  # Skip loci with no data

        allele_freqs = allele_counts / allele_counts.sum()

        # Create DataFrame for this locus
        df_locus = pd.DataFrame(
            {"frequency": allele_freqs, "locus": locus}
        ).reset_index()
        df_locus = df_locus.rename(columns={locus: "allele"})

        df = pd.concat([df, df_locus], ignore_index=True)

    # Format alleles as 3-digit strings
    df["allele"] = df["allele"].apply(lambda x: f"{x:03d}")

    # Save frequencies
    output_file = f"{config.output_path}/raw_data/frequencies_{config.project_name}_{selection_name}.csv"
    df.to_csv(output_file, index=False, sep=CSV_SEPARATOR)

    return df


def heterozygosity_calculation(
    config: FileConfiguration,
    df_frequencies: pd.DataFrame,
    selection_name: str,
) -> Tuple[pd.DataFrame, float]:
    """
    Calculate expected heterozygosity for each locus and globally.

    Args:
        config: FileConfiguration object
        df_frequencies: DataFrame containing allele frequencies
        selection_name: Name identifier for the selection

    Returns:
        Tuple of (DataFrame with locus heterozygosity, global heterozygosity value)
    """
    he_distr = pd.DataFrame()
    for locus in df_frequencies["locus"].unique():
        locus_data = df_frequencies[df_frequencies["locus"] == locus][
            "frequency"
        ].values
        he_loc_sums = np.sum([np.power(freq, 2) for freq in locus_data])
        he_value = 1 - he_loc_sums

        he_distr = pd.concat(
            [he_distr, pd.DataFrame({"locus": locus, "he": he_value}, index=[0])]
        )

    # Save heterozygosity table
    output_file = (
        f"{config.output_path}/heterozygosity/TableHeterozygosity_{selection_name}.csv"
    )
    he_distr.to_csv(output_file, sep=CSV_SEPARATOR, index=False)

    # Calculate global heterozygosity
    n_loci = len(he_distr.locus.unique())
    he_global = he_distr.he.sum() / n_loci
    return he_distr, he_global


def plot_heterozygosity(
    config: FileConfiguration, he_distr: pd.DataFrame, selection_name: str
) -> None:
    """
    Plot expected heterozygosity for each locus.

    Args:
        config: FileConfiguration object
        he_distr: DataFrame containing heterozygosity values by locus
        selection_name: Name identifier for the selection
    """
    he_distr = he_distr.sort_values(by="he", ascending=False)

    plt.figure(figsize=DEFAULT_FIGURE_SIZE)
    plt.bar(he_distr.locus, he_distr.he, color="skyblue", ec="k")
    plt.xlabel("Locus")
    plt.ylabel("Expected Heterozygosity (H0)")
    plt.title(f"Expected Heterozygosity for each Locus in {selection_name}")
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, 1)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()

    output_file = (
        f"{config.output_path}/heterozygosity/heterozygosity_{selection_name}.pdf"
    )
    plt.savefig(output_file)
    plt.close()


def frequency_and_heterozygosity_all_samples(
    config: FileConfiguration, selection: pd.DataFrame, selection_name: str
) -> pd.DataFrame:
    df_selection = select_and_concat_all_genotypes(
        selection, config.genotypes_two_lines
    )
    df_frequencies = calculate_frequencies(
        config, df_selection, selection_name + "_all"
    )
    he_distr, he_value = heterozygosity_calculation(
        config, df_frequencies, selection_name
    )
    df_heterozygosity = pd.DataFrame(
        {"Group": selection_name, "Heterozygosity": he_value}, index=[0]
    )
    plot_heterozygosity(config, he_distr, selection_name)
    return df_heterozygosity


def frequency_and_heterozygosity_pops(
    config: FileConfiguration, selection: pd.DataFrame, selection_name: str
) -> pd.DataFrame:
    df_heterozygosity = pd.DataFrame()
    for _pop in selection.Population.unique():
        selec = select_a_given_pop(config.genotypes_two_lines, _pop)
        df_frequencies = calculate_frequencies(
            config, selec, config.selection_name + f"_{_pop}"
        )
        he_distr, he_value = heterozygosity_calculation(
            config, df_frequencies, selection_name
        )
        df = pd.DataFrame({"Population": _pop, "Heterozygosity": he_value}, index=[0])
        df_heterozygosity = pd.concat([df, df_heterozygosity], ignore_index=True)
        plot_heterozygosity(config, he_distr, _pop)
    return df_heterozygosity


def frequency_and_heterozygosity_pop_years(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    df_heterozygosity = pd.DataFrame()
    selection = selection[["Population", "Year"]].drop_duplicates()
    for _pop, _year in zip(selection.Population, selection.Year):
        selec = select_a_given_pop_year(config.genotypes_two_lines, _pop, _year)
        df_frequencies = calculate_frequencies(
            config, selec, config.selection_name + f"_{_pop}_{_year}"
        )
        he_distr, he_value = heterozygosity_calculation(
            config, df_frequencies, f"{_pop}_{_year}"
        )

        df = pd.DataFrame(
            {"Population": _pop, "Year": _year, "Heterozygosity": he_value}, index=[0]
        )
        df_heterozygosity = pd.concat([df, df_heterozygosity], ignore_index=True)
        plot_heterozygosity(config, he_distr, f"{_pop}_{_year}")
    return df_heterozygosity


def frequency_and_heterozygosity_subcat(
    config: FileConfiguration, selection: pd.DataFrame
) -> pd.DataFrame:
    df_heterozygosity = pd.DataFrame()
    for _pop, _year, _sub in zip(
        selection.Population, selection.Year, selection.Subcategory
    ):
        selec = select_a_given_pop_year_subcat(
            config.genotypes_two_lines, _pop, _year, _sub
        )
        _ext = get_extension_if_subcat(_sub)
        df_frequencies = calculate_frequencies(
            config, selec, config.selection_name + f"_{_pop}_{_year}{_ext}"
        )
        he_distr, he_value = heterozygosity_calculation(
            config, df_frequencies, f"{_pop}_{_year}{_ext}"
        )

        df = pd.DataFrame(
            {
                "Population": _pop,
                "Year": _year,
                "Subcategory": _sub,
                "Heterozygosity": he_value,
            },
            index=[0],
        )
        df_heterozygosity = pd.concat([df, df_heterozygosity], ignore_index=True)
        plot_heterozygosity(config, he_distr, f"{_pop}_{_year}{_ext}")
    return df_heterozygosity


def plot_frequencies(
    config: FileConfiguration,
    selection: pd.DataFrame,
    selection_name: str,
    agg_type: str,
    order_list: list[str],
) -> None:

    df_pops = pd.DataFrame()
    if agg_type == "all":
        df_pops = pd.read_csv(
            f"{config.output_path}/raw_data/frequencies_{config.project_name}_{selection_name}_all.csv",
            sep=";",
        )
        df_pops["population"] = "all"
        df_pops["legend"] = selection_name

    elif agg_type == "pops":
        for _pop in selection.Population.unique():
            df = pd.read_csv(
                f"{config.output_path}/raw_data/frequencies_{config.project_name}_{selection_name}_{_pop}.csv",
                sep=";",
            )
            df["population"] = _pop
            df["legend"] = _pop
            df_pops = pd.concat([df_pops, df], ignore_index=True)

    elif agg_type == "pop_years":
        selection = selection[["Population", "Year"]].drop_duplicates()
        for _pop, _year in zip(selection.Population, selection.Year):
            df = pd.read_csv(
                f"{config.output_path}/raw_data/frequencies_{config.project_name}_{selection_name}_{_pop}_{_year}.csv",
                sep=";",
            )
            df["population"] = _pop
            df["year"] = _year
            df["legend"] = f"{_pop} {_year}"
            df_pops = pd.concat([df_pops, df], ignore_index=True)

    else:
        for _pop, _year, _sub in zip(
            selection.Population, selection.Year, selection.Subcategory
        ):
            _ext = get_extension_if_subcat(_sub)
            df = pd.read_csv(
                f"{config.output_path}/raw_data/frequencies_{config.project_name}_{selection_name}_{_pop}_{_year}{_ext}.csv",
                sep=";",
            )
            df["population"] = _pop
            df["year"] = _year
            df["subcategory"] = _sub
            df["legend"] = f"{_pop} {_year}{_ext}"
            df_pops = pd.concat([df_pops, df], ignore_index=True)

    df_pops = df_pops.sort_values(by=["locus", "allele", "population"])
    plot_index_row, plot_index_col = 0, 0
    nb_plots = len(df_pops.locus.unique())
    nbr_rows_plots = int(np.ceil(nb_plots / 3)) if int(np.ceil(nb_plots / 3)) > 1 else 2
    _, ax = plt.subplots(nbr_rows_plots, 3, figsize=(30, 10 * nbr_rows_plots))

    for _loc in df_pops.locus.unique():
        df_loc = df_pops[(df_pops.locus == _loc)][
            ["legend", "allele", "frequency"]
        ].copy()
        df_loc = df_loc.pivot_table(
            index=["allele"], values=["frequency"], columns=["legend"], sort=False
        )
        df_loc.columns = df_loc.columns.droplevel(0)
        if len(order_list) > 0:
            df_loc = df_loc[order_list]
        df_loc.plot.bar(ax=ax[plot_index_col, plot_index_row])
        ax[plot_index_col, plot_index_row].set_title(_loc)

        if plot_index_row < 2:
            plot_index_row += 1
        else:
            plot_index_row = 0
            plot_index_col += 1

    # remove axis of empty plots
    if nb_plots % 3 != 0:
        while plot_index_row < 3 and plot_index_col < nbr_rows_plots:
            ax[plot_index_col, plot_index_row].set_axis_off()
            plot_index_row += 1

    plt.tight_layout()
    plt.savefig(f"{config.output_path}/raw_data/frequencies_{selection_name}.pdf")
    plt.close()

    # Save summary table
    df_pops = df_pops.pivot_table(
        index=["locus", "allele"], values=["frequency"], columns=["legend"], sort=False
    )
    df_pops.columns = df_pops.columns.droplevel(0)
    if len(order_list) > 0:
        df_pops = df_pops[order_list]
    df_pops = df_pops.reset_index()
    df_pops.to_csv(
        f"{config.output_path}/raw_data/frequencies_summary_table_{selection_name}.csv",
        sep=";",
        index=False,
    )


def compute_frequency_and_heterozygosity(
    config: FileConfiguration, pop_order_list: list[str]
) -> None:
    selection = config.pops_to_select.copy()
    selection_name = config.selection_name
    if config.agg_type == "all":
        df_heterozygosity = frequency_and_heterozygosity_all_samples(
            config, selection, selection_name
        )

    elif config.agg_type == "pops":
        df_heterozygosity = frequency_and_heterozygosity_pops(
            config, selection, selection_name
        )

    elif config.agg_type == "pop_years":
        df_heterozygosity = frequency_and_heterozygosity_pop_years(config, selection)

    else:  # subcat
        df_heterozygosity = frequency_and_heterozygosity_subcat(config, selection)

    # selection_name += "_" + config.agg_type
    plot_frequencies(config, selection, selection_name, config.agg_type, pop_order_list)

    # save
    df_heterozygosity.to_csv(
        f"{config.output_path}/heterozygosity/heterozygosity_values_{selection_name}.csv",
        sep=";",
        index=False,
    )
