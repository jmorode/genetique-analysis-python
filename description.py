from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from config import FileConfiguration
from utils import (
    select_a_given_pop,
    select_a_given_pop_year,
    select_a_given_pop_year_subcat,
    select_and_concat_all_genotypes,
    get_extension_if_subcat
)


def compute_proba_same_ind(config: FileConfiguration, df_allele: pd.DataFrame, name: str) -> float:
    df_allele = df_allele.drop(columns=["Population", "Sample", "Year", "Subcategory"])
    # Initialization of probability
    proba = 1.0
    # keep info for each locus
    df_by_locus = pd.DataFrame()
    # Loop over loci
    for loci in list(df_allele):
        # Count alleles (exclude missing values)
        loci_alleles = df_allele[loci].value_counts().dropna()

        # Calculate term 1: Sum of allele frequencies to the power of 4
        term_1 = sum((loci_alleles / loci_alleles.sum()) ** 4)

        # Calculate term 2: Sum of squared products of allele frequencies
        term_2 = 0
        if len(loci_alleles) > 1:
            freqs = loci_alleles / loci_alleles.sum()

            for i in range(len(freqs) - 1):
                for j in range(i + 1, len(freqs)):
                    term_2 += (freqs.iloc[i] * freqs.iloc[j]) ** 2

        # Update probability for this locus
        proba *= term_1 + 4 * term_2
        # Keep info by locus
        df_by_locus = pd.concat([df_by_locus,
                                 pd.DataFrame({'pop':name, 'locus': loci, 'proba': term_1 + 4 * term_2}, index=[0])],
                                ignore_index=True)
    df_by_locus.to_csv(
        f"{config.output_path}/heterozygosity/probability_same_inds_by_loci_{name}.csv",
        index=False,
        sep=";"
    )
    return proba


def save_proba_same_ind(config: FileConfiguration) -> None:
    selection = config.pops_to_select.copy()
    selection_name = config.selection_name

    if config.agg_type == "all":
        df_selection = select_and_concat_all_genotypes(
            selection, config.genotypes_two_lines
        )
        df_probas = pd.DataFrame(
            {
                "Group": selection_name,
                "ProbabilitySameInd": compute_proba_same_ind(config, df_selection, "all"),
            },
            index=[0],
        )
        selection_name += "_all"

    elif config.agg_type == "pops":
        df_probas = pd.DataFrame()
        for _pop in selection.Population.unique():
            selec = select_a_given_pop(config.genotypes_two_lines, _pop)
            df = pd.DataFrame(
                {
                    "Population": _pop,
                    "ProbabilitySameInd": compute_proba_same_ind(config, selec, _pop),
                },
                index=[0],
            )
            df_probas = pd.concat([df, df_probas], ignore_index=True)
        selection_name += "_pops"

    elif config.agg_type == "pop_years":
        df_probas = pd.DataFrame()
        selection = selection[["Population", "Year"]].drop_duplicates()
        for _pop, _year in zip(selection.Population, selection.Year):
            selec = select_a_given_pop_year(config.genotypes_two_lines, _pop, _year)
            df = pd.DataFrame(
                {
                    "Population": _pop,
                    "Year": _year,
                    "ProbabilitySameInd": compute_proba_same_ind(config, selec, f"{_pop}_{_year}"),
                },
                index=[0],
            )
            df_probas = pd.concat([df, df_probas], ignore_index=True)
        selection_name += "_pops_years"

    else:  # subcat
        df_probas = pd.DataFrame()

        for _pop, _year, _sub in zip(
            selection.Population, selection.Year, selection.Subcategory
        ):
            selec = select_a_given_pop_year_subcat(
                config.genotypes_two_lines, _pop, _year, _sub
            )
            df = pd.DataFrame(
                {
                    "Population": _pop,
                    "Year": _year,
                    "Subcategory": _sub,
                    "ProbabilitySameInd": compute_proba_same_ind(config, selec, f"{_pop}_{_year}_{_sub}"),
                },
                index=[0],
            )
            df_probas = pd.concat([df, df_probas], ignore_index=True)
        selection_name += "_subcat"

    # save
    df_probas.to_csv(
        f"{config.output_path}/heterozygosity/probability_same_inds_{selection_name}.csv",
        index=False,
        sep=";",
    )


def calculate_frequencies(
    config: FileConfiguration,
    selection_genotypes_two_lines: pd.DataFrame,
    selection_name: str,
) -> pd.DataFrame:

    df = pd.DataFrame()
    for loci in config.loci_list:
        df2 = pd.DataFrame()
        df2["frequency"] = selection_genotypes_two_lines[loci].value_counts().dropna()
        df2["frequency"] = df2["frequency"] / df2["frequency"].sum()
        df2["locus"] = loci
        df2 = df2.reset_index()
        df2 = df2.rename(columns={loci: "allele"})
        df = pd.concat([df, df2], ignore_index=True)

    # Filter out rows with frequency 1 (fixed genotypes)
    # df = df[df["frequency"] != 1].copy()

    # Assure format 3 digit
    df["allele"] = df["allele"].apply(lambda x: f"{x:03d}")

    # Save
    df.to_csv(
        f"{config.output_path}/raw_data/frequencies_{selection_name}.csv",
        index=False,
        sep=";",
    )
    return df


def heterozygosity_calculation(
    config: FileConfiguration,
    df_frequencies: pd.DataFrame,
    selection_genotypes_two_lines: pd.DataFrame,
    selection_name: str,
) -> Tuple[pd.DataFrame, str]:

    he_distr = pd.DataFrame()
    for _locus in df_frequencies["locus"].unique():
        locus_data = df_frequencies[df_frequencies["locus"] == _locus][
            "frequency"
        ].values
        he_loc_sums = np.sum([np.power(i, 2) for i in locus_data])
        he_distr = pd.concat(
            [
                he_distr,
                pd.DataFrame({"locus": _locus, "he": 1 - he_loc_sums}, index=[0]),
            ]
        )

    he_distr.to_csv(
        f"{config.output_path}/heterozygosity/TableHeterozygosity_{selection_name}.csv",
        sep=";",
        index=False,
    )
    he_global = he_distr.he.sum() / (
        len(selection_genotypes_two_lines.columns) - 4
    )  # Remove metadata columns
    return he_distr, he_global


def plot_heterozygosity(
    config: FileConfiguration, he_distr: pd.DataFrame, selection_name: str
) -> None:
    he_distr = he_distr.sort_values(by="he", ascending=False)
    plt.figure(figsize=(8, 6))  # Adjust figure size as needed
    plt.bar(he_distr.locus, he_distr.he, color="skyblue", ec="k")
    plt.xlabel("Locus")
    plt.ylabel("Expected Heterozygosity (H0)")
    plt.title(f"Expected Heterozygosity for each Locus in {selection_name}")
    plt.xticks(rotation=45, ha="right")  # Rotate locus labels for better readability
    plt.ylim(0, 1)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(
        f"{config.output_path}/heterozygosity/heterozygosity_{selection_name}.pdf"
    )
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
        config, df_frequencies, selection, selection_name
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
        df_frequencies = calculate_frequencies(config, selec, config.selection_name + f"_{_pop}")
        he_distr, he_value = heterozygosity_calculation(
            config, df_frequencies, selec, selection_name
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
        df_frequencies = calculate_frequencies(config, selec, config.selection_name +f"_{_pop}_{_year}")
        he_distr, he_value = heterozygosity_calculation(
            config, df_frequencies, selec, f"{_pop}_{_year}"
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
        df_frequencies = calculate_frequencies(config, selec, config.selection_name + f"_{_pop}_{_year}{_ext}")
        he_distr, he_value = heterozygosity_calculation(
            config, df_frequencies, selec, f"{_pop}_{_year}{_ext}"
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
            f"{config.output_path}/raw_data/frequencies_{selection_name}.csv", sep=";"
        )
        df_pops["population"] = "all"
        df_pops["legend"] = selection_name

    elif agg_type == "pops":
        for _pop in selection.Population.unique():
            df = pd.read_csv(
                f"{config.output_path}/raw_data/frequencies_{selection_name}_{_pop}.csv", sep=";"
            )
            df["population"] = _pop
            df["legend"] = _pop
            df_pops = pd.concat([df_pops, df], ignore_index=True)

    elif agg_type == "pop_years":
        selection = selection[["Population", "Year"]].drop_duplicates()
        for _pop, _year in zip(selection.Population, selection.Year):
            df = pd.read_csv(
                f"{config.output_path}/raw_data/frequencies_{selection_name}_{_pop}_{_year}.csv", sep=";"
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
                f"{config.output_path}/raw_data/frequencies_{selection_name}_{_pop}_{_year}{_ext}.csv",
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
    nbr_rows_plots = int(np.ceil(nb_plots) / 3) if int(np.ceil(nb_plots) / 3) > 1 else 2
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
    if (nb_plots % 3 != 0) | (int(np.ceil(nb_plots / 3)) == 1):
        while plot_index_row < 3:
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
