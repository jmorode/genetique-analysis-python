import os
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from ..core.config import FileConfiguration
from ..core.constants import rename_parents
from ..utils.utils import (
    conversion_two_lines_to_one_lines_genotypes,
    get_extension_if_subcat,
    set_random_seed,
)

# Set random seed for reproducibility
set_random_seed(12)

from .pairwise_differences import calculate_pairwise_differences


def generate_unrelated_individuals_following_frequencies(
    config: FileConfiguration, frequencies: pd.DataFrame, num_individuals: int
) -> pd.DataFrame:
    inds_list = [f"Ind{i}" for i in range(1, num_individuals + 1)]
    ind = pd.DataFrame(index=np.sort(inds_list + inds_list), columns=config.loci_list)

    for locus in config.loci_list:
        locus_data = frequencies[frequencies["locus"] == locus]
        allele_frequencies = list(locus_data["frequency"])
        alleles = list(locus_data["allele"])

        ind[locus] = random.choices(
            alleles, weights=allele_frequencies, k=2 * num_individuals
        )
    ind = ind.reset_index().rename(columns={"index": "POP"})
    return ind


def get_child_allele_at_locus(genotype: pd.DataFrame, ind1: str, ind2: str, locus: str):
    # Get random allele from ind1
    ind1_locus = genotype[genotype.POP == ind1][locus].to_list()
    allele1 = random.choice(ind1_locus)

    # Get random allele from ind2
    ind2_locus = genotype[genotype.POP == ind2][locus].to_list()
    allele2 = random.choice(ind2_locus)
    return [allele1, allele2]


def generate_child(
    config: FileConfiguration,
    parents_genotypes: pd.DataFrame,
    ind1: str,
    ind2: str,
    name: str,
) -> pd.DataFrame:
    df_child = pd.DataFrame(index=[name, name], columns=config.loci_list)

    for _locus in config.loci_list:
        list_allele = get_child_allele_at_locus(parents_genotypes, ind1, ind2, _locus)
        df_child[_locus] = list_allele

    return df_child.reset_index(names=["POP"])


def generate_related_individuals_following_frequencies(
    config: FileConfiguration, selection_name: str
) -> pd.DataFrame:
    frequencies = pd.read_csv(
        f"{config.output_path}/raw_data/frequencies_{config.project_name}_{selection_name}.csv",
        sep=";",
    )
    # Generate parents : P1 (ind1), M (ind2), P2 (ind3)
    df_parents = generate_unrelated_individuals_following_frequencies(
        config, frequencies, num_individuals=3
    )
    df_parents["POP"] = df_parents.POP.apply(lambda x: rename_parents[x])

    # Generate progeny : F1, F2, DF
    df_f1 = generate_child(config, df_parents, "P1", "M", "F1")
    df_f2 = generate_child(config, df_parents, "P1", "M", "F2")
    df_df = generate_child(config, df_parents, "M", "P2", "DF")

    df_family = pd.concat([df_parents, df_f1, df_f2, df_df], ignore_index=True)
    return df_family


def generate_n_families_following_frequencies_for_ml_relate(
    config: FileConfiguration, nb_families: int, selection_name: str
) -> pd.DataFrame:
    df_families = pd.DataFrame()
    for i in range(nb_families):
        family = generate_related_individuals_following_frequencies(
            config, selection_name
        )
        # Zero-pad family number to 4 digits (e.g., 0 -> 0000, 1 -> 0001)
        family["POP"] = family.POP.apply(lambda x: x + f"_{i:04d},")
        df_families = pd.concat([df_families, family], ignore_index=True)

    # reformat for ml_relate
    df_families = conversion_two_lines_to_one_lines_genotypes(
        df_families, ["POP"], config.loci_list
    )
    for locus in config.loci_list:
        # Format alleles as 3-digit numbers with leading zeros (e.g., 5 -> "005", 42 -> "042")
        df_families[locus] = (
            df_families[f"{locus}_1"].apply(lambda x: f"{int(x):03d}" if pd.notna(x) else "000")
            + df_families[f"{locus}_2"].apply(lambda x: f"{int(x):03d}" if pd.notna(x) else "000")
        )
    df_families = df_families[["POP"] + list(config.loci_list)]
    
    # Sort to keep families together and in numeric order
    # Extract family number and individual type for sorting
    # POP format: "DF_0000," -> extract "0000" and "DF"
    df_families["_family_num"] = df_families["POP"].apply(
        lambda x: int(x.rstrip(',').split('_')[1]) if len(x.rstrip(',').split('_')) == 2 else 999999
    )
    df_families["_ind_type"] = df_families["POP"].apply(
        lambda x: x.rstrip(',').split('_')[0] if len(x.rstrip(',').split('_')) == 2 else "ZZZ"
    )
    # Define order within family: DF, M, F1, F2, P1, P2
    ind_order = {"DF": 0, "M": 1, "F1": 2, "F2": 3, "P1": 4, "P2": 5}
    df_families["_ind_order"] = df_families["_ind_type"].apply(lambda x: ind_order.get(x, 99))
    
    # Sort by family number first, then by individual order within family
    df_families = df_families.sort_values(by=["_family_num", "_ind_order"]).reset_index(drop=True)
    
    # Remove temporary sorting columns
    df_families = df_families.drop(columns=["_family_num", "_ind_type", "_ind_order"])
    
    return df_families


def _process_single_family_for_pairwise_distances(
    config: FileConfiguration, frequencies: pd.DataFrame, family_id: int
) -> pd.DataFrame:
    """
    Process a single family to generate pairwise differences.

    This is a helper function designed to be used in parallel processing.
    It generates one family and extracts the specific pairwise differences
    needed for relationship analysis.

    Args:
        config: FileConfiguration object containing project settings
        frequencies: Pre-loaded frequency data to avoid repeated file I/O
        family_id: Unique identifier for the family (used for random seed)

    Returns:
        DataFrame with pairwise differences for one family
    """
    # Set unique random seed for each family to ensure reproducibility
    # while maintaining independence between families
    set_random_seed(12 + family_id)

    # Generate parents : P1 (ind1), M (ind2), P2 (ind3)
    df_parents = generate_unrelated_individuals_following_frequencies(
        config, frequencies, num_individuals=3
    )
    df_parents["POP"] = df_parents.POP.apply(lambda x: rename_parents[x])

    # Generate progeny : F1, F2, DF
    df_f1 = generate_child(config, df_parents, "P1", "M", "F1")
    df_f2 = generate_child(config, df_parents, "P1", "M", "F2")
    df_df = generate_child(config, df_parents, "M", "P2", "DF")

    family = pd.concat([df_parents, df_f1, df_f2, df_df], ignore_index=True)
    family = conversion_two_lines_to_one_lines_genotypes(
        family, ["POP"], config.loci_list
    )

    # Add missing metadata columns to match expected format for calculate_pairwise_differences
    # Expected format: Sample, Population, Year, Subcategory, then loci
    family_formatted = family.copy()
    family_formatted.insert(1, "Population", family_formatted["POP"])
    family_formatted.insert(2, "Year", "")
    family_formatted.insert(3, "Subcategory", "")
    # Rename POP to Sample to match expected column name
    family_formatted = family_formatted.rename(columns={"POP": "Sample"})

    # Temporarily add generated family labels to dict_pop_samples if it exists
    # to avoid KeyError when mapping in calculate_pairwise_differences
    original_dict = None
    dict_was_none = False
    if hasattr(config, 'dict_pop_samples'):
        if config.dict_pop_samples is not None:
            original_dict = config.dict_pop_samples.copy()
            for sample in family_formatted["Sample"].unique():
                if sample not in config.dict_pop_samples:
                    config.dict_pop_samples[sample] = "Generated"
        else:
            # Initialize dict_pop_samples if it's None
            dict_was_none = True
            config.dict_pop_samples = {}
            for sample in family_formatted["Sample"].unique():
                config.dict_pop_samples[sample] = "Generated"
    else:
        # Create dict_pop_samples if it doesn't exist
        config.dict_pop_samples = {}
        for sample in family_formatted["Sample"].unique():
            config.dict_pop_samples[sample] = "Generated"
    
    try:
        diff_family = calculate_pairwise_differences(
            config, config.selection_name + f"_temp", family_formatted
        )

    finally:
        # Restore original dict_pop_samples if we modified it
        if original_dict is not None:
            config.dict_pop_samples = original_dict
        elif dict_was_none:
            config.dict_pop_samples = None

    return pd.DataFrame(
        {
            "U": diff_family.loc["P1", "P2"],
            "PO": diff_family.loc["P1", "F1"],
            "FS": diff_family.loc["F1", "F2"],
            "HS": diff_family.loc["F1", "DF"],
        },
        index=[0],
    )


def generate_n_families_following_frequencies_for_pairwise_distances(
    config: FileConfiguration,
    nb_families: int,
    selection_name: str,
    max_workers: int = None,
    use_parallel: bool = True,
) -> pd.DataFrame:
    """
    Generate multiple families and calculate pairwise distances using parallel processing.

    This optimized version uses parallel processing to generate families concurrently,
    significantly improving performance for large numbers of families. It also loads
    the frequency data once to avoid repeated file I/O operations.

    Args:
        config: FileConfiguration object containing project settings
        nb_families: Number of families to generate
        selection_name: Name of the selection for frequency data
        max_workers: Maximum number of parallel workers (default: CPU count)
        use_parallel: Whether to use parallel processing (fallback to sequential if False)

    Returns:
        DataFrame with pairwise differences for all families
    """
    # Load frequency data once to avoid repeated file I/O
    frequencies = pd.read_csv(
        f"{config.output_path}/raw_data/frequencies_{config.project_name}_{selection_name}.csv",
        sep=";",
    )

    if not use_parallel or nb_families < 2:
        # Fallback to sequential processing for small numbers or when parallel is disabled
        return _generate_families_sequential(config, frequencies, nb_families)

    if max_workers is None:
        max_workers = min(os.cpu_count() or 1, nb_families)

    try:
        # Use parallel processing for family generation
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Create partial function with fixed config and frequencies
            process_func = partial(
                _process_single_family_for_pairwise_distances, config, frequencies
            )

            # Submit all family processing tasks
            future_to_family = {
                executor.submit(process_func, family_id): family_id
                for family_id in range(nb_families)
            }

            # Collect results as they complete
            results = []
            for future in as_completed(future_to_family):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as exc:
                    family_id = future_to_family[future]
                    print(f"Family {family_id} generated an exception: {exc}")
                    raise

        # Concatenate all results at once (more efficient than loop concatenation)
        if results:
            df_pairwise = pd.concat(results, ignore_index=True)
        else:
            df_pairwise = pd.DataFrame(columns=["U", "PO", "FS", "HS"])

        return df_pairwise

    except Exception as e:
        print(f"Parallel processing failed: {e}")
        print("Falling back to sequential processing...")
        return _generate_families_sequential(config, frequencies, nb_families)


def _generate_families_sequential(
    config: FileConfiguration, frequencies: pd.DataFrame, nb_families: int
) -> pd.DataFrame:
    """
    Sequential fallback for family generation when parallel processing fails.

    Args:
        config: FileConfiguration object containing project settings
        frequencies: Pre-loaded frequency data
        nb_families: Number of families to generate

    Returns:
        DataFrame with pairwise differences for all families
    """
    results = []
    for family_id in range(nb_families):
        result = _process_single_family_for_pairwise_distances(
            config, frequencies, family_id
        )
        results.append(result)

    if results:
        return pd.concat(results, ignore_index=True)
    else:
        return pd.DataFrame(columns=["U", "PO", "FS", "HS"])



def get_j_i_position_for_two_by_two_plots(i: int) -> Tuple[int, int]:
    if i < 2:
        j = 0
    else:
        j = 1
        i -= 2
    return j, i


def plot_pairwise_distances_per_relationships_per_selection(
    config: FileConfiguration, nb_families: int, selection_name: str
) -> None:
    # get pairwise for each relationship
    df_pairwise = generate_n_families_following_frequencies_for_pairwise_distances(
        config, nb_families, selection_name
    )

    # Plot distributions
    fig, ax = plt.subplots(
        nrows=2, ncols=2, figsize=(14, 10)
    )  # Adjust figure size as needed


    df_freq_all = pd.DataFrame()
    for i, rel in enumerate(df_pairwise.columns):
        j, i = get_j_i_position_for_two_by_two_plots(i)
        df_freq = (
            df_pairwise[rel]
            .value_counts(normalize=True)
            .reset_index()
            .rename(columns={rel: "nb_diff"})
        )
        ax[j, i].bar(
            df_freq["nb_diff"], df_freq["proportion"], width=0.8, color="darkgrey"
        )
        ax[j, i].set_xlabel("Number of differences")
        ax[j, i].set_ylabel("Frequency")
        ax[j, i].set_title(f"{rel} Distribution")
        ax[j, i].set_xlim(
            0, 2 * len(config.loci_list) + 1
        )  # Adjust x-axis limits based on data
        ax[j, i].set_ylim(0, 0.25)  # Adjust y-axis limits based on data
        ax[j, i].axhline(
            y=0.05, color="blue", linestyle="--", label="Threshold"
        )  # Add threshold line

        # Save number of individuals simulated per relationship and number of differences
        df_freq_individuals = (
            df_pairwise[rel]
            .value_counts()
            .reset_index()
            .rename(columns={rel: "nb_diff"})
        )
        df_freq_individuals["Relationship"] = rel
        df_freq_all = pd.concat([df_freq_all, df_freq_individuals], ignore_index=True)

    plt.tight_layout()
    plt.savefig(
        f"{config.output_path}/pairwise_differences/plots/plot_pairwise_frequencies_by_relationships_{selection_name}.pdf"
    )
    plt.close()

    df_freq_all.sort_values(by=["Relationship", "nb_diff"]).to_csv(
        f"{config.output_path}/pairwise_differences/pairwise_frequencies_by_relationships_{selection_name}.csv",
        index=False,
        sep=";",
    )


def plot_pairwise_distances_per_relationships(
    config: FileConfiguration, nb_families: int
) -> None:
    if config.agg_type == "all":
        plot_pairwise_distances_per_relationships_per_selection(
            config, nb_families, config.selection_name + "_all"
        )

    elif config.agg_type == "pops":
        for _pop in config.pops_to_select.Population.unique():
            plot_pairwise_distances_per_relationships_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}"
            )

    elif config.agg_type == "pop_years":
        for _pop, _year in zip(
            config.pops_to_select.Population, config.pops_to_select.Year
        ):
            plot_pairwise_distances_per_relationships_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}_{_year}"
            )

    else:  # subcat
        for _pop, _year, _sub in zip(
            config.pops_to_select.Population,
            config.pops_to_select.Year,
            config.pops_to_select.Subcategory,
        ):
            _ext = get_extension_if_subcat(_sub)
            plot_pairwise_distances_per_relationships_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}_{_year}{_ext}"
            )


def plot_cumulated_and_not_frequencies_for_simu_relationships_and_sample_by_pop(
    config: FileConfiguration, nb_families: int
) -> None:
    assert config.agg_type == "pops"
    pdf = PdfPages(
        f"{config.output_path}/pairwise_differences/plots/freq_cum_pairwise_distances_relationships_{config.selection_name + '_' + config.agg_type}.pdf"
    )
    plot_index_row = 0

    for name_sample in list(config.pops_to_select.Population.unique()):

        if (plot_index_row % 3) == 0:
            if plot_index_row != 0:
                plt.tight_layout()
                pdf.savefig()
                plt.close()
            fig, ax = plt.subplots(3, 2, figsize=(30, 30))
            plot_index_row = 0

        # get pairwise diff
        raw_data = generate_n_families_following_frequencies_for_pairwise_distances(
            config,
            nb_families=nb_families,
            selection_name=config.selection_name + f"_{name_sample}",
        )

        # recover pairwise distrib compute frequency distribution
        data_ech = pd.read_csv(
            f"{config.output_path}/pairwise_differences/pairwise_differences_{config.selection_name + '_' + name_sample}.csv",
            sep=";",
        )
        data_ech = data_ech.set_index("Sample")
        vals = np.array(data_ech)
        vals = vals[np.tril_indices(len(data_ech), k=-1)]
        df_sample = pd.DataFrame({"Sample": vals})
        unique, counts = np.unique(vals, return_counts=True)
        add_row = pd.DataFrame(dict(zip(unique, counts / sum(counts))), index=[0])
        add_row["Sample"] = "Sample"
        add_row = add_row.set_index("Sample")
        add_row = add_row.T
        add_row = add_row.reset_index().dropna()
        add_row = add_row.rename(columns={"index": "nb_diff"})

        # Recover relationship simu values
        data_sim = pd.DataFrame(columns=["nb_diff"])
        for _rel in ["PO", "FS", "HS", "U"]:
            prop = (
                raw_data[_rel]
                .value_counts(normalize=True)
                .reset_index()
                .rename(columns={_rel: "nb_diff", "proportion": _rel})
            )
            data_sim = pd.merge(data_sim, prop, on="nb_diff", how="outer")

        # Concat sim & sample frequencies
        data = pd.merge(data_sim, add_row, on="nb_diff", how="outer")
        data = data.fillna(0).set_index("nb_diff")

        # get stat for mean value of distribution
        raw_data["Sample"] = df_sample["Sample"].astype("float").copy()

        stats = raw_data.describe().reset_index()

        # Plot
        for col, _c in zip(data.columns, ["blue", "orange", "green", "red", "purple"]):
            data[col].plot(
                ax=ax[plot_index_row, 0], legend=True, alpha=0.8, color=_c, fontsize=20
            )
            xx = np.array(data[col].index.values).astype("int")
            yy = np.array(data[col].values)
            ax[plot_index_row, 0].fill_between(x=xx, y1=yy, alpha=0.2, color=_c)
            ax[plot_index_row, 0].set_ylabel("Frequency", fontsize=20)
            ax[plot_index_row, 0].set_xlabel("Distances", fontsize=20)
            ax[plot_index_row, 0].set_xlim(0, len(config.loci_list) * 2)
            ax[plot_index_row, 0].set_ylim(0, 0.4)
            ax[plot_index_row, 0].vlines(
                x=stats[stats["index"] == "mean"][col].iloc[0],
                ymin=0,
                ymax=0.4,
                color=_c,
            )
            ax[plot_index_row, 0].set_title(name_sample, fontsize=20)
            ax[plot_index_row, 0].legend(loc="best", prop={"size": 20})

        # CumData
        # Cumsum of distrib
        cumdata = data.cumsum(axis=0)

        # Plot
        for col in cumdata.columns:
            cumdata[col].plot(ax=ax[plot_index_row, 1], legend=True, fontsize=20)
        ax[plot_index_row, 1].set_ylabel("Cumulative Frequency", fontsize=20)
        ax[plot_index_row, 1].set_xlabel("Distances", fontsize=20)
        ax[plot_index_row, 1].set_xlim(0, len(config.loci_list) * 2)
        ax[plot_index_row, 1].set_title(name_sample, fontsize=20)
        ax[plot_index_row, 1].legend(loc="lower right", prop={"size": 20})
        plot_index_row += 1

    # remove axis of empty plots
    while plot_index_row < 3:
        ax[plot_index_row, 0].set_axis_off()
        ax[plot_index_row, 1].set_axis_off()
        plot_index_row += 1

    pdf.savefig()
    plt.close()
    pdf.close()
