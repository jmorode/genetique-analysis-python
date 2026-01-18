import os
from typing import Tuple

import numpy as np
import pandas as pd

from ..core.config import FileConfiguration
from ..core.constants import NB_IND_RELATION, REFERENCE_ARRAY
from ..utils.utils import (
    conversion_two_lines_to_one_lines_genotypes,
    get_extension_if_subcat,
)
from .generation import generate_n_families_following_frequencies_for_ml_relate


def write_ml_relate_input_file(
    config: FileConfiguration, df: pd.DataFrame, legend: str
) -> None:
    out_filename = f"{config.output_path}/ml_relate/input_ml_relate_{legend}.txt"
    if os.path.exists(out_filename):
        os.remove(out_filename)
    # Incipit
    with open(out_filename, "a") as file:
        file.write(f"Microsatellite genotypes - {legend}\n")
        file.write("\n".join(config.loci_list) + "\n")
        file.write("POP\n")
    # Genotypes
    df.to_csv(out_filename, index=False, sep=" ", header=False, mode="a")


def reformat_genotypes_selection_for_ml_relate(
    config: FileConfiguration,
) -> pd.DataFrame:
    df_genotypes = config.genotypes_two_lines.copy()
    df_genotypes = conversion_two_lines_to_one_lines_genotypes(
        df_genotypes, ["Sample", "Population", "Year", "Subcategory"], config.loci_list
    )

    for locus in config.loci_list:
        # Format alleles as 3-digit numbers with leading zeros (e.g., 5 -> "005", 42 -> "042")
        df_genotypes[locus] = (
            df_genotypes[f"{locus}_1"].apply(lambda x: f"{int(x):03d}" if pd.notna(x) else "000")
            + df_genotypes[f"{locus}_2"].apply(lambda x: f"{int(x):03d}" if pd.notna(x) else "000")
        )
    df_genotypes = df_genotypes[["Sample"] + list(config.loci_list)]
    df_genotypes["Sample"] = df_genotypes.Sample.apply(lambda x: x + ",")
    return df_genotypes


def write_ml_relate_input_file_genotypes(config: FileConfiguration) -> None:
    df = reformat_genotypes_selection_for_ml_relate(config)
    write_ml_relate_input_file(config, df, config.selection_name)


def write_ml_relate_input_file_simu_families_per_selection(
    config: FileConfiguration, nb_families: int, selection_name: str
) -> None:
    df = generate_n_families_following_frequencies_for_ml_relate(
        config, nb_families, selection_name
    )
    write_ml_relate_input_file(config, df, "simulation_families_" + selection_name)


def write_ml_relate_input_file_simu_families(
    config: FileConfiguration, nb_families: int
) -> None:
    if config.agg_type == "all":
        write_ml_relate_input_file_simu_families_per_selection(
            config, nb_families, config.selection_name + "_all"
        )

    elif config.agg_type == "pops":
        for _pop in config.pops_to_select.Population.unique():
            write_ml_relate_input_file_simu_families_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}"
            )

    elif config.agg_type == "pop_years":
        for _pop, _year in zip(
            config.pops_to_select.Population, config.pops_to_select.Year
        ):
            write_ml_relate_input_file_simu_families_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}_{_year}"
            )

    else:  # subcat
        for _pop, _year, _sub in zip(
            config.pops_to_select.Population,
            config.pops_to_select.Year,
            config.pops_to_select.Subcategory,
        ):
            _ext = get_extension_if_subcat(_sub)
            write_ml_relate_input_file_simu_families_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}_{_year}{_ext}"
            )


def count_nbr_relationships_from_ml_relate_output(
    config: FileConfiguration,
) -> pd.DataFrame:
    df = pd.read_csv(
        f"{config.input_path}/ml_relate_output_{config.selection_name}.csv", sep=";"
    )
    po_count, fs_count, hs_count, un_count = 0, 0, 0, 0
    for col in df.columns:
        relationship_counts = df[col].value_counts()
        po_count += relationship_counts.get("PO", 0)
        fs_count += relationship_counts.get("FS", 0)
        hs_count += relationship_counts.get("HS", 0)
        un_count += relationship_counts.get("U", 0)

    df_counts = pd.DataFrame(
        {"PO": po_count, "FS": fs_count, "HS": hs_count, "U": un_count}, index=[0]
    )
    df_counts.to_csv(f"{config.output_path}/ml_relate/nbr_relationships_{config.selection_name}.csv",
                     sep = ";",
                     index=False)
    return df_counts


def compare_numpy_arrays(
    array1: np.array, array_ref: np.array
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    assert array1.shape == array_ref.shape, "Arrays must have the same shape"

    match_mask = array1 == array_ref

    # get matching relationships
    matching_indices = np.where(match_mask)
    matching = list(zip(array1[matching_indices], array_ref[matching_indices]))
    matching_counts = (
        pd.Series([list(i)[0] for i in matching]).value_counts().reset_index()
    )
    matching_counts = matching_counts.rename(columns={"index": "true_relation"})
    matching_counts = matching_counts[matching_counts.true_relation != "-"]

    # get mismatching relationships
    mismatch_indices = np.where(~match_mask)
    mismatches_counts = (
        pd.DataFrame(
            {
                "found_relation": array1[mismatch_indices],
                "true_relation": array_ref[mismatch_indices],
            }
        )
        .value_counts()
        .reset_index()
    )

    return matching_counts, mismatches_counts


def add_missing_relationships_matches(df: pd.DataFrame) -> pd.DataFrame:
    for _rel in ["PO", "FS", "HS", "U"]:
        if _rel not in df["true_relation"].to_list():
            df = pd.concat(
                [df, pd.DataFrame({"true_relation": _rel, "count": 0}, index=[0])]
            )
    return df


def add_missing_relationships_mismatches(df: pd.DataFrame) -> pd.DataFrame:
    df["temp_duet"] = df.apply(lambda x: x.found_relation + x.true_relation, axis=1)
    for _rel_f in ["PO", "FS", "HS", "U"]:
        for _rel_t in ["PO", "FS", "HS", "U"]:
            if _rel_t == _rel_f:
                continue
            if _rel_f + _rel_t not in df["temp_duet"].to_list():
                df = pd.concat(
                    [
                        df,
                        pd.DataFrame(
                            {
                                "found_relation": _rel_f,
                                "true_relation": _rel_t,
                                "count": 0,
                            },
                            index=[0],
                        ),
                    ]
                )
    df = df.drop(columns="temp_duet")
    return df.sort_values(by=["found_relation", "true_relation"])


def get_percentage_and_count_relationships_families(
    df_stats: pd.DataFrame, meta_col: list[str], nb_families: int
) -> pd.DataFrame:
    df_stats = df_stats.groupby(meta_col)["count"].sum().reset_index()
    df_stats = add_missing_relationships_matches(df_stats)
    df_stats["percentage"] = df_stats.apply(
        lambda x: x["count"] / (NB_IND_RELATION[x[meta_col[0]]] * nb_families) * 100,
        axis=1,
    )
    return df_stats


def get_total_percentage_over_families(
    df_stats: pd.DataFrame, nb_families: int
) -> float:
    return (
        df_stats["count"].sum()
        / (np.sum(list(NB_IND_RELATION.values())) * nb_families)
        * 100
    )


def get_sum_percentage_by_category(df_stats: pd.DataFrame, category: str) -> float:
    return df_stats[df_stats.found_relation == category]["percentage"].sum()


def overall_stats_relationships_simulated_families(
    df_matching_stats: pd.DataFrame,
    df_mismatching_stats: pd.DataFrame,
    nb_families: int,
) -> pd.DataFrame:
    # Stats over all families
    df_matching_stats = get_percentage_and_count_relationships_families(
        df_matching_stats, ["true_relation"], nb_families
    )
    df_mismatching_stats = get_percentage_and_count_relationships_families(
        df_mismatching_stats, ["found_relation", "true_relation"], nb_families
    )

    # Format output file
    # Global stats
    df_reliability = pd.DataFrame(
        {
            "category": [
                "Total percentage of correct relationships",
                "Total percentage of incorrect relationships",
                "Average percentage of correct relationships",
                "Percentage of incorrect PO",
                "Percentage of incorrect FS",
                "Percentage of incorrect HS",
                "Percentage of incorrect U",
            ],
            "percentage": [
                get_total_percentage_over_families(df_matching_stats, nb_families),
                get_total_percentage_over_families(df_mismatching_stats, nb_families),
                df_matching_stats.percentage.mean(),
                get_sum_percentage_by_category(df_mismatching_stats, "PO"),
                get_sum_percentage_by_category(df_mismatching_stats, "FS"),
                get_sum_percentage_by_category(df_mismatching_stats, "HS"),
                get_sum_percentage_by_category(df_mismatching_stats, "U"),
            ],
        }
    )
    # Correct by relationships
    df_matching_stats["true_relation"] = df_matching_stats["true_relation"].apply(
        lambda x: "Percentage of correct " + x
    )
    df_matching_stats = df_matching_stats.rename(columns={"true_relation": "category"})

    # Incorrect by relationships
    df_mismatching_stats["category"] = df_mismatching_stats.apply(
        lambda x: "Percentage of incorrect "
        + x.found_relation
        + " that should be "
        + x.true_relation,
        axis=1,
    )
    df_global_stats = pd.concat(
        [
            df_reliability,
            df_matching_stats,
            df_mismatching_stats.drop(columns=["true_relation", "found_relation"]),
        ],
        ignore_index=True,
    )

    return df_global_stats


def reliability_ml_relate_based_simulated_families(
    config: FileConfiguration, nb_families: int
) -> pd.DataFrame:
    df = pd.read_csv(
        f"{config.input_path}/ml_relate_output_simulated_families_{config.selection_name}.csv",
        sep=";",
    )
    df_matching_stats = pd.DataFrame()
    df_mismatching_stats = pd.DataFrame()

    for i in range(nb_families):
        # Select family
        fam = f"{i:04d}"
        # ML-Relate output is already ordered correctly; only extract the family submatrix
        family_cols = [
            c
            for c in df.columns
            if isinstance(c, str) and c != "Unnamed: 0" and c.endswith(f"_{fam}")
        ]
        if len(family_cols) != REFERENCE_ARRAY.shape[0]:
            raise ValueError(
                f"Expected {REFERENCE_ARRAY.shape[0]} individuals for family {fam}, "
                f"found {len(family_cols)} columns: {family_cols}")

        df_family = df.loc[df["Unnamed: 0"].isin(family_cols), ["Unnamed: 0"] + family_cols].copy()
        df_family = df_family.set_index("Unnamed: 0")

        # Standardize labels to match REFERENCE_ARRAY (e.g., "DF_0000" -> "DF")
        df_family.index = df_family.index.map(
            lambda x: x.split("_")[0] if isinstance(x, str) else x
        )
        df_family.columns = [col.split("_")[0] for col in df_family.columns]

        # Compute stats
        matching_counts, mismatching_counts = compare_numpy_arrays(
            df_family.to_numpy(), REFERENCE_ARRAY
        )

        # Append stats of the family
        df_matching_stats = pd.concat(
            [df_matching_stats, matching_counts], ignore_index=True
        )
        df_mismatching_stats = pd.concat(
            [df_mismatching_stats, mismatching_counts], ignore_index=True
        )

    df_global_stats = overall_stats_relationships_simulated_families(
        df_matching_stats, df_mismatching_stats, nb_families
    )
    df_global_stats.to_csv(
        f"{config.output_path}/ml_relate/reliability_stats_ml_relate_{config.selection_name}.csv",
        sep=";",
        index=False,
    )
    return df_global_stats
