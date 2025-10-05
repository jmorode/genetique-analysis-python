import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from ..core.config import FileConfiguration
from ..utils.utils import get_extension_if_subcat
from .generation import generate_n_families_following_frequencies_for_pairwise_distances


def compute_ks_test_between_pairwise_distance_and_relationships_simu_per_selection(
    config: FileConfiguration, nb_families: int, selection_name: str
) -> pd.DataFrame:
    # Recover sample values from file and compute frequency distribution
    data_ech = pd.read_csv(
        f"{config.output_path}/pairwise_differences/pairwise_differences_{selection_name}.csv",
        sep=";",
    )
    data_ech = data_ech.set_index("Sample")
    vals = np.array(data_ech)

    # take only lower triangular without diagonal
    vals = vals[np.tril_indices(len(data_ech), k=-1)]
    df_sample = pd.DataFrame({"Sample": vals})
    df_sample["Sample"] = df_sample["Sample"].astype("float")

    # Get relationship simu values
    data_sim = generate_n_families_following_frequencies_for_pairwise_distances(
        config, nb_families, selection_name
    )

    # Concat both
    if len(data_sim) > len(df_sample):
        data_ks = data_sim.copy()
        data_ks["Sample"] = df_sample["Sample"].copy()
    else:
        data_ks = df_sample.copy()
        for _rel in ["U", "PO", "FS", "HS"]:
            data_ks[_rel] = data_sim[_rel].copy()

    # Kolmogorov Smirnov
    df_ks = pd.DataFrame(
        columns=data_ks.columns.tolist(), index=data_ks.columns.tolist()
    )
    for col in data_ks.columns:
        for ind in data_ks.columns:
            df_ks.loc[ind][col] = ks_2samp(
                data_ks[ind].dropna(), data_ks[col].dropna()
            ).pvalue
    df_ks.to_csv(
        f"{config.output_path}/test_stats/table_ks_test_pairwise_vs_relationships_simu_{selection_name}.csv",
        sep=";",
        index=True,
    )
    return data_ks


def get_mean_std_for_ks_test_between_pairwise_distance_and_relationships_simu(
    raw_data: pd.DataFrame, name_sample: str
) -> pd.DataFrame:
    stats = raw_data.describe().reset_index()
    stats = stats[stats["index"].isin(["mean", "std"])]
    piv_stats = stats.pivot(columns="index", values=["Sample", "U", "PO", "FS", "HS"])
    piv_stats["name"] = name_sample
    return piv_stats


def compute_ks_test_between_pairwise_distance_and_relationships_simu(
    config: FileConfiguration, nb_families: int
) -> None:
    df_stats = pd.DataFrame()

    if config.agg_type == "all":
        data_ks = compute_ks_test_between_pairwise_distance_and_relationships_simu_per_selection(
            config, nb_families, config.selection_name + "_all"
        )
        df_stats = (
            get_mean_std_for_ks_test_between_pairwise_distance_and_relationships_simu(
                data_ks, config.selection_name + "_all"
            )
        )

    elif config.agg_type == "pops":
        for _pop in config.pops_to_select.Population.unique():
            data_ks = compute_ks_test_between_pairwise_distance_and_relationships_simu_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}"
            )
            sample_stat = get_mean_std_for_ks_test_between_pairwise_distance_and_relationships_simu(
                data_ks, config.selection_name + f"_{_pop}"
            )
            df_stats = pd.concat([df_stats, sample_stat], ignore_index=True)

    elif config.agg_type == "pop_years":
        for _pop, _year in zip(
            config.pops_to_select.Population, config.pops_to_select.Year
        ):
            data_ks = compute_ks_test_between_pairwise_distance_and_relationships_simu_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}_{_year}"
            )
            sample_stat = get_mean_std_for_ks_test_between_pairwise_distance_and_relationships_simu(
                data_ks, config.selection_name + f"_{_pop}_{_year}"
            )
            df_stats = pd.concat([df_stats, sample_stat], ignore_index=True)

    else:  # subcat
        for _pop, _year, _sub in zip(
            config.pops_to_select.Population,
            config.pops_to_select.Year,
            config.pops_to_select.Subcategory,
        ):
            _ext = get_extension_if_subcat(_sub)
            data_ks = compute_ks_test_between_pairwise_distance_and_relationships_simu_per_selection(
                config, nb_families, config.selection_name + f"_{_pop}_{_year}{_ext}"
            )
            sample_stat = get_mean_std_for_ks_test_between_pairwise_distance_and_relationships_simu(
                data_ks, config.selection_name + f"_{_pop}_{_year}{_ext}"
            )
            df_stats = pd.concat([df_stats, sample_stat], ignore_index=True)

    df_stats = df_stats.groupby("name").max()
    df_stats.to_csv(
        f"{config.output_path}/test_stats/stats_data_pairwise_sample_and_relationships_simu_{config.selection_name + '_' + config.agg_type}.csv",
        sep=";",
        index=True,
    )
