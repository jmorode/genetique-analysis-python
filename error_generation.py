import multiprocessing as mp
import random
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from config import FileConfiguration
from description import calculate_frequencies
from generation import generate_unrelated_individuals_following_frequencies
from joblib import Parallel, delayed
from pairwise_differences import calculate_pairwise_differences
from recaptures import get_recaptures_sample
from tqdm import tqdm
from utils import conversion_two_lines_to_one_lines_genotypes, select_a_given_pop


def get_another_allele_in_pop(
    locus: str, current_allele: float, df_freqs: pd.DataFrame
) -> float:
    alleles = df_freqs[df_freqs.locus == locus]["allele"].to_list()
    # do not catch the same allele as current - unless only one available
    if len(alleles) > 1:
        alleles = list(set(alleles) - {current_allele})
    return random.choice(alleles)


def introduce_errors(
    genotypes: pd.DataFrame,
    tx_error: float,
    df_frequency: pd.DataFrame,
    type_err: str,
    type_hezy_err: str,
) -> pd.DataFrame:
    individus = genotypes.POP.unique()
    loci = genotypes.drop(columns=["POP"]).columns
    nb_errors = len(individus) * len(loci) * tx_error

    if nb_errors != 0:
        # Get random ind / loci depending on type of err
        individus_loci_impacted = get_random_ind_loci(
            genotypes, individus, loci, type_err, nb_errors
        )

        # for each couple of loci/ind
        for _ind, _locus in individus_loci_impacted:
            # recover allele
            # get alleles at loci and its indexes
            ind_locus = genotypes[genotypes.POP == _ind][_locus].to_list()
            indexes = genotypes.index[genotypes.POP == _ind].tolist()

            # if loci homozygote - random allele of pop
            if len(np.unique(ind_locus)) == 1:
                current_alelle = float(list(np.unique(ind_locus))[0])
                new_allele = get_another_allele_in_pop(
                    _locus, current_alelle, df_frequency
                )
                # change it randomly on the loci
                genotypes.loc[random.choice(indexes), _locus] = new_allele

            # if heterozygote - loci become homozygote, 2 possibilities
            else:
                # largest allele lost
                if type_hezy_err == "hetero-biased":
                    index_to_change = indexes[np.argmax(ind_locus)]
                    genotypes.loc[index_to_change, _locus] = min(ind_locus)

                # random choice of allele to drop
                else:
                    index_to_change = random.choice(indexes)
                    index_to_keep = list(set(indexes) - {index_to_change})[0]
                    new_allele = genotypes.loc[index_to_keep, _locus]
                    genotypes.loc[index_to_change, _locus] = new_allele

    # genotypes.to_csv(f'inputs/genotypes_with_simu_errors_{pop_name}.csv', sep=";")
    return genotypes


def map_zygotes_loci(
    genotypes: pd.DataFrame, individuals: list[str], loci: list[str]
) -> Tuple[list[list[str]], list[list[str]], list[list[str]]]:
    homozygotes = []
    heterozygotes = []
    all_zygotes = []
    for _ind in individuals:
        for _locus in loci:
            ind_locus = genotypes[genotypes.POP == _ind][_locus].to_list()
            if len(np.unique(ind_locus)) == 1:
                homozygotes.append([_ind, _locus])
            else:
                heterozygotes.append([_ind, _locus])
            all_zygotes.append([_ind, _locus])
    return homozygotes, heterozygotes, all_zygotes


def get_random_choice_lists(
    zygotes: list[list[str]], nb_errors: float
) -> list[list[str]]:
    choice_index_zygotes = np.random.choice(
        np.arange(len(zygotes)), size=int(np.ceil(nb_errors)), replace=False
    )
    zygotes = [zygotes[int(i)] for i in choice_index_zygotes]
    return zygotes


def get_random_ind_loci(
    genotypes: pd.DataFrame,
    individus: list[str],
    loci: list[str],
    type_err: str,
    nb_errors: float,
) -> list[list[str]]:
    if type_err == "heterozygotes":
        _, heterozygotes, all_zygotes = map_zygotes_loci(genotypes, individus, loci)
        correction_heterozygotes = len(heterozygotes) / len(all_zygotes)
        individus_loci_impacted = get_random_choice_lists(
            heterozygotes, nb_errors * correction_heterozygotes
        )

    elif type_err == "homozygotes":
        homozygotes, _, all_zygotes = map_zygotes_loci(genotypes, individus, loci)
        correction_homozygotes = len(homozygotes) / len(all_zygotes)
        individus_loci_impacted = get_random_choice_lists(
            homozygotes, nb_errors * correction_homozygotes
        )

    else:  # type_err == "random"
        _, _, all_zygotes = map_zygotes_loci(genotypes, individus, loci)
        individus_loci_impacted = get_random_choice_lists(all_zygotes, nb_errors)

    return individus_loci_impacted


def get_nb_recaptures_following_error_introduction(
    config: FileConfiguration,
    pop_name: str,
    num_inds: int,
    error_rate: float,
    error_type: str,
    type_hezy_err: str,
) -> float:
    # 1/ Get N genotypes randomly following allelic frequencies of a given population
    df_genotype = select_a_given_pop(config.genotypes_two_lines, pop_name)
    df_frequency = calculate_frequencies(
        config, df_genotype, config.selection_name + f"_{pop_name}_error"
    )
    df_generated = generate_unrelated_individuals_following_frequencies(
        config, df_frequency, num_inds
    )

    # 2/ Introduce errors randomly
    genotype_with_errors = introduce_errors(
        df_generated, error_rate, df_frequency, error_type, type_hezy_err
    )

    # 3/ Compute recaptures
    df_2lgn_genotypes_err = conversion_two_lines_to_one_lines_genotypes(
        genotype_with_errors, ["POP"], config.loci_list
    )
    df_pairwise = calculate_pairwise_differences(
        config, f"_{pop_name}_error", df_2lgn_genotypes_err
    )
    recaptures = get_recaptures_sample(config, df_pairwise, "", False)

    return len(recaptures)


def get_distrib_recap_errors(
    config: FileConfiguration,
    pop_name: str,
    num_inds: int,
    error_rate: float,
    error_type: str,
    type_hezy_err: str,
    nb_evaluation_recaptures: int,
) -> list[float]:

    distrib_recap = []
    for _ in tqdm(range(nb_evaluation_recaptures)):
        nb_recap = get_nb_recaptures_following_error_introduction(
            config, pop_name, num_inds, error_rate, error_type, type_hezy_err
        )
        distrib_recap.append(nb_recap)

    return distrib_recap


def plot_distrib_recap_errors(
    config: FileConfiguration,
    distrib_recap: list[float],
    pop_name: str,
    error_rate: float,
    error_type: str,
    hezy_type: str,
    nb_simulations: int,
) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    labels, counts = np.unique(distrib_recap, return_counts=True)
    plt.bar(labels, counts / nb_simulations, align="center")
    plt.gca().set_xticks(labels)
    ax.set_xlabel("Nombre de recaptures")
    ax.set_ylabel("Fréquence")
    hezy_text = f" / {hezy_type}" if hezy_type != "" else ""
    ax.set_title(
        f"Population: {pop_name} - Error rate: {error_rate} - Type: {error_type}"
        + hezy_text
    )
    ax.set_xlim(0, 20)
    fig.tight_layout()
    plt.savefig(
        f"{config.output_path}/error_recapture_tests/distribution_recap_{pop_name}_error_rate_{str(error_rate).replace('.','-')}_{error_type}_type_{hezy_type}.pdf"
    )
    plt.close()
    fig.clear()


def get_average_recaptures(list_recaptures: list[float]) -> float:
    return np.mean(list_recaptures).astype(float)


def get_recap_at_q95(list_recaptures: list[float]) -> float:
    return np.quantile(list_recaptures, q=0.95)


def compute_avg_q95_simu(
    config: FileConfiguration,
    name_pop: str,
    nb_generated_individuals: int,
    nb_simu: int,
    err_rate: float,
    err_type: str,
):
    # same seed of each sub process
    np.random.seed(42)
    random.seed(64)

    results = pd.DataFrame()
    # Hezy type only use for heterozygotes allocation errors (so heterozygotes & random err type)
    if (err_type == "homozygotes") | (err_rate == 0):
        list_hezy = [""]
    else:
        list_hezy = ["hetero-random", "hetero-biased"]

    for _hezy_err_type in list_hezy:
        list_recap = get_distrib_recap_errors(
            config,
            name_pop,
            nb_generated_individuals,
            err_rate,
            err_type,
            _hezy_err_type,
            nb_simu,
        )
        print(
            name_pop,
            nb_generated_individuals,
            err_rate,
            err_type,
            _hezy_err_type,
            nb_simu,
        )

        plot_distrib_recap_errors(
            config, list_recap, name_pop, err_rate, err_type, _hezy_err_type, nb_simu
        )
        q95 = get_recap_at_q95(list_recap)
        avg = get_average_recaptures(list_recap)
        print(avg)

        res = pd.DataFrame(
            {
                "Error rate": err_rate,
                "Error type": err_type,
                "Heterozygote Change Type": _hezy_err_type,
                "Mean recapture": avg,
                "q95 recapture": q95,
            },
            index=[0],
        )
        results = pd.concat([results, res], ignore_index=True)

    return results


def test_error_impact_scenarios_on_given_pop(
    config: FileConfiguration,
    name_population: str,
    nb_generated_inds: int,
    nb_evaluation_recap: int,
) -> None:
    results_parallel = Parallel(n_jobs=mp.cpu_count())(
        delayed(compute_avg_q95_simu)(
            config,
            name_population,
            nb_generated_inds,
            nb_evaluation_recap,
            _err_rate,
            _err_type,
        )
        for _err_rate in [0, 0.01, 0.1]
        for _err_type in ["homozygotes", "heterozygotes", "random"]
    )
    results_parallel = pd.concat(results_parallel)
    print(results_parallel)

    # Save
    mean_vals = pd.pivot(
        results_parallel,
        index="Error rate",
        columns=["Error type", "Heterozygote Change Type"],
        values="Mean recapture",
    )
    mean_vals.to_csv(
        f"{config.output_path}/error_recapture_tests/simulation_errors_{name_population}_average_recaptures.csv",
        sep=";",
        index=True,
    )

    q95_vals = pd.pivot(
        results_parallel,
        index="Error rate",
        columns=["Error type", "Heterozygote Change Type"],
        values="q95 recapture",
    )
    q95_vals.to_csv(
        f"{config.output_path}/error_recapture_tests/simulation_errors_{name_population}_q95_recaptures.csv",
        sep=";",
        index=True,
    )
