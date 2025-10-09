"""
Utility functions for genetics analysis.

This module contains helper functions for file operations, data selection,
and genotype format conversions used throughout the genetics analysis pipeline.
"""

import os
import random
from typing import List, Optional

import numpy as np
import pandas as pd


def set_random_seed(seed: int = 12) -> None:
    """
    Set random seed for reproducibility across all random number generators.

    This function ensures that both Python's random module and NumPy's random
    number generator use the same seed for consistent, reproducible results.

    Args:
        seed: Random seed value (default: 12)
    """
    random.seed(seed)
    np.random.seed(seed)

    
def create_folder_if_necessary(folder_path: str) -> None:
    """
    Create a folder if it doesn't exist.

    Args:
        folder_path: Path to the folder to create
    """
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def select_and_concat_all_genotypes(
    selection: pd.DataFrame, genotypes: pd.DataFrame
) -> pd.DataFrame:
    """
    Select and concatenate genotypes for all population-year combinations in selection.

    Args:
        selection: DataFrame containing population and year selections
        genotypes: DataFrame containing genotype data

    Returns:
        DataFrame with concatenated genotypes for selected populations
    """
    df_selec = pd.DataFrame()
    selection = selection[["Population", "Year"]].drop_duplicates()

    for pop, year in zip(selection.Population, selection.Year):
        df = genotypes[(genotypes.Population == pop) & (genotypes.Year == year)].copy()
        df_selec = pd.concat([df_selec, df], ignore_index=True)

    return df_selec


def select_a_given_pop(genotypes: pd.DataFrame, pop_name: str) -> pd.DataFrame:
    """
    Select genotypes for a specific population.

    Args:
        genotypes: DataFrame containing genotype data
        pop_name: Name of the population to select

    Returns:
        DataFrame with genotypes for the specified population
    """
    return genotypes[genotypes.Population == pop_name].copy()


def select_a_given_pop_year(
    genotypes: pd.DataFrame, pop_name: str, year: str
) -> pd.DataFrame:
    """
    Select genotypes for a specific population and year.

    Args:
        genotypes: DataFrame containing genotype data
        pop_name: Name of the population to select
        year: Year to select

    Returns:
        DataFrame with genotypes for the specified population and year
    """
    return genotypes[
        (genotypes.Population == pop_name) & (genotypes.Year == year)
    ].copy()


def select_a_given_pop_year_subcat(
    genotypes: pd.DataFrame, pop_name: str, year: str, subcat: str
) -> pd.DataFrame:
    """
    Select genotypes for a specific population, year, and subcategory.

    Args:
        genotypes: DataFrame containing genotype data
        pop_name: Name of the population to select
        year: Year to select
        subcat: Subcategory to select

    Returns:
        DataFrame with genotypes for the specified population, year, and subcategory
    """
    selection = select_a_given_pop_year(genotypes, pop_name, year)
    if pd.notna(subcat):
        selection = selection[selection.Subcategory == subcat]
    return selection


def conversion_two_lines_to_one_lines_genotypes(
    input_data: pd.DataFrame, meta_cols: List[str], loci_cols: List[str]
) -> pd.DataFrame:
    """
    Convert two-line genotype format to one-line format for pairwise analysis.

    This function transforms genotypes from a two-line format (where each individual
    has two rows, one for each allele) to a one-line format (where each individual
    has one row with two columns per locus: locus_1 and locus_2).

    Args:
        input_data: DataFrame with two-line genotype format
        meta_cols: List of metadata column names
        loci_cols: List of locus column names

    Returns:
        DataFrame with one-line genotype format
    """
    dfs = []
    for col in loci_cols:
        df = input_data[meta_cols + [col]].copy()
        df = (
            df.groupby(meta_cols[0], as_index=False)
            .apply(lambda x: x.sort_values(col))
            .reset_index(drop=True)
        )
        dfs.append(df)

    data = dfs[0].copy()
    for df, col in zip(dfs[1:], loci_cols[1:]):
        data[col] = df[col]

    # Split into even and odd rows (representing the two alleles)
    data_even = data.iloc[::2, :].copy()
    data_even = data_even.rename(
        columns={col: col + "_1" for col in loci_cols}
    ).reset_index(drop=True)

    data_odd = data.iloc[1::2, :].copy()
    data_odd = data_odd.rename(
        columns={col: col + "_2" for col in loci_cols}
    ).reset_index(drop=True)

    # Combine metadata with both alleles
    data = data_even[meta_cols].copy()
    for col in loci_cols:
        data[col + "_1"] = data_even[col + "_1"].copy()
        data[col + "_2"] = data_odd[col + "_2"].copy()

    return data


def get_extension_if_subcat(subcategory: str) -> str:
    """
    Get file extension string for subcategory if it exists.

    Args:
        subcategory: Subcategory value (may be NaN)

    Returns:
        Extension string (empty if subcategory is NaN, otherwise "_subcategory")
    """
    if pd.notna(subcategory):
        return f"_{subcategory}"
    return ""
