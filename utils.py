import os

import pandas as pd


def create_folder_if_necessary(folder_path: str) -> None:
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def select_and_concat_all_genotypes(
    selection: pd.DataFrame, genotypes: pd.DataFrame
) -> pd.DataFrame:
    df_selec = pd.DataFrame()
    selection = selection[["Population", "Year"]].drop_duplicates()
    for _pop, _year in zip(selection.Population, selection.Year):
        df = genotypes[
            (genotypes.Population == _pop) & (genotypes.Year == _year)
        ].copy()
        df_selec = pd.concat([df_selec, df], ignore_index=True)
    return df_selec


def select_a_given_pop(genotypes: pd.DataFrame, pop_name: str) -> pd.DataFrame:
    return genotypes[genotypes.Population == pop_name].copy()


def select_a_given_pop_year(
    genotypes: pd.DataFrame, pop_name: str, year: str
) -> pd.DataFrame:
    return genotypes[
        (genotypes.Population == pop_name) & (genotypes.Year == year)
    ].copy()


def select_a_given_pop_year_subcat(
    genotypes: pd.DataFrame, pop_name: str, year: str, subcat: str
) -> pd.DataFrame:
    selection = select_a_given_pop_year(genotypes, pop_name, year)
    if pd.notna(subcat):
        selection = selection[selection.Subcategory == subcat]
    return selection


def conversion_two_lines_to_one_lines_genotypes(
    input_data: pd.DataFrame, meta_cols: list[str], loci_cols: list[str]
) -> pd.DataFrame:
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

    data_even = data.iloc[::2, :].copy()
    data_even = data_even.rename(
        columns={col: col + "_1" for col in loci_cols}
    ).reset_index(drop=True)

    data_odd = data.iloc[1::2, :].copy()
    data_odd = data_odd.rename(
        columns={col: col + "_2" for col in loci_cols}
    ).reset_index(drop=True)
    data = data_even[meta_cols].copy()

    for col in loci_cols:
        data[col + "_1"] = data_even[col + "_1"].copy()
        data[col + "_2"] = data_odd[col + "_2"].copy()
    return data


def get_extension_if_subcat(_sub: str) -> str:
    _ext = ""
    if pd.notna(_sub):
        _ext = f"_{_sub}"
    return _ext