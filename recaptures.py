import pandas as pd
from config import FileConfiguration
from loguru import logger


def get_recaptures(config: FileConfiguration) -> None:
    selection = config.pops_to_select.copy()

    if config.agg_type == "all":
        df_pairwise = pd.read_csv(
            f"{config.output_path}/pairwise_differences/pairwise_differences_{config.selection_name}_all.csv",
            sep=";",
        )
        df_pairwise = df_pairwise.set_index("Sample")
        get_recaptures_sample(config, df_pairwise, config.agg_type, True)

    elif config.agg_type == "pops":
        for _pop in selection.Population.unique():
            df_pairwise = pd.read_csv(
                f"{config.output_path}/pairwise_differences/pairwise_differences_{config.selection_name}_{_pop}.csv",
                sep=";",
            )
            df_pairwise = df_pairwise.set_index("Sample")
            get_recaptures_sample(config, df_pairwise, _pop, True)

    elif config.agg_type == "pop_years":
        selection = selection[["Population", "Year"]].drop_duplicates()
        for _pop, _year in zip(selection.Population, selection.Year):
            df_pairwise = pd.read_csv(
                f"{config.output_path}/pairwise_differences/pairwise_differences_{config.selection_name}_{_pop}_{_year}.csv",
                sep=";",
            )
            df_pairwise = df_pairwise.set_index("Sample")
            get_recaptures_sample(config, df_pairwise, f"{_pop}_{_year}", True)

    else:  # subcat
        for _pop, _year, _sub in zip(
            selection.Population, selection.Year, selection.Subcategory
        ):
            df_pairwise = pd.read_csv(
                f"{config.output_path}/pairwise_differences/pairwise_differences_{config.selection_name}_{_pop}_{_year}_{_sub}.csv",
                sep=";",
            )
            df_pairwise = df_pairwise.set_index("Sample")
            get_recaptures_sample(config, df_pairwise, f"{_pop}_{_year}_{_sub}", True)


def get_recaptures_sample(
    config: FileConfiguration, df_pairwise: pd.DataFrame, name: str, save_file: False
) -> pd.DataFrame:

    cols = list(df_pairwise.columns)
    res = pd.DataFrame()
    for _col in cols:
        recap = df_pairwise.index[df_pairwise[_col] == 0].tolist()
        if len(recap) > 0:
            for _re in recap:
                res = pd.concat(
                    [res, pd.DataFrame({"ind1": _col, "ind2": _re}, index=[0])],
                    ignore_index=True,
                )

    if len(res) > 0:
        res["combination"] = res.apply(lambda x: set([x.ind1, x.ind2]), axis=1)
        df = pd.DataFrame(
            res["combination"].drop_duplicates().reset_index(drop=True).to_list(),
            columns=["ind1", "ind2"],
        )

        if save_file:
            df["ind1_pop"] = df.ind1.apply(lambda x: config.dict_pop_samples[x])
            df["ind2_pop"] = df.ind2.apply(lambda x: config.dict_pop_samples[x])

            df.to_csv(
                f"{config.output_path}/recaptures/recaptures_{name}.csv",
                sep=";",
                index=False,
            )
            logger.info("Recaptures. File saved.")
        return df
    else:
        if save_file:
            logger.info("No recaptures.")
        return pd.DataFrame()
