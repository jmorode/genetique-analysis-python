import json
import os
import sys

import pandas as pd
from loguru import logger
from utils import (
    conversion_two_lines_to_one_lines_genotypes,
    create_folder_if_necessary,
)


class FileConfiguration:
    def __init__(self, project_name: str, agg_type: str, selection_name: str):
        assert agg_type in ["all", "pops", "pop_years", "subcat"]

        self.project_name = project_name
        self.agg_type = agg_type
        self.selection_name = selection_name

        # INIT
        self.config = {}

        self.input_path = f"./projects/{self.project_name}/inputs"
        self.path_to_config_file = f"{self.input_path}/config.json"
        self.output_path = f"./projects/{self.project_name}/outputs"
        self.genotypes_two_lines = None
        self.genotypes_data_pairwise = None
        self.bool_config = True
        self.pops_to_select = None
        self.dict_pop_samples = None
        self.loci_list = None

    def _no_config(self) -> None:
        self.bool_config = False
        logger.info("Config file does not exist!")

    def _check_config_file(self) -> None:
        if not os.path.exists(self.path_to_config_file):
            self._no_config()

    def _retrieve_config(self) -> None:
        with open(self.path_to_config_file, encoding="utf-8") as f:
            data = json.load(f)

        for key in data.keys():
            self.config[key] = data[key]

    def initialize_folders(self) -> None:
        for _sub in [
            "raw_data",
            "heterozygosity",
            "pairwise_differences",
            "pairwise_differences/plots/",
            "recaptures",
            "ml_relate",
            "test_stats",
            "error_recapture_tests"
        ]:
            create_folder_if_necessary(f"{self.output_path}/{_sub}/")

    def check_population_year(self, selection: pd.DataFrame) -> None:
        pop_year_geno = self.genotypes_two_lines[
            ["Population", "Year"]
        ].drop_duplicates()
        for _pop, _year in zip(selection.Population, selection.Year):
            assert (
                len(
                    pop_year_geno[
                        (pop_year_geno.Population == _pop)
                        & (pop_year_geno.Year == _year)
                    ]
                )
                != 0
            ), (
                f"{_pop} / {_year} in selection.csv doesn't exist in" + " genotypes.csv"
            )

    def check_population_year_subcat(self, selection: pd.DataFrame) -> None:
        pop_year_geno = self.genotypes_two_lines[
            ["Population", "Year", "Subcategory"]
        ].drop_duplicates()
        for _pop, _year, _sub in zip(
            selection.Population, selection.Year, selection.Subcategory
        ):
            if pd.notna(_sub):
                assert (
                    len(
                        pop_year_geno[
                            (pop_year_geno.Population == _pop)
                            & (pop_year_geno.Year == _year)
                            & (pop_year_geno.Subcategory == _sub)
                        ]
                    )
                    != 0
                ), (
                    f"{_pop} / {_year} / {_sub} in selection.csv doesn't exist"
                    + " in genotypes.csv"
                )

    def check_loci_names(self, conversion: pd.DataFrame) -> None:
        list_loci_genotypes = self.genotypes_two_lines.columns
        list_loci_genotypes = list_loci_genotypes.drop(
            ["Population", "Sample", "Year", "Subcategory"]
        )
        list_loci_conversion = conversion["Loci"].unique()
        assert set(list_loci_conversion) == set(list_loci_genotypes), (
            f"Loci in conversion.csv {set(list_loci_conversion)}"
            + " and in genotypes.csv {set(list_loci_genotypes)} are different"
        )

    def convert_raw_genotypes(self, conversion_table):
        for _loci in conversion_table.Loci.unique():
            subconv_table = conversion_table[conversion_table.Loci == _loci]
            subconv_table = subconv_table.drop("Loci", axis=1)
            subconv_dict = subconv_table.set_index("allele_raw").to_dict("dict")

            self.genotypes_two_lines[_loci] = self.genotypes_two_lines[_loci].apply(
                lambda x: (
                    int(subconv_dict["allele_corrected"][x]) if pd.notnull(x) else None
                )
            )
        self.genotypes_two_lines.to_csv(
            f"{self.output_path}/raw_data/genotypes.csv", sep=";", index=False
        )

    def convert_genotypes_two_lines_to_one_line(self) -> None:
        input_data = self.genotypes_two_lines.copy()
        loci_cols = input_data.columns
        meta_cols = ["Sample", "Population", "Year", "Subcategory"]
        loci_cols = loci_cols.drop(meta_cols)

        data = conversion_two_lines_to_one_lines_genotypes(
            input_data, meta_cols, loci_cols
        )
        data.to_csv(
            f"{self.output_path}/raw_data/genotypes_pairwise_diff.csv",
            sep=";",
            index=False,
            header=False,
        )
        self.genotypes_data_pairwise = data.copy()

    def retrieve_and_prepare_inputs(self, selection_name: str):
        # get genotypes full file
        print(self.input_path)
        print(os.getcwd())
        if os.path.exists(f"{self.input_path}/genotypes.csv"):
            self.genotypes_two_lines = pd.read_csv(
                f"{self.input_path}/genotypes.csv", sep=";"
            )
        else:
            logger.error("genotypes.csv file does not exist ! Exiting...")
            sys.exit()

        # convert genotypes
        if os.path.exists(f"{self.input_path}/conversion.csv"):
            conversion = pd.read_csv(f"{self.input_path}/conversion.csv", sep=";")
            self.check_loci_names(conversion)
            self.convert_raw_genotypes(conversion)
        else:
            logger.info(
                "No conversion.csv file was loaded, no allele conversion will be applied."
            )

        # prepare one line genotypes
        self.convert_genotypes_two_lines_to_one_line()

        # save all pop/year/subcat/combinaisons
        self.save_all_pop_year_subcategory_combinations()

        # get pops to select
        self.pops_to_select = pd.read_csv(
            f"{self.input_path}/selection_{selection_name}.csv", sep=";"
        )
        # check selection variables
        self.check_population_year(self.pops_to_select)
        if self.agg_type == "subcat":
            self.check_population_year_subcat(self.pops_to_select)

    def initialize_global(self):
        self._check_config_file()
        if self.bool_config:
            self._retrieve_config()
        self.initialize_folders()
        self.retrieve_and_prepare_inputs(self.selection_name)
        self.init_dict_pop_samples()
        self.initialize_loci_list()

    def save_all_pop_year_subcategory_combinations(self):
        combinations = self.genotypes_two_lines[
            ["Population", "Year", "Subcategory"]
        ].drop_duplicates()
        combinations.to_csv(
            f"{self.input_path}/select_everything.csv", sep=";", index=False
        )

    def init_dict_pop_samples(self):
        pop_samples = self.genotypes_two_lines.copy()

        if self.agg_type == "pops":
            pop_samples["legend"] = pop_samples["Population"]
        elif self.agg_type in ["all", "pop_years"]:
            pop_samples["legend"] = (
                pop_samples["Population"] + " - " + pop_samples["Year"].astype(str)
            )
        else:
            pop_samples["legend"] = (
                pop_samples["Population"]
                + " - "
                + pop_samples["Year"].astype(str)
                + " - "
                + pop_samples["Subcategory"].astype(str)
            )

        pop_samples = pop_samples[["Sample", "legend"]].drop_duplicates()
        self.dict_pop_samples = dict(zip(pop_samples["Sample"], pop_samples["legend"]))

    def initialize_loci_list(self):
        loci_list = self.genotypes_two_lines.columns
        self.loci_list = loci_list.drop(["Population", "Sample", "Year", "Subcategory"])
