import os
import re
from typing import Dict, List, Optional
from pathlib import Path

import pandas as pd
import sh

from src.utils import detect_delimiter


def run_magma(no_arg_flags: List, yes_arg_flags: Dict):
    """
    run MAGMA with command-line arguments

    :param no_arg_flags: MAGMA flags that do not take an argument (e.g. --make-bed)
    :param yes_arg_flags: MAGMA flags that do take an argument (e.g. --out)
    :return:
    """
    assert "/magma" in os.environ.get(
        "PATH"
    ), f"Couldn't find a magma executable in the $PATH variable"

    magma = sh.Command("magma")
    command_list = [*no_arg_flags]
    for key, value in yes_arg_flags.items():
        if re.search(r"^--[a-z]", key):
            command_list += [f"{str(key)}", value]
        else:
            command_list += [f"{str(key)}={value}"]

    magma(command_list)


def munge_gwas(
    gwas_pval: Path,
    variant_id: Optional[str] = None,
    pval: Optional[str] = None,
    n: Optional[str] = None,
    chromosome: Optional[str] = None,
    bp: Optional[str] = None,
    ref_allele: Optional[str] = None,
    effect_allele: Optional[str] = None,
):
    """
    format a GWAS p-value file such that it's useable by MAGMA

    :param gwas_pval: a GWAS p-value file
    :param variant_id: the variant ID column in the GWAS file. if no value is provided, attempt to infer
    :param pval: the p-value column in the GWAS file. if no value is provided, attempt to infer
    :param n: the sample size column in the GWAS file. if no value is provided, attempt to infer
    :param chromosome: the chromosome column. used to generate a filler ID if a variant is missing its ID.
        if no value is provided, attempt to infer
    :param bp: the base pair column. used to generate a filler ID if a variant is missing its ID.
        if no value is provided, attempt to infer
    :param ref_allele: the reference allele column. used to generate a filler ID if a variant is
        missing its ID. if no value is provided, attempt to infer
    :param effect_allele: the effect allele column. used to generate a filler ID if a variant is
        missing its ID. if no value is provided, attempt to infer
    :return:
    """
    # find the columns in the p-value file that MAGMA actually uses
    with open(gwas_pval, "r") as infile:
        header = infile.readline()

    delimiter = detect_delimiter(gwas_pval)
    header = header.strip().split(delimiter)

    column_patterns = {
        "pval": pval if pval else re.compile(r"^p.?(val)?", flags=re.IGNORECASE),
        "variant_id": variant_id
        if variant_id
        else re.compile(r"^(snp)|(variant)", flags=re.IGNORECASE),
        "n": n if n else re.compile(r"^n$", flags=re.IGNORECASE),
        "chr": chromosome if chromosome else re.compile(r"^chr", flags=re.IGNORECASE),
        "bp": bp if bp else re.compile(r"^(bp)|(base.?pair)"),
        "ref_allele": ref_allele
        if ref_allele
        else re.compile(r"^(ref(erence)?)|(other)|(major)"),
        "effect_allele": effect_allele
        if effect_allele
        else re.compile(r"^(effect)|(minor)"),
    }

    for idx, column in enumerate(header):
        for key in column_patterns:
            if isinstance(column_patterns[key], re.Pattern) and column_patterns[
                key
            ].search(column):
                column_patterns[key] = column

    # assert that columns for the variant ID and p-value exist in the file
    assert isinstance(column_patterns["pval"], str)
    assert isinstance(column_patterns["variant_id"], str)
    assert isinstance(column_patterns["effect_allele"], str)
    core_cols = [
        column_patterns["pval"],
        column_patterns["variant_id"],
        column_patterns["effect_allele"],
    ]
    core_dtypes = {
        column_patterns["pval"]: float,
        column_patterns["variant_id"]: str,
        column_patterns["effect_allele"]: str,
    }

    # if a sample size column is included per SNP, include that too for parsing the file
    if isinstance(column_patterns["n"], str):
        core_cols.append(column_patterns["n"])
        core_cols[column_patterns["n"]] = "Int64"

    # try to lazily load the gwas p-value file by not loading the chromosome, base pair,
    # reference allele, or effect allele information unless variant IDs are missing
    gwas = pd.read_csv(
        gwas_pval, sep=delimiter, usecols=core_cols, dtype=core_dtypes
    )

    # if there are any missing variant IDs, fill using the hg19 format:
    # [CHR]_[BP]_[REF ALLELE]_[EFFECT ALLELE]
    if any(gwas[variant_id].isna()):
        filler_cols = [
            column_patterns["chr"],
            column_patterns["bp"],
            column_patterns["ref_allele"],
        ]
        filler_dtypes = {
            column_patterns["chr"]: str,
            column_patterns["bp"]: "Int64",
            column_patterns["ref_allele"]: str,
        }

        extension = pd.read_csv(
            gwas_pval, sep=delimiter, usecols=filler_cols, dtype=filler_dtypes
        )

        gwas = pd.concat([gwas, extension], axis=1)
        null_id_mask = gwas[column_patterns["variant_id"]].isna()
        gwas_null_id = gwas[null_id_mask]

        filler_cols.append(column_patterns["effect_allele"])
        gwas.loc[null_id_mask, column_patterns["variant_id"]] = gwas_null_id[
            filler_cols
        ].apply(lambda x: "_".join(x.values.astype(str)), axis=1)

    # need to remove variants where there isn't an effect allele
    gwas = gwas[~gwas[column_patterns["effect_allele"]].isna()]

    # MAGMA doesn't accept p-values <= 1e-308
    # convert p-values below this value to the floor value
    gwas.loc[
        gwas[column_patterns["pval"]] < 9.9999999999999999e-307, column_patterns["pval"]
    ] = 9.9999999999999999e-307

    # the annotation step requires the first three columns of the SNP location file
    # to be the SNP ID, chromosome, and base pair location. any additional columns are ignored
    gwas = gwas[[
        column_patterns["variant_id"],
        column_patterns["chr"],
        column_patterns["bp"],
        column_patterns["pval"],
    ]]

    return gwas
