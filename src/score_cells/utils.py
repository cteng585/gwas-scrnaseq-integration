from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd

from src.utils import detect_delimiter


def load_gene_name_map(
        gene_name_map: Union[str, Path],
        header: bool,
        delimiter: Optional[str] = None,
        id_column: Union[str, int] = 0,
        name_column: Union[str, int] = 1,
) -> Dict:
    """
    read in a CSV or TSV to map the gene IDs used during the MAGMA annotation step to
    HGNC gene names

    if an ID column index and a name column index aren't explicitly passed, expect the first
    column to be the IDs to convert from and expect the second column to be the names to
    convert to

    :param gene_name_map: the CSV to TSV that maps gene IDs to gene names
    :param header: whether the gene_name_map file has a header
    :param delimiter: the column delimiter of the gene_name_map
    :param id_column: which column in the gene_name_map file contains the IDs to convert from
    :param name_column: which column in the gene_name_map file contains the names to convert to
    :return: a dict mapping gene IDs to gene names
    """
    if not isinstance(gene_name_map, Path):
        gene_name_map = Path(gene_name_map)

    if not delimiter:
        delimiter = detect_delimiter(gene_name_map)

    if header:
        gene_name_map = pd.read_csv(gene_name_map, sep=delimiter)
    else:
        gene_name_map = pd.read_csv(gene_name_map, sep=delimiter, header=None)

    try:
        if isinstance(id_column, str):
            gene_ids = gene_name_map.loc[:, id_column]
        elif isinstance(id_column, int):
            gene_ids = gene_name_map.iloc[:, id_column]
    except Exception as e:
        raise e

    try:
        if isinstance(name_column, str):
            gene_names = gene_name_map.loc[:, name_column]
        elif isinstance(name_column, int):
            gene_names = gene_name_map.iloc[:, name_column]
    except Exception as e:
        raise e


    gene_name_map = (
        pd.concat([gene_ids, gene_names], axis=1)
        .set_index(gene_name_map.columns[id_column])
        .to_dict(orient="dict")
        .get(gene_name_map.columns[name_column])
    )

    return gene_name_map


def munge_magma(
        associated_genes: Union[str, Path],
        trait: str,
        num_genes: int,
        work_dir: Path,
        delimiter: Optional[str] = r"\s+",
        header: bool = True,
        gene_column: Union[str, int] = "GENE",
        pvalue_column: Union[str, int] = "P",
        gene_name_map: Optional[Dict] = None
) -> Path:
    """
    munge a MAGMA .genes.out file from SNP-wise gene analysis to be useable with scDRS

    :param associated_genes: a MAGMA .genes.out file containing at least gene names/IDs and their p-values for the significance of association of the gene with a trait of interest
    :param trait: the trait that the genes are associated with
    :param num_genes: the number of top genes to consider when constructing a gene set for the trait
        of interest
    :param work_dir: the working directory to output the munged gene set to
    :param delimiter: the column/field delimiter of associated_genes
    :param gene_column: which column of the MAGMA .genes.out file contains the gene IDs/names
    :param pvalue_column: which column of the MAGMA .genes.out file contains the gene significance values
    :param header: whether the MAGMA .genes.out file has a header
    :param gene_name_map: optional dict to convert gene IDs to gene names
    :return: the path to the output file
    """
    if not isinstance(associated_genes, Path):
        associated_genes = Path(associated_genes)

    if not delimiter:
        delimiter = detect_delimiter(associated_genes)

    if header:
        trait_gene_set = pd.read_csv(associated_genes, sep=delimiter)
        trait_gene_set.rename(
            columns={
                gene_column: "GENE",
                pvalue_column: "P",
            }, inplace=True
        )
    else:
        trait_gene_set = pd.read_csv(associated_genes, header=None, sep=delimiter)
        assert isinstance(gene_column, int), "If the trait-associated gene file has no header, " \
                                             "requires an int be passed as the gene column"
        assert isinstance(pvalue_column, int), "If the trait-associated gene file has no header, " \
                                               "requires an int be passed as the p-value column"
        trait_gene_set.rename(
            columns={
                trait_gene_set.columns[gene_column]: "GENE",
                trait_gene_set.columns[pvalue_column]: "P",
            }
        )

    if gene_name_map:
        trait_gene_set[gene_column].replace(gene_name_map, inplace=True)

    trait_gene_set.sort_values(by=["P"], inplace=True)
    trait_gene_set = ",".join(trait_gene_set.loc[:num_genes, "GENE"])

    gene_set_out = Path(f"{work_dir}/{associated_genes.stem}.gs")
    header = ["TRAIT", "GENESET"]
    with open(gene_set_out, "w") as outfile:
        outfile.write("\t".join(header) + "\n")
        outfile.write("\t".join([trait, trait_gene_set]))

    return gene_set_out
