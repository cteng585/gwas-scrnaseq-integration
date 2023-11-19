from pathlib import Path
from typing import Optional, Tuple, Union

from utils import run_magma


def annotate_variants(
    gene_loc: Path,
    snp_loc: Path,
    output_prefix: Path,
    annotation_window: Optional[Tuple[int, int]],
) -> Tuple[Path, Path]:
    """
    from the MAGMA documentation - "the annotation step is a pre-processing step prior to the
    actual analysis, in which SNPs are mapped to genes. The mapping is based on genomic location,
    assigning an SNP to a gene if the SNP’s location falls inside the region provided for each gene"

    :param gene_loc: path to the gene location file
    :param snp_loc: path to the SNP location file
    :param output_prefix: the prefix for each output of the MAGMA annotation step
    :param annotation_window: window to use around each gene to associate SNPs outside the main gene body
        with the gene
    :return: the output annotation files
    """

    no_arg_flags = ["--annotate"]
    yes_arg_flags = {
        "window": f"{annotation_window[0]},{annotation_window[1]}",
        "--gene-loc": gene_loc,
        "--snp-loc": snp_loc,
        "--out": output_prefix,
    }

    # this comes after since the order of arguments in `yes_arg_flags` is respected when passed to sh.Command
    if not annotation_window:
        del yes_arg_flags["window"]

    run_magma(no_arg_flags, yes_arg_flags)
    annotation_output = Path(f"{str(output_prefix)}.genes.annot")
    annotation_log = Path(f"{str(output_prefix)}.log")

    assert (
        annotation_output.exists()
    ), f"Something went wrong in the annotation step; no annotation output was produced"
    assert (
        annotation_log.exists()
    ), f"Something went wrong in the annotation step; no annotation log was produced"

    return annotation_output, annotation_log


def gene_analysis(
    bfile: Path,
    gene_annot: Path,
    gwas: Path,
    variant_id: str,
    pval: str,
    n: Union[int, str],
    output_prefix: Path,
) -> Tuple[Path, Path, Path, Path]:
    """
    wrapper for conducting gene analysis on GWAS results in the form of an SNP-wise p-value file

    :param bfile: the reference data used to estimate LD between SNPs. passed as the prefix of a PLINK
        bfile set
    :param gene_annot: the gene annotation file linking variant IDs to genes
    :param gwas: the GWAS p-value file. should have at least a column with variant IDs and a column with
        p-values per variant ID
    :param variant_id: the column in the GWAS p-value file containing variant IDs
    :param pval: the column in the GWAS p-value file containin variant p-values
    :param n: either an int or a str. if an int, n is the total sample size of the GWAS study for
        case-control studies. if a str, specifies a column in the p-value file that contains the sample
        size used per SNP
    :param output_prefix: where output from this command should be output
    :return: the output of the MAGMA gene-analysis step. should be
    """
    no_arg_flags = []
    yes_arg_flags = {
        "--bfile": bfile,
        "--gene-annot": gene_annot,
        "--pval": gwas,
        "snp-id": variant_id,
        "pval": pval,
    }

    if isinstance(n, int):
        yes_arg_flags["N"] = n
    elif isinstance(n, str):
        yes_arg_flags["ncol"] = n
    else:
        raise NotImplementedError

    yes_arg_flags["--out"] = output_prefix

    run_magma(no_arg_flags, yes_arg_flags)

    gene_analysis_log = Path(f"{str(output_prefix)}.log")
    gene_analysis_supplemental_log = Path(f"{str(output_prefix)}.log.suppl")
    gene_analysis_raw = Path(f"{str(output_prefix)}.genes.raw")
    gene_analysis_out = Path(f"{str(output_prefix)}.genes.out")

    assert (
        gene_analysis_log.exists()
    ), f"Something went wrong in the gene analysis step; no analysis log was produced"
    assert (
        gene_analysis_supplemental_log.exists()
    ), f"Something went wrong in the gene analysis step; no analysis supplemental log was produced"
    assert (
        gene_analysis_raw.exists()
    ), f"Something went wrong in the gene analysis step; no analysis raw output was produced"
    assert (
        gene_analysis_out.exists()
    ), f"Something went wrong in the gene analysis step; no analysis final output was produced"

    return (
        gene_analysis_out,
        gene_analysis_raw,
        gene_analysis_log,
        gene_analysis_supplemental_log,
    )
