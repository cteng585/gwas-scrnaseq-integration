import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

import sh

import utils
from classes import BFileSet


def prune_problem_variants(bfiles: List, work_dir: str) -> None:
    """
    extract error-causing variants from the .missnp file and create pruned versions
    of the bfile sets that do not contain those error-causing variants

    :param bfiles: the list of bfile sets
    :param work_dir: the working directory for the reference-making process
    :return:
    """
    with open(f"{work_dir}/merge-merge.missnp", "r") as infile:
        merge_problem_variants = list()
        for line in infile.readlines():
            variant_search = re.search(r"rs\d+", line)
            if variant_search:
                merge_problem_variants.append(variant_search.group(0))

    os.rename(f"{work_dir}/merge-merge.missnp", f"{work_dir}/failed_merge.missnp")
    os.rename(f"{work_dir}/merge.log", f"{work_dir}/failed_merge.log")

    with open(f"{work_dir}/exclude_merge_variants.txt", "w") as outfile:
        outfile.write("\n".join(merge_problem_variants))

    flags = dict()
    for bfile in bfiles:
        print(f"Pruning problem variants from {bfile}")
        flags["--bfile"] = f"{bfile}"
        flags["--exclude"] = f"{work_dir}/exclude_merge_variants.txt"
        flags["--out"] = f"{work_dir}/{os.path.basename(bfile)}.pruned"

        utils.run_plink(no_arg_flags=["--make-bed"], yes_arg_flags=flags)


def keep_bfiles(bfile_sets: Dict[str, BFileSet], keep: List[str]) -> None:
    """
    only keep bfile sets that will merged

    :params keep: expect string patterns corresponding to bfile sets that should be kept
    :return: None
    """
    to_keep = re.compile("|".join([f"({pattern})" for pattern in keep]))

    to_remove = list()
    for bfile_set in bfile_sets:
        if not re.search(to_keep, bfile_set):
            to_remove.append(bfile_set)

    for bfile_set in to_remove:
        del bfile_sets[bfile_set]


def merge_bfiles(
    bfile_sets: Dict[str, BFileSet],
    work_dir: str,
    output_prefix: str = "merge",
    keep: List[str] = None,
) -> Tuple[Path, Path, Path]:
    """
    merge a list of bfile sets

    :param bfile_sets: a Dict of available bfile sets to merge
    :param work_dir: the working directory for the reference-making process
    :param output_prefix: what the path stem of the merged bfile set should be
    :param keep: expect string patterns corresponding to bfile sets that should be kept
    :return:
    """
    keep_bfiles(bfile_sets, keep)
    to_merge = list(bfile_sets.keys())

    # make a temp file containing bfiles to merge into the main bfile
    # per PLINK's workflow for merging bfile sets
    with open(f"{work_dir}/merge_list.txt", "w") as outfile:
        outfile.write("\n".join(to_merge[1:]))

    try:
        print(f"Attempting to merge bfile sets {to_merge}")
        flags = dict()
        flags["--bfile"] = to_merge[0]
        flags["--merge-list"] = f"{work_dir}/merge_list.txt"
        flags["--out"] = f"{work_dir}/{output_prefix}"
        utils.run_plink(no_arg_flags=["--make-bed"], yes_arg_flags=flags)

    # handle the case where the merge doesn't work due to a small-subset of variants
    except sh.ErrorReturnCode:
        print(
            f"Failed to merge bfile sets {to_merge}. Attempting to prune problem variants before "
            f"trying again"
        )
        prune_problem_variants(to_merge, work_dir)
        pruned_bfiles = [
            f"{work_dir}/{os.path.basename(bfile)}.pruned" for bfile in bfile_sets
        ]
        flags["--bfile"] = pruned_bfiles[0]
        print(f"Retrying bfile set merge...")
        utils.run_plink(no_arg_flags=["--make-bed"], yes_arg_flags=flags)
        print(f"...success")

    # check that the merged files exist
    merged_bed = Path(f"{work_dir}/{output_prefix}.bed")
    merged_bim = Path(f"{work_dir}/{output_prefix}.bim")
    merged_fam = Path(f"{work_dir}/{output_prefix}.fam")

    assert merged_bed.exists(), f"Something went wrong; {merged_bed} doesn't exist"
    assert merged_bim.exists(), f"Something went wrong; {merged_bim} doesn't exist"
    assert merged_fam.exists(), f"Something went wrong; {merged_fam} doesn't exist"

    return merged_bed, merged_bim, merged_fam
