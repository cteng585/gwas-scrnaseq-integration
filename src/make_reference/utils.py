import os
import re
from pathlib import Path
from typing import Dict, List, Union

import sh
from src.make_reference.classes import Ancestry, BFileType, BFileSet

# TODO: change strings to pathlib objects
# TODO: logging


def use_bfiles(bfile_dir: Union[Path, str], *args: Ancestry) -> Dict[str, BFileSet]:
    """
    make a list of bfiles to attempt to merge

    :param bfile_dir: path to where the reference bfiles are stored
    :param args: expects a list of Ancestry enums
    :return: list of bfiles to attempt to merge
    """
    bfile_sets = dict()
    for root, _, files in os.walk(bfile_dir):
        for file in files:
            bfile_path = Path(f"{root}/{file}")
            try:
                # validate that the bfile type is a standard PLINK bfile type by trying to create
                # a BFileType enum based on the file extension
                bfile_type = BFileType(bfile_path.suffix)
                bfile_prefix = f"{bfile_path.parent}/{bfile_path.stem}"

                # validate that the bfile ancestry is supported by trying to create
                # an Ancestry enum based on a regex search for standard 1k genome ancestry strings
                ancestry = re.search(
                    r"(?<=[^A-Za-z])((amr)|(eas)|(eur)|(sas)|(afr))(?=[^A-Za-z])",
                    str(bfile_path),
                    flags=re.IGNORECASE,
                )
                ancestry = Ancestry(ancestry.group().lower())

                # only add bfile paths for the ancestries that need to be merged
                if ancestry not in args:
                    continue

            # only handle recognized PLINK bfile types (bed, bim fam)
            except ValueError:
                continue

            if bfile_prefix not in bfile_sets:
                bfile_sets[bfile_prefix] = {
                    "BED": None,
                    "BIM": None,
                    "FAM": None,
                    "ANCESTRY": None,
                }
            elif bfile_prefix in bfile_sets and all(bfile_sets[bfile_prefix].values()):
                bfile_files_w_prefix = "\n".join(bfile_sets[bfile_prefix].values())
                raise KeyError(
                    f"All bfiles with stem {bfile_prefix} have already been found:\n"
                    f"---\n{bfile_files_w_prefix}\n--- "
                    f"but an additional bfile: \n{bfile_path}\n was found that seems to be for "
                    f"the same ancestry. Check that only one bfile set for {bfile_prefix} "
                    f"exists in the passed reference bfile directory"
                )

            match bfile_type:
                case BFileType.BED:
                    bfile_sets[bfile_prefix]["BED"] = bfile_path

                case BFileType.BIM:
                    bfile_sets[bfile_prefix]["BIM"] = bfile_path

                case BFileType.FAM:
                    bfile_sets[bfile_prefix]["FAM"] = bfile_path

                case _:
                    raise ValueError(f"Unexpected bfile type for bfile {bfile_path}")

            bfile_sets[bfile_prefix]["ANCESTRY"] = ancestry

    # convert each sub-dictionary to a BFileSet for more comprehensible access to bfile paths
    # also validates that each PLINK bfile sets is a complete set of 3 bfiles + an ancestry
    for bfile_stem in bfile_sets:
        bfile_sets[bfile_stem] = BFileSet.model_validate(bfile_sets[bfile_stem])

    print("Available bfile sets that match ancestry: ")
    print("\n".join(bfile_sets))

    return bfile_sets


def run_plink(no_arg_flags: List, yes_arg_flags: Dict):
    """
    run PLINK with command-line arguments

    :param no_arg_flags: PLINK flags that do not take an argument (e.g. --make-bed)
    :param yes_arg_flags: PLINK flags that do take an argument (e.g. --out)
    :return:
    """
    assert "/plink" in os.environ.get(
        "PATH"
    ), f"Couldn't find a plink executable in the $PATH variable"

    plink = sh.Command("plink")
    command_list = [*no_arg_flags]
    for key, value in yes_arg_flags.items():
        if re.search(r"^--[a-z]", key):
            command_list += [f"{str(key)}", value]
        else:
            command_list += [f"{str(key)}={value}"]

    plink(command_list)
