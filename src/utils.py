import csv
import os
from pathlib import Path
from typing import Tuple, Union

import secrets


def detect_delimiter(infile: Path) -> str:
    """
    find the column/field delimiter in a given file

    :param infile: file to detect the column/field delimiter for
    :return: the column/field delimiter as a string
    """
    with open(infile, "r") as infile:
        delimiter = str(csv.Sniffer().sniff(infile.readline()).delimiter)

    return delimiter


def make_dir(parent_dir: str, dir_name: str = "") -> str:
    """
    make a directory under a given parent directory if it doesn't exist. makes a
    child directory with a random hash name if a directory name isn't provided

    :param parent_dir: the parent directory to create a child directory under
    :param dir_name: the name of the child directory
    :return: name of the directory that was (attempted to be) created
    """
    if not os.path.exists(f"{parent_dir}/{dir_name}"):
        print(f"making directory {parent_dir}/{dir_name}")
        os.makedirs(f"{parent_dir}/{dir_name}")
    else:
        print(f"directory {parent_dir}/{dir_name} already exists")

    return f"{parent_dir}/{dir_name}"


def make_work_dir(parent_dir: str) -> str:
    """
    make a temporary working directory with a randomly generated hex name

    :param parent_dir: the parent directory to make the working directory under
    :return: name of the directory that was (attempted to be) created
    """
    work_dir = secrets.token_hex(nbytes=16)
    return make_dir(f"{parent_dir}", work_dir)


def move_output(output_dir: Union[str, Path], *args) -> None:
    """
    move output of a command to the correct output directory

    :param output_dir: location to move final output files to
    :return:
    """
    for file_path in args:
        if isinstance(file_path, str):
            file_path = Path(file_path)

        assert file_path.exists(), f"Expected output file {file_path} to exist but it can't be found"

        os.rename(file_path, f"{output_dir}/{file_path.name}")


def setup(output_dir: str) -> Tuple[str, str]:
    """
    setup directories that will be used for workflow commands

    :param output_dir: where the final output should be written to
    :return: the final output directory and the working directory for this module
    """
    if "/" not in output_dir:
        parent_dir = os.getcwd()
        output_dir = make_dir(parent_dir, output_dir)
    else:
        parent_dir = os.path.dirname(output_dir)
        output_dir = make_dir(parent_dir, os.path.basename(output_dir))

    tmp_dir = make_dir(parent_dir, "tmp")

    return output_dir, tmp_dir
