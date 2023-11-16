from enum import Enum
from pathlib import Path
from pydantic import BaseModel


class Ancestry(Enum):
    """
    define the 1k genome ancestries that are supported for merging into a custom
    reference panel
    """

    AFRICAN = "afr"
    AMERICAN = "amr"
    EAST_ASIAN = "eas"
    EUROPEAN = "eur"
    SOUTH_ASIAN = "sas"


class BFileType(Enum):
    """
    defines the useable bfile types in a 1k genome data directory
    """

    BED = ".bed"
    BIM = ".bim"
    FAM = ".fam"


class BFileSet(BaseModel):
    """
    defines the expected structure of a PLINK bfile set
    """

    BED: Path
    BIM: Path
    FAM: Path
    ANCESTRY: Ancestry
