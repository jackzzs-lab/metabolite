from dataclasses import dataclass
from typing import Tuple

@dataclass(frozen=True)
class Pathway:
    name: str
    smpdb_id: str = None
    kegg_id: str = None

@dataclass(frozen=True)
class Enzyme:
    name: str
    uniprot_id: str
    gene_name: str = None

@dataclass(frozen=True)
class Metabolite:
    name: str
    hmdb_id: str
    smiles: str
    pathways: Tuple[Pathway]
    enzymes: Tuple[Enzyme]
    cas_registry_number: str = None
