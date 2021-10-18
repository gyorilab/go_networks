"""
BaseModels carrying data and properties
"""
from typing import List, Dict, Tuple

from pydantic import BaseModel, constr

NonEmptyStr = constr(strip_whitespace=True, min_length=1)


class Entity(BaseModel):
    """Data model to set for each entity"""

    ns: NonEmptyStr
    id: NonEmptyStr
    name: NonEmptyStr


class PairProperty(BaseModel):
    """Data model for each pair (A, B) in the sif dump"""

    a: Entity
    b: Entity
    order: Tuple[str, str]  # The original order of A,B
    hashes: Dict[str, List[int]]  # (stmt_type, hash)
    directed: bool  # if directed A->B statement exists
    reverse_directed: bool  # true if directed B->A statements exists
    directed_evidence_count: Dict[str, int]  # ev count per statement type
    reverse_directed_evidence_count: Dict[str, int]  # ev count per statement type
    undirected_evidence_count: Dict[str, int]  # ev count per statement type
