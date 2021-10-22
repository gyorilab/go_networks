"""
BaseModels carrying data and properties
"""
from typing import List, Dict, Tuple

from pydantic import BaseModel, constr

from indra.statements import Statement


# Derived types
NonEmptyStr = constr(strip_whitespace=True, min_length=1)


class Entity(BaseModel):
    """Data model to set for each entity"""

    ns: NonEmptyStr
    id: NonEmptyStr
    name: NonEmptyStr


class StmtsByDirectness(BaseModel):
    """Group statements based on of they are directed or not"""

    directed: Dict[int, Statement]
    undirected: Dict[int, List[Statement]]

    def is_empty(self):
        """Shortcut to check if model contains any data"""
        return len(self.directed) == 0 and len(self.undirected) == 0


class PairProperty(BaseModel):
    """Data model for each pair (A, B) in the sif dump"""

    a: Entity
    b: Entity
    statements: Dict[int, List[Statement]]
    directed: bool  # if directed A->B statement exists
    reverse_directed: bool  # true if directed B->A statements exists
    directed_evidence_count: Dict[str, int]  # ev count per statement type
    reverse_directed_evidence_count: Dict[str, int]  # ev count per statement type
    undirected_evidence_count: Dict[str, int]  # ev count per statement type

    def get_rev_pair(self) -> Tuple[str, str]:
        """Get the reverse pair"""
        return self.b.name, self.a.name
