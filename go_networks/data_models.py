"""
BaseModels carrying data and properties
"""
from typing import Dict, Tuple, DefaultDict, List

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
    undirected: DefaultDict[int, list[Statement]]

    def has_data(self):
        """Shortcut to check if model contains any data"""
        return len(self.undirected) > 0 or len(self.directed) > 0

    class Config:
        arbitrary_types_allowed = True


class PairProperty(BaseModel):
    """Data model for each pair (A, B) in the sif dump"""

    a: Entity
    b: Entity
    statements: StmtsByDirectness
    directed: bool  # if directed A->B statement exists
    reverse_directed: bool  # true if directed B->A statements exists
    directed_evidence_count: Dict[str, int]  # ev count per statement type
    reverse_directed_evidence_count: Dict[str, int]  # ev count per statement type
    undirected_evidence_count: Dict[str, int]  # ev count per statement type

    def get_rev_pair(self) -> Tuple[str, str]:
        """Get the reverse pair"""
        return self.b.name, self.a.name
