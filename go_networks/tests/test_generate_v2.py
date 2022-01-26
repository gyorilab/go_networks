from go_networks.tests import _gen_df
from go_networks.util import DIRECTED_TYPES, UNDIRECTED_TYPES, set_directed


def test_set_directed():
    sif = _gen_df()
    directed_count = sif.stmt_type.isin(DIRECTED_TYPES).sum()
    undirected_count = sif.stmt_type.isin(UNDIRECTED_TYPES).sum()
    set_directed(sif)

    assert sif.directed.sum() == directed_count
    assert (sif.directed == False).sum() == undirected_count
    assert sif.directed.isna().sum() == 0
