import logging
from pathlib import Path
from typing import List, Dict, Optional

from ndex2 import Ndex2

from indra.databases import ndex_client
from indra.statements import *

__all__ = [
    "stmts_by_directedness",
    "DIRECTED_TYPES",
    "UNDIRECTED_TYPES",
    "NDEX_ARGS",
    "get_ndex_web_client",
    "get_networks_in_set",
]


logger = logging.getLogger(__name__)

HERE = Path(__file__).absolute().parent


def stmts_by_directedness(directed: bool) -> List[str]:
    """Get a list of statement types that are directed or undirected

    Parameters
    ----------
    directed:
        If True, get directed statement types, otherwise get undirected types

    Returns
    -------
    :
        A list of statement type names
    """
    if directed:
        # Modifications: enz-sub
        # RegulateActivitity: sub-obj
        # Gef: gef-ras
        # Gap: gap-ras
        # RegulateAmount: subj-obj
        stmt_types = {
            Conversion.__name__,
            Modification.__name__,
            Gef.__name__,
            Gap.__name__,
        }
        stmt_types.update(s.__name__ for s in get_all_descendants(Modification))
        stmt_types.update(s.__name__ for s in get_all_descendants(RegulateActivity))
        stmt_types.update(s.__name__ for s in get_all_descendants(RegulateAmount))
    else:
        # Complex
        stmt_types = {Complex.__name__}
    return sorted(stmt_types)


# statement types by directedness
DIRECTED_TYPES = stmts_by_directedness(True)
UNDIRECTED_TYPES = stmts_by_directedness(False)


NDEX_BASE_URL = "http://public.ndexbio.org"


def _get_default_ndex_args() -> Dict[str, str]:
    usr, pwd = ndex_client.get_default_ndex_cred(ndex_cred=None)
    return {"username": usr, "password": pwd, "server": NDEX_BASE_URL}


NDEX_ARGS = _get_default_ndex_args()


def get_ndex_web_client() -> Ndex2:
    return Ndex2(**{(k if k != "server" else "host"): v for k, v in NDEX_ARGS.items()})


def get_networks_in_set(
    network_set_id: str, client: Optional[Ndex2] = None
) -> List[str]:
    """Get the UUIDs of the networks in a network set."""
    if client is None:
        client = get_ndex_web_client()
    network_set = client.get_networkset(network_set_id)
    return network_set["networks"]
