import logging
import pickle
from collections import defaultdict
from copy import deepcopy
from itertools import combinations
from pathlib import Path
from typing import List, Set, Dict

import boto3
import pandas as pd
from tqdm import tqdm

from go_networks.data_models import StmtsByDirectness
from indra.sources import indra_db_rest
from indra.sources.indra_db_rest import DBQueryStatementProcessor
from indra.statements import *
from indra.util import batch_iter
from indra_db.cli.dump import Sif, S3Path, get_latest_dump_s3_path

__all__ = [
    "load_latest_sif",
    "stmts_by_directedness",
    "DIRECTED_TYPES",
    "UNDIRECTED_TYPES",
    #"set_directed",
    #"set_reverse_directed",
    #"set_pair",
    #"get_stmts",
]


logger = logging.getLogger(__name__)

HERE = Path(__file__).absolute().parent


def load_latest_sif() -> pd.DataFrame:
    sif_s3p: S3Path = get_latest_dump_s3_path(Sif.name)
    if sif_s3p is None:
        raise ValueError("No sif file found on S3")

    s3 = boto3.client("s3")
    fileio = sif_s3p.get(s3)
    sif_df = pickle.loads(fileio["Body"].read())
    assert isinstance(sif_df, pd.DataFrame)
    return sif_df


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

# def set_directed(sif_df: pd.DataFrame):
#     """Extend dataframe with column "directed" indicting directedness

#     "directed": true if directed A->B statement exists otherwise false

#     Parameters
#     ----------
#     sif_df:
#         The dataframe to manipulate
#     """

#     # Extend with "directed" column: check rows that have directed stmt types
#     directed_rows: pd.Series = sif_df.stmt_type.isin(DIRECTED_TYPES)
#     undirected_rows: pd.Series = sif_df.stmt_type.isin(UNDIRECTED_TYPES)
#     assert directed_rows.sum() + undirected_rows.sum() == sif_df.shape[0]

#     # Set new column
#     sif_df["directed"] = pd.NA

#     # Set directed
#     sif_df.loc[directed_rows, "directed"] = True

#     # Set undirected
#     sif_df.loc[undirected_rows, "directed"] = False

#     # Check that no rows were left unset
#     assert sif_df.directed.isna().sum() == 0


# def set_reverse_directed(sif_df: pd.DataFrame):
#     """Set reverse directed property for each pair (A, B)

#     "reverse directed": true if directed B->A statement exists otherwise false

#     Since each new value per row depends on looking up other rows, it's hard
#     to implement some vectorized version, use temporary columns instead.

#     Parameters
#     ----------
#     sif_df:
#         DataFrame to set column reverse_directed in
#     """
#     # Set directed column if it doesn't exist
#     if "directed" not in sif_df.columns:
#         set_directed(sif_df)

#     sif_df["reverse_directed"] = False

#     # Set temporary column AB f'{agA_name}_{agB_name}'
#     sif_df["AB"] = sif_df.agA_name + "_" + sif_df.agB_name
#     # Set temporary column BA f'{agB_name}_{agA_name}'
#     sif_df["BA"] = sif_df.agB_name + "_" + sif_df.agA_name

#     # Set column reverse_directed:
#     # "if each AB in BA & is directed"
#     sif_df.loc[
#         (sif_df.AB.isin(sif_df.BA.values) & sif_df.directed), "reverse_directed"
#     ] = True
#     assert (~sif_df.directed & sif_df.reverse_directed).sum() == 0

#     # Drop temporary columns
#     sif_df.drop(columns=["AB", "BA"], inplace=True)


# def set_pair(sif_df: pd.DataFrame):
#     """Set pair column in DataFrame

#     The pair is constructed for each ordered pair (A, B) such that all
#     interactions, directed or undirected, can be grouped together

#     Parameters
#     ----------
#     sif_df :
#         DataFrame to set column pair in
#     """
#     def pair(r):
#         return '|'.join(sorted([r.agA_name, r.agB_name]))
#     sif_df["pair"] = sif_df.apply(pair, axis=1)


# STMTS_PKL: str = HERE.joinpath("_stmts_cache.pkl").absolute().as_posix()
# STMT_CACHE: Dict[int, Statement] = {}

# def download_statements(hashes: Set[int]) -> Dict[int, Statement]:
#     """Get statements from a set of hashes

#     Parameters
#     ----------
#     hashes :
#         The set of hashes for the statements to download

#     Returns
#     -------
#     :
#         A dictionary mapping hash to statement
#     """
#     stmts_by_hash = {}
#     logger.info(f"Downloading {len(hashes)} statements")
#     for hash_iter in tqdm(
#         batch_iter(hashes, batch_size=1000), total=len(hashes) // 1000 + 1
#     ):
#         idbp: DBQueryStatementProcessor = indra_db_rest.get_statements_by_hash(
#             list(hash_iter), ev_limit=10
#         )
#         for stmt in idbp.statements:
#             stmts_by_hash[stmt.get_hash()] = stmt
#     return stmts_by_hash


# def expand_complex(complex_stmt: Complex) -> List[Complex]:
#     """Replace a Complex statement with a list of binary Complex statements

#     Parameters
#     ----------
#     complex_stmt :
#         A Complex statement to expand out to a list of statements

#     Returns
#     -------
#     :
#         A list of complexes with only two memebers
#     """
#     stmts = []
#     added = set()
#     for m1, m2 in combinations(complex_stmt.members, 2):
#         keys = (m1.entity_matches_key(), m2.entity_matches_key())
#         if keys in added:
#             continue
#         if len(set(keys)) == 1:
#             continue
#         ordered = sorted([m1, m2], key=lambda x: x.entity_matches_key())
#         c = Complex(ordered, evidence=deepcopy(complex_stmt.evidence))
#         stmts.append(c)
#         added.add(keys)
#     return stmts


# def get_stmts(sif_df: pd.DataFrame) -> Dict[str, StmtsByDirectness]:
#     # Get hashes by pair
#     hashes_by_pair: Dict[str, List[int]] = list(
#         sif_df.groupby("pair")
#         .aggregate({"stmt_hash": lambda x: x.tolist()})
#         .to_dict()
#         .values()
#     )[0]

#     hashes = set(sif_df.stmt_hash)

#     # Load cached statements
#     if Path(STMTS_PKL).is_file():
#         logger.info("Loading pickled statements")
#         with Path(STMTS_PKL).open('rb') as f:
#             pkl_stmts: Dict[int, Statement] = pickle.load(f)
#         STMT_CACHE.update(pkl_stmts)
#     else:
#         pkl_stmts = {}

#     # Get current stmts - check loaded stmts, then local cache
#     stmt_by_hash = {h: st for h, st in pkl_stmts.items() if h in hashes}
#     for h in hashes:
#         if h in STMT_CACHE and h not in stmt_by_hash:
#             stmt_by_hash[h] = STMT_CACHE[h]

#     # stmt_by_hash.update({h: st for h, st in STMT_CACHE.items() if h in hashes})
#     missing_hashes = hashes.difference(set(stmt_by_hash))

#     if missing_hashes:
#         # Download missing statements
#         logger.info(f"Downloading {len(missing_hashes)} missing statements")
#         new_stmts_by_hash = download_statements(missing_hashes)

#         # Merge dicts
#         stmt_by_hash.update(new_stmts_by_hash)

#         # Update global
#         STMT_CACHE.update(stmt_by_hash)

#         # Check if new statements were added overall and write to cache
#         if set(STMT_CACHE.keys()).difference(pkl_stmts.keys()):
#             logger.info(f"Dumping new statements cache to {STMTS_PKL}")
#             with Path(STMTS_PKL).open('wb') as f:
#                 pickle.dump(obj=STMT_CACHE, file=f)

#     # Combine statements with pairs, expand Complexes
#     stmts_by_pair = {}
#     for pair, hash_list in tqdm(hashes_by_pair.items()):
#         stmt_by_dir = StmtsByDirectness(directed={}, undirected=defaultdict(list))
#         for h, st in zip(hash_list, (stmt_by_hash.get(h) for h in hash_list)):
#             if st is None:
#                 logger.warning(f"No statement for hash {h}")
#                 continue

#             if isinstance(st, Complex):
#                 cplx_stmts = expand_complex(st)
#                 stmt_by_dir.undirected[h].extend(cplx_stmts)
#             else:
#                 stmt_by_dir.directed[h] = st
#         if stmt_by_dir.has_data():
#             stmts_by_pair[pair] = stmt_by_dir
#     return stmts_by_pair
