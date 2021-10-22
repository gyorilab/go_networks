import logging
import pickle
from typing import List, Set, Dict

import boto3
import pandas as pd
from tqdm import tqdm

from indra.sources import indra_db_rest
from indra.sources.indra_db_rest import DBQueryStatementProcessor
from indra.statements import *
from indra.util import batch_iter
from indra_db.cli.dump import Sif, S3Path, get_latest_dump_s3_path

__all__ = [
    "load_latest_sif",
    "stmts_by_directedness",
    "set_directed",
    "set_reverse_directed",
    "DIRECTED_TYPES",
    "UNDIRECTED_TYPES",
    "set_pair",
]


logger = logging.getLogger(__name__)


def load_latest_sif() -> pd.DataFrame:
    sif_s3p: S3Path = get_latest_dump_s3_path(Sif.name)
    if sif_s3p is None:
        raise ValueError("No sif file found on S3")

    s3 = boto3.client("s3")
    fileio = sif_s3p.get(s3)
    sif_df = pickle.loads(fileio["Body"].read())
    assert isinstance(sif_df, pd.DataFrame)
    return sif_df


def download_statements(hashes: Set[int]) -> Dict[int, Statement]:
    stmts_by_hash = {}
    for hash_iter in tqdm(batch_iter(hashes, batch_size=1000),
                          total=len(hashes)//1000 + 1):
        idbp: DBQueryStatementProcessor = \
            indra_db_rest.get_statements_by_hash(list(hash_iter), ev_limit=10)
        for stmt in idbp.statements:
            stmts_by_hash[stmt.get_hash()] = stmt
    return stmts_by_hash


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
        # Association
        stmt_types = {Complex.__name__, Association.__name__}
    return sorted(stmt_types)


def set_directed(sif_df: pd.DataFrame):
    """Extend dataframe with column "directed" indicting directedness

    "directed": true if directed A->B statement exists otherwise false

    Parameters
    ----------
    sif_df:
        The dataframe to manipulate
    """

    # Extend with "directed" column: check rows that have directed stmt types
    directed_rows: pd.Series = sif_df.stmt_type.isin(DIRECTED_TYPES)
    undirected_rows: pd.Series = sif_df.stmt_type.isin(UNDIRECTED_TYPES)
    assert directed_rows.sum() + undirected_rows.sum() == sif_df.shape[0]

    # Set new column
    sif_df["directed"] = pd.NA

    # Set directed
    sif_df.loc[directed_rows, "directed"] = True

    # Set undirected
    sif_df.loc[undirected_rows, "directed"] = False

    # Check that no rows were left unset
    assert sif_df.directed.isna().sum() == 0


def set_reverse_directed(sif_df: pd.DataFrame):
    """Set reverse directed property for each pair (A, B)

    "reverse directed": true if directed B->A statement exists otherwise false

    Since each new value per row depends on looking up other rows, it's hard
    to implement some vectorized version, use temporary columns instead.

    Parameters
    ----------
    sif_df:
        DataFrame to set column reverse_directed in
    """
    # Set directed column if it doesn't exist
    if "directed" not in sif_df.columns:
        set_directed(sif_df)

    sif_df["reverse_directed"] = False

    # Set temporary column AB f'{agA_name}_{agB_name}'
    sif_df["AB"] = sif_df.agA_name + "_" + sif_df.agB_name
    # Set temporary column BA f'{agB_name}_{agA_name}'
    sif_df["BA"] = sif_df.agB_name + "_" + sif_df.agA_name

    # Set column reverse_directed:
    # "if each AB in BA & is directed"
    sif_df.loc[
        (sif_df.AB.isin(sif_df.BA.values) & sif_df.directed), "reverse_directed"
    ] = True
    assert (~sif_df.directed & sif_df.reverse_directed).sum() == 0

    # Drop temporary columns
    sif_df.drop(columns=["AB", "BA"], inplace=True)


def set_pair(sif_df: pd.DataFrame):
    """Set pair column in DataFrame

    The pair is constructed for each ordered pair (A, B) such that all
    interactions, directed or undirected, can be grouped together

    Parameters
    ----------
    sif_df :
        DataFrame to set column pair in
    """
    sif_df['pair'] = sif_df.apply(lambda r: f'{r.agA_name}|{r.agb_name}',
                                  axis=1)


# statement types by directedness
DIRECTED_TYPES = stmts_by_directedness(True)
UNDIRECTED_TYPES = stmts_by_directedness(False)
