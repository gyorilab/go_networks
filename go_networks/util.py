import logging
import pickle
from typing import List

import boto3
import pandas as pd
from tqdm import tqdm

from indra.statements import *
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

    The pair is constructed for each unordered pair (A, B) such that all
    interactions A->B, B->A, A-B can be grouped together

    Parameters
    ----------
    sif_df :
        DataFrame to set column pair in
    """
    # Created list of pair that has the same length as number of rows in
    # DataFrame, allowing
    pairs = []
    added = set()
    for a_name, b_name in tqdm(
        zip(sif_df.agA_name, sif_df.agB_name), total=sif_df.shape[0]
    ):
        pair = f"{a_name}|{b_name}"
        rev_pair = f"{b_name}|{a_name}"

        # If rev_pair already exists that means (B,A) appeared previously as
        # (A,B) and we don't need to add pair, as rev_pair will cover the
        # same group. Otherwise add pair.
        if rev_pair in added:
            pairs.append(rev_pair)
        else:
            pairs.append(pair)
            added.add(pair)

    assert len(pairs) == sif_df.shape[0]

    sif_df["pair"] = pairs


# statement types by directedness
DIRECTED_TYPES = stmts_by_directedness(True)
UNDIRECTED_TYPES = stmts_by_directedness(False)
