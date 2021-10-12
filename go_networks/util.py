import pickle
from typing import List

import boto3
import pandas as pd

from indra.statements import *
from indra_db.cli.dump import Sif, S3Path, get_latest_dump_s3_path

__all__ = ["load_lastest_sif"]


def load_lastest_sif() -> pd.DataFrame:
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
