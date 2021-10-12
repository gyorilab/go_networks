import pickle

import boto3
import pandas as pd

from indra_db.cli.dump import Sif, S3Path, get_latest_dump_s3_path

__all__ = ["load_lastest_sif"]


def load_lastest_sif() -> pd.DataFrame:
    sif_s3p: S3Path = get_latest_dump_s3_path(Sif.name)
    if sif_s3p is None:
        raise ValueError("No sif file found on S3")

    s3 = boto3.client("s3")
    fileio = sif_s3p.get(s3)
    sif_df = pickle.loads(fileio['Body'].read())
    assert isinstance(sif_df, pd.DataFrame)
    return sif_df
