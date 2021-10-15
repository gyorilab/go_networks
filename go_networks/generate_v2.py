"""
Generate GO Networks
"""
import logging
import pickle
from pathlib import Path
from typing import Optional, Dict, Any, Set

import numpy as np
import pandas as pd
from tqdm import tqdm

from go_networks.util import load_latest_sif, set_directed, \
    set_reverse_directed
# Derived types
from indra.databases import uniprot_client

Go2Genes = Dict[str, Set[str]]

# Constants
GO_PATH = Path(__file__).absolute().parent.joinpath("goa_human.gaf")


logger = logging.getLogger(__name__)


def get_sif(local_sif: Optional[str] = None) -> pd.DataFrame:
    if local_sif and Path(local_sif).is_file():
        with open(local_sif, "rb") as bf:
            sif_df = pickle.load(file=bf)

    else:
        sif_df = load_latest_sif()
    assert isinstance(sif_df, pd.DataFrame)
    return sif_df


def filter_to_hgnc(sif: pd.DataFrame) -> pd.DataFrame:
    """Filter sif dataframe to HGNC pairs only"""
    return sif.query("agA_ns == 'HGNC' & agB_ns == 'HGNC'")


def generate_props(sif_df: pd.DataFrame) -> Dict[str, Any]:
    """Generate properties per pair

    For each pair of genes (A,B) (excluding self loops), generate the
    following properties:
        - "directed": true if directed A->B statement exists otherwise false
        - "reverse directed": true if directed B->A statement exists
          otherwise false
        - "SOURCE => TARGET": aggregate number of evidences by statement
          type for A->B statements
        - "TARGET => SOURCE": aggregate number of evidences by statement
          type for B->A statements
        - "SOURCE - TARGET": aggregate number of evidences by statement type
          for A-B undirected statements
    """
    # Set directed/undirected column
    set_directed(sif_df)

    # Set reverse directed column
    set_reverse_directed(sif_df)

    # Do group-by on stmt_type and get:
    #   - A->B aggregated evidence counts
    #   - B->A aggregated evidence counts
    #   - A-B (undirected) aggregated evidence counts

    dir_count = (
        sif_df[sif_df.directed].groupby("stmt_type").aggregate(np.sum).evidence_count
    ).to_dict()
    rev_dir_count = (
        sif_df[sif_df.reverse_directed]
        .groupby("stmt_type")
        .aggregate(np.sum)
        .evidence_count
    ).to_dict()
    undir_count = (
        sif_df[sif_df.directed == False]
        .groupby("stmt_type")
        .aggregate(np.sum)
        .evidence_count
    ).to_dict()

    pair_props = {
        (A, B): {"directed": d, "reverse_directed": r}
        for A, B, d, r in sif_df[
            ["agA_name", "agB_name", "directed", "reverse_directed"]
        ].values
    }

    return {
        "dir_count": dir_count,
        "rev_dir_count": rev_dir_count,
        "undir_count": undir_count,
        "pair_props": pair_props,
    }


def genes_by_go_id(go_path: str = GO_PATH) -> Go2Genes:
    """For each GO id, get the associated genes

    Parameters
    ----------
    go_path :
        If provided, load the go file from here, otherwise assume the file
        is in the directory of this file with file name goa_human.gaf

    Returns
    -------
    :
        A mapping of GO ids to genes
    """
    goa_df = pd.read_csv(
        go_path,
        sep="\t",
        comment="!",
        dtype=str,
        header=None,
        names=[
            "DB",
            "DB_ID",
            "DB_Symbol",
            "Qualifier",
            "GO_ID",
            "DB_Reference",
            "Evidence_Code",
            "With_From",
            "Aspect",
            "DB_Object_Name",
            "DB_Object_Synonym",
            "DB_Object_Type",
            "Taxon",
            "Date",
            "Assigned",
            "Annotation_Extension",
            "Gene_Product_Form_ID",
        ],
    )
    # Filter out negative evidence
    goa_df = goa_df[goa_df.Qualifier.str.startswith("NOT")]
    goa_df["entity"] = list(zip(goa_df.DB, goa_df.DB_ID, goa_df.DB_Symbol))
    up_mapping = goa_df.groupby("GO_ID").agg({"DB_ID": lambda x: x.tolist()}).to_dict()
    logger.info("Translating genes from UP to HGNC")
    mapping = {}
    for go_id, gene_list in tqdm(up_mapping.items()):
        gene_names = set()
        for up_id in gene_list:
            gene_name = uniprot_client.get_gene_name(up_id)
            if gene_name is not None:
                gene_names.add(gene_name)
        if gene_names:
            mapping[go_id] = gene_names

    return mapping


def build_networks(go2genes_map: Go2Genes, pair_props):
    """Build networks per gene

    Iterate by GO ID and for each list of genes, build a network:
        - Make a node for each gene with metadata representing db_refs
        - Add all the edges between the list of genes with metadata coming
          from the pre-calculated data in (4), i.e. the pair properties
            - When generating links to different statement types, use query
              links instead of hash-based links e.g., instead of
              https://db.indra.bio/statements/from_hash/-26724098956735262?format=html,
              link to
              https://db.indra.bio/statements/from_agents?subject=HBP1&object=CDKN2A&type=IncreaseAmount&format=html
        - Should we include any sentences?
            - Maybe choose the statement involving A and B with the highest
              evidence count / belief score and choose one example sentence
              to add to the edge?
        - How to generate the edge label (the general pattern is
          "A (interaction type) B")?
            - Make these as specific as possible, the least specific being
              "A (interacts with) B"

    go2genes_map:
        A dict mapping GO IDs to lists of genes
    pair_props:

    """
    # Apply Layout after creation
    pass


def generate(local_sif: Optional[str] = None):
    """Generate new GO networks from INDRA statements"""
    # Load the latest INDRA SIF dump
    sif_df = get_sif(local_sif)

    # Filter to HGNC-only rows
    sif_df = filter_to_hgnc(sif_df)

    # Generate properties
    sif_props = generate_props(sif_df)

    # Make genes by GO ID dict
    go2genes_map = genes_by_go_id()

    # Iterate by GO ID and for each list of genes, build a network
    build_networks(go2genes_map, sif_props)


if __name__ == "__main__":
    # Todo: allow local file to be passed
    generate()
