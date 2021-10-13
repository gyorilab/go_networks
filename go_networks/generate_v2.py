"""
Generate GO Networks
"""
import pickle
from pathlib import Path
from typing import Optional, Dict, List

import pandas as pd

from go_networks.util import load_latest_sif, set_directed, \
    set_reverse_directed

# Derived types
Go2Genes = Dict[str, List[str]]


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


def generate_props(sif_df: pd.DataFrame):
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
        - "SOURCE -TARGET": aggregate number of evidences by statement type
          for A-B undirected statements
    """
    # Set directed/undirected column
    set_directed(sif_df)

    # Set reverse directed column
    set_reverse_directed(sif_df)

    # Do group-by statement types and get counts


def genes_by_go_id() -> Go2Genes:
    """For each GO id get the associated genes"""
    pass


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
