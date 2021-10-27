"""
Generate GO Networks
"""
from collections import defaultdict
import logging
import pickle
import obonet
from itertools import combinations
from pathlib import Path
from typing import Optional, Dict, Set, Tuple, List, Union

import networkx as nx
import numpy as np
import pandas as pd
from networkx import MultiDiGraph
from tqdm import tqdm

from go_networks.data_models import PairProperty, Entity, StmtsByDirectness
from go_networks.util import (
    DIRECTED_TYPES,
    UNDIRECTED_TYPES,
    load_latest_sif,
)
from go_networks.network_assembly import GoNetworkAssembler
from indra.databases import uniprot_client, hgnc_client

# Derived types
Go2Genes = Dict[str, Set[str]]
EvCountDict = Dict[str, Dict[str, int]]
NameEntityMap = Dict[str, Tuple[str, str]]

# Constants
HERE = Path(__file__).absolute().parent.parent
GO_ANNOTS_PATH = HERE.joinpath("goa_human.gaf").absolute().as_posix()
GO_OBO_PATH = HERE.joinpath("go.obo").absolute().as_posix()
LOCAL_SIF = '/Users/ben/.data/indra/db/sif.pkl'
PROPS_FILE = HERE.joinpath("props.pkl").absolute().as_posix()
NETWORKS_FILE = HERE.joinpath("networks.pkl").absolute().as_posix()

min_gene_count = 5
max_gene_count = 200

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


def generate_props(
    sif_file: str, props_file: Optional[str] = None
) -> Dict[str, PairProperty]:
    """Generate properties per pair from the Sif dump

    For each pair of genes (A,B) (excluding self loops), generate the
    following properties:
        - "SOURCE => TARGET": aggregate number of evidences by statement
          type for A->B statements
        - "TARGET => SOURCE": aggregate number of evidences by statement
          type for B->A statements
        - "SOURCE - TARGET": aggregate number of evidences by statement type
          for A-B undirected statements
    """

    if props_file is not None and Path(props_file).is_file():
        logger.info(f"Loading property lookup from {props_file}")
        with Path(props_file).open(mode="rb") as fr:
            props_by_pair = pickle.load(fr)
    else:
        # Load the latest INDRA SIF dump
        sif_df = get_sif(sif_file)

        # Filter to HGNC-only rows
        sif_df = filter_to_hgnc(sif_df)

        # Filter out self-loops
        sif_df = filter_self_loops(sif_df)

        hashes_by_pair = defaultdict(set)
        props_by_hash = {}

        def get_direction(row, pair):
            directed = (row.stmt_type in DIRECTED_TYPES)
            direction = (row.agA_name == pair[0])
            if directed:
                if direction:
                    return "forward"
                else:
                    return "reverse"
            else:
                return "undirected"

        for _, row in tqdm(sif_df.iterrows(), total=len(sif_df)):
            pair = tuple(sorted([row.agA_name, row.agB_name]))
            hashes_by_pair[pair].add(row.stmt_hash)
            if row.stmt_hash not in props_by_hash:
                props_by_hash[row.stmt_hash] = {
                    "ev_count": row.evidence_count,
                    "stmt_type": row.stmt_type,
                    "direction": get_direction(row, pair)
                }
        hashes_by_pair = dict(hashes_by_pair)

        def aggregate_props(props):
            ev_forward = defaultdict(int)
            ev_reverse = defaultdict(int)
            ev_undirected = defaultdict(int)
            for prop in props:
                if prop["direction"] == "forward":
                    ev_forward[prop["stmt_type"]] += prop["ev_count"]
                elif prop["direction"] == "reverse":
                    ev_reverse[prop["stmt_type"]] += prop["ev_count"]
                else:
                    ev_undirected[prop["stmt_type"]] += prop["ev_count"]
            return dict(ev_forward), dict(ev_reverse), dict(ev_undirected)

        props_by_pair = {}
        for pair, hashes in hashes_by_pair.items():
            props_by_pair[pair] = aggregate_props([props_by_hash[h]
                                                   for h in hashes])

        # Write to file if provided
        if props_file:
            logger.info(f"Saving property lookup to {props_file}")
            Path(props_file).absolute().parent.mkdir(exist_ok=True, parents=True)
            with Path(props_file).open(mode="wb") as fo:
                pickle.dump(obj=props_by_pair, file=fo)

    return props_by_pair


def genes_by_go_id():
    """Load the gene/GO annotations as a pandas data frame."""
    go_dag = obonet.read_obo(GO_OBO_PATH)

    goa = pd.read_csv(GO_ANNOTS_PATH, sep='\t',
                      comment='!', dtype=str,
                      header=None,
                      names=['DB',
                             'DB_ID',
                             'DB_Symbol',
                             'Qualifier',
                             'GO_ID',
                             'DB_Reference',
                             'Evidence_Code',
                             'With_From',
                             'Aspect',
                             'DB_Object_Name',
                             'DB_Object_Synonym',
                             'DB_Object_Type',
                             'Taxon',
                             'Date',
                             'Assigned',
                             'Annotation_Extension',
                             'Gene_Product_Form_ID'])
    # Filter out all "NOT" negative evidences
    goa['Qualifier'].fillna('', inplace=True)
    goa = goa[~goa['Qualifier'].str.startswith('NOT')]
    # We can filter to just GO terms in the ontology since
    # obsolete terms are not included in the GO DAG
    goa = goa[goa['GO_ID'].isin(go_dag)]

    genes_by_go_id = defaultdict(set)
    for go_id, up_id in zip(goa.GO_ID, goa.DB_ID):
        if go_dag.nodes[go_id]['namespace'] != 'biological_process':
            continue
        hgnc_id = uniprot_client.get_hgnc_id(up_id)
        if hgnc_id:
            gene_name = hgnc_client.get_hgnc_name(hgnc_id)
            genes_by_go_id[go_id] = genes_by_go_id[go_id] | {gene_name}

    for go_id in set(genes_by_go_id):
        for child_go_id in nx.ancestors(go_dag, go_id):
            genes_by_go_id[go_id] |= genes_by_go_id[child_go_id]

    return genes_by_go_id


def build_networks(go2genes_map: Go2Genes,
                   pair_props: Dict[Tuple, List[Dict[str, int]]],
                   ) -> Dict[str, GoNetworkAssembler]:
    """Build networks per go-id associated genes

    Parameters
    ----------
    go2genes_map :
        A dict mapping GO IDs to lists of genes
    pair_props :
        Lookup for edges

    Returns
    -------
    :
        Dict of assembled networks by go id
    """
    networks = {}
    # Only pass the relevant parts of the pair_props dict
    for go_id, gene_set in tqdm(go2genes_map.items(), total=len(go2genes_map)):
        def _key(g1, g2):
            return tuple(sorted([g1, g2]))
        # Get relevant pairs from pair_properties
        prop_dict = {_key(g1, g2): pair_props[_key(g1, g2)]
                     for g1, g2 in combinations(gene_set, 2)
                     if _key(g1, g2) in pair_props}

        if not prop_dict:
            logger.info(f"No statements for ID {go_id}")
            continue

        gna = GoNetworkAssembler(
            identifier=go_id,
            entity_list=list(gene_set),
            pair_properties=prop_dict,
        )
        gna.assemble()
        networks[go_id] = gna
    return networks


def filter_self_loops(df):
    """Remove self-loops from a dataframe

    Parameters
    ----------
    df :
        The dataframe to filter

    Returns
    -------
    :
        The filtered dataframe
    """
    return df[df.agA_name != df.agB_name]


def filter_go_ids(go2genes_map):
    return {go_id: genes for go_id, genes in go2genes_map.items()
            if min_gene_count <= len(genes) <= max_gene_count}


def generate(sif_file: Optional[str] = None, props_file: Optional[str] = None):
    """Generate new GO networks from INDRA statements

    Parameters
    ----------
    local_sif :
        If provided, load sif dump from this file. Default: load from S3.
    props_file :
        If provided, load property lookup from this file. Default: generate
        from sif dump.
    """
    # Make genes by GO ID dict
    go2genes_map = genes_by_go_id()

    # Filter GO IDs
    go2genes_map = filter_go_ids(go2genes_map)

    # Generate properties
    sif_props = generate_props(sif_file, props_file)

    # Iterate by GO ID and for each list of genes, build a network
    return build_networks(go2genes_map, sif_props)


if __name__ == "__main__":
    # Todo: allow local file to be passed
    networks = generate(LOCAL_SIF, PROPS_FILE)
    with open(NETWORKS_FILE, 'wb') as f:
        pickle.dump(networks, f)
