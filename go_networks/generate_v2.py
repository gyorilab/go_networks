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
    #set_directed,
    #set_reverse_directed,
    #set_pair,
    #get_stmts,
)
from go_networks.network_assembly import GoNetworkAssembler
from indra.databases import uniprot_client

# Derived types
Go2Genes = Dict[str, Set[str]]
EvCountDict = Dict[str, Dict[str, int]]
NameEntityMap = Dict[str, Tuple[str, str]]

# Constants
HERE = Path(__file__).absolute().parent
GO_PATH = HERE.joinpath("goa_human.gaf").absolute().as_posix()
GO_OBO_PATH = HERE.joinpath("go.obo").absolute().as_posix()

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


# def get_pair_properties(
#     dir_dict: Dict[str, Dict[str, bool]],
#     dir_ev_count: EvCountDict,
#     rev_dir_ev_count: EvCountDict,
#     undir_ev_count: EvCountDict,
#     entity_dict: NameEntityMap,
# ) -> Dict[str, PairProperty]:
#     pair_properties = {}
#     for pair, stmts_by_dir in stmts_by_pair.items():
#         # Get properties per pair
#         is_dir = dir_dict[pair]["directed"]
#         is_rev_dir = dir_dict[pair]["reverse_directed"]
#         dir_ec = dir_ev_count.get(pair, {})  # Empty dict = no counts
#         r_dir_ec = rev_dir_ev_count.get(pair, {})  # Empty dict = no counts
#         u_dir_ec = undir_ev_count.get(pair, {})  # Empty dict = no counts

#         # Get name
#         a_name, b_name = pair.split("|")
#         a_ns, a_id = entity_dict[a_name]
#         b_ns, b_id = entity_dict[b_name]

#         # Set entity data
#         a = Entity(ns=a_ns, id=a_id, name=a_name)
#         b = Entity(ns=b_ns, id=b_id, name=b_name)

#         # Add to output dict
#         pair_properties[pair] = PairProperty(
#             a=a,
#             b=b,
#             statements=stmts_by_dir,
#             directed=is_dir,
#             reverse_directed=is_rev_dir,
#             directed_evidence_count=dir_ec,
#             reverse_directed_evidence_count=r_dir_ec,
#             undirected_evidence_count=u_dir_ec,
#         )

#     return pair_properties


def generate_props(
    sif_df: pd.DataFrame, props_file: Optional[str] = None
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

        for _, row in tqdm.tqdm(sif_df.iterrows(), total=len(sif_df)):
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


def genes_by_go_id(
    go_path: str = GO_PATH, go_obo_dag: Optional[MultiDiGraph] = None
) -> Go2Genes:
    """For each GO id, get the associated genes

    Parameters
    ----------
    go_path :
        If provided, load the go file from here, otherwise assume the file
        is in the directory of this file with file name goa_human.gaf
    go_obo_dag :
        If provided, the dag representing the GO ontology hierarchy. If not
        provided, load it from its default

    Returns
    -------
    :
        A mapping of GO ids to genes
    """
    logger.info(f"Reading GO annotations file from {go_path}")
    goa_df: pd.DataFrame = pd.read_csv(
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

    # Get go dag
    if go_obo_dag is None:
        logger.info(f"Reading GO OBO file from {GO_OBO_PATH}")
        go_dag = obonet.read_obo(GO_OBO_PATH)
    else:
        go_dag = go_obo_dag

    # Filter out negative evidence
    goa_df = goa_df[goa_df.Qualifier.str.startswith("NOT")]
    goa_df["entity"] = list(zip(goa_df.DB, goa_df.DB_ID, goa_df.DB_Symbol))
    up_mapping = dict(
        *(
            goa_df.groupby("GO_ID").agg({"DB_ID": lambda x: x.tolist()}).to_dict()
        ).values()
    )
    logger.info("Translating genes from UP to HGNC")
    mapping = {}
    for go_id, gene_list in tqdm(up_mapping.items(), total=len(up_mapping)):
        gene_names = set()
        for up_id in gene_list:
            gene_name = uniprot_client.get_gene_name(up_id)
            if gene_name is not None:
                gene_names.add(gene_name)
        if gene_names:
            mapping[go_id] = gene_names

    logger.info("Extending go to gene mapping with child ids")
    for go_id, gene_set in tqdm(mapping.items(), total=len(mapping)):
        for child_go_id in nx.ancestors(go_dag, go_id):
            gene_set.update(mapping.get(child_go_id, set()))

    return mapping


def build_networks(go2genes_map: Go2Genes, pair_props: Dict[str,
                                                            PairProperty],
                   go_dag) -> Dict[str, GoNetworkAssembler]:
    """Build networks per go-id associated genes

    Parameters
    ----------
    go2genes_map :
        A dict mapping GO IDs to lists of genes
    pair_props :
        Lookup for edges
    go_dag :
        The ontology hierarchy represented in a graph

    Returns
    -------
    :
        Dict of assembled networks by go id
    """
    networks = {}
    # Only pass the relevant parts of the pair_props dict
    for go_id, gene_set in tqdm(go2genes_map.items(), total=len(go2genes_map)):
        # Get relevant pairs from pair_properties
        prop_dict: Dict[str, List[Dict[str, int]]] = {}
        for g1, g2 in combinations(gene_set, 2):
            # Get pair and property for it
            pair = sorted([g1, g2])
            prop = pair_props.get(pair)
            if prop is not None:
                prop_dict[pair] = prop

        if not prop_dict:
            logger.info(f"No statements for ID {go_id}")
            continue

        gna = GoNetworkAssembler(
            identifier=go_id,
            identifier_description=go_dag.nodes[go_id]["name"],
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


def generate(local_sif: Optional[str] = None, props_file: Optional[str] = None):
    """Generate new GO networks from INDRA statements

    Parameters
    ----------
    local_sif :
        If provided, load sif dump from this file. Default: load from S3.
    props_file :
        If provided, load property lookup from this file. Default: generate
        from sif dump.
    """
    # Load the latest INDRA SIF dump
    sif_df = get_sif(local_sif)

    # Filter to HGNC-only rows
    sif_df = filter_to_hgnc(sif_df)

    # Filter out self-loops
    sif_df = filter_self_loops(sif_df)

    # Generate properties
    sif_props = generate_props(sif_df, props_file)

    # Make genes by GO ID dict
    go_dag = obonet.read_obo(GO_OBO_PATH)
    go2genes_map = genes_by_go_id(go_path=GO_PATH, go_obo_dag=go_dag)

    # Iterate by GO ID and for each list of genes, build a network
    return build_networks(go2genes_map, sif_props, go_dag)


if __name__ == "__main__":
    # Todo: allow local file to be passed
    generate()
