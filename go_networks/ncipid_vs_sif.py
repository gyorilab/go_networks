"""
Compare  the NCIPID CX network on NDEx with the INDRA Sif dump in order to
find the sets of interactions for each of them and compare intersection and
differences.
"""
import csv
import gzip
import json
import shutil
from contextlib import closing
from pathlib import Path
from urllib import parse, request
import argparse
import logging
import os
import pickle
from typing import Tuple, Optional, List, Dict, Union
from urllib.error import URLError

import numpy as np
import pandas as pd
from ndex2 import create_nice_cx_from_server
from ndex2.client import Ndex2
from ndex2.nice_cx_network import NiceCXNetwork
from lxml import etree
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from tqdm import tqdm

from go_networks.util import get_ndex_web_client, NDEX_ARGS
from indra.statements import Statement
from indra.statements.agent import default_ns_order
from indra.ontology.bio import bio_ontology
from indra.preassembler.grounding_mapper.gilda import get_grounding
from indra.util.statement_presentation import reader_sources, db_sources
from indra.sources.biopax import process_owl

from protmapper.uniprot_client import get_primary_id

from ndex2 import create_nice_cx_from_file

logger = logging.getLogger(__name__)


NDEX_BASE_URL_V2 = "http://public.ndexbio.org/v2"
NETWORK_ENDPOINT = f"{NDEX_BASE_URL_V2}/network"
NETWORKSET_ENDPOINT = f"{NDEX_BASE_URL_V2}/networkset"
NCIPID_SET = "8a2d7ee9-1513-11e9-bb6a-0ac135e8bacf"
NDEX_FTP_BASE = (
    "ftp://ftp.ndexbio.org/NCI_PID_BIOPAX_2016-06-08-PC2v8-API/{pathway_name}.owl.gz"
)


proteins = ['FPLX', 'UPPRO', 'HGNC', 'UP']
small_molecules = ['CHEBI', 'CHEMBL', 'PUBCHEM']


def _rank(list_of_names: List[str], name) -> int:
    return list_of_names.index(name) if name in list_of_names else len(list_of_names)


def _belief_filter(bl: Union[List[float], float], bc: float) -> bool:
    if isinstance(bl, list):
        return False if all(pd.isna(bl)) else any(b > bc for b in bl)
    return False if pd.isna(bl) else bl > bc


def _source_count_filter(
    scl: List[Dict[str, int]], minc: int, sources: List[str]
) -> bool:
    for scd in scl:
        if scd is not None and len(set(scd) & set(sources)) >= minc:
            return True
    return False


def _reader_count_filter(scl: List[Dict[str, int]], minc: int) -> bool:
    return _source_count_filter(scl, minc, reader_sources)


def _db_count_filter(scl: List[Dict[str, int]], minc: int) -> bool:
    return _source_count_filter(scl, minc, db_sources)


def _normalize(db_name: str) -> str:
    # Normalize the ns string to be the same as the INDRA SIF dump

    # Handle uniprot
    if db_name == "uniprot":
        return "UP"
    # Todo: Handle other namespaces as they are added
    return db_name


def _ndex_ftp_owl_url(pathway_name: str) -> str:
    """Get the url for the NDEx FTP owl file for a given pathway."""
    # URL encode the pathway name
    url_encoded = parse.quote(pathway_name)
    return NDEX_FTP_BASE.format(pathway_name=url_encoded)


def _get_ndex_graph_info(network_id: str, client: Optional[Ndex2]) -> Tuple[str, str]:
    """Get name and ftp url for a network

    Parameters
    ----------
    client : ndex2.client.Ndex2
        An NDEx client
    network_id :
        The UUID of the network

    Returns
    -------
    :
        Tuple of name and ftp url
    """
    if client is None:
        client = get_ndex_web_client()

    info_dict = client.get_network_summary(network_id)
    link_tag = None
    for prop_dict in info_dict["properties"]:
        if prop_dict["predicateString"] == "prov:wasDerivedFrom":
            link_tag = prop_dict["value"]
            break
    if link_tag is None:
        raise ValueError("Network has no 'wasDerivedFrom' value")

    # Extract the link from the tag
    ftp_url = etree.fromstring(link_tag).values()[0]

    # Extract the network name
    network_name = info_dict["name"]

    return network_name, ftp_url


def _get_cx_graph_from_server(
    network_id: str, out_dir: Union[Path, str]
) -> NiceCXNetwork:
    """Get a NiceCXNetwork from an NDEx network."""
    # Check if the network is already in the local cache
    cx_file = Path(out_dir).joinpath(f"{network_id}.cx")
    if Path(cx_file).exists():
        return create_nice_cx_from_file(cx_file)

    # Otherwise, download the network from the server and dump it to the local cache
    cx = create_nice_cx_from_server(**NDEX_ARGS, uuid=network_id)
    with cx_file.open("w") as f:
        json.dump(cx.to_cx(), f)
    return cx


def _download_extract_owl_file(ftp_url: str, out_dir: str) -> str:
    """Download the owl file for a network from NDEx. and return the name"""
    owl_gz = ftp_url.split("/")[-1]
    assert owl_gz.endswith("owl.gz"), f"Expected owl.gz file, got {owl_gz}"
    gz_path = Path(out_dir).joinpath(owl_gz)
    if not gz_path.exists():
        # Download the file at the url for the ftp server to the given file path
        with closing(request.urlopen(ftp_url)) as r:
            with gz_path.open("wb") as f:
                shutil.copyfileobj(r, f)

    # Extract the owl file from the gz file
    owl_path = Path(out_dir).joinpath(owl_gz.replace(".gz", ""))
    if not owl_path.exists():
        with gzip.open(gz_path, "rb") as f, owl_path.open("wb") as f2:
            shutil.copyfileobj(f, f2)
    else:
        logger.info(f"{owl_path} already exists, skipping download")

    return owl_path.absolute().as_posix()


def _get_nci_statements(owl_file: str) -> List[Statement]:
    """Get statements from and OWL file using the biopax api"""
    bp = process_owl(owl_file)
    return bp.statements


def _get_grounding(name, original_ns=None):
    """Get the grounding from the gilda"""
    gilda_res = get_grounding(name)
    if gilda_res and gilda_res[1]:  # gilda_res == ({...}, [...])
        # Get first match that also corresponds to the original namespace
        if original_ns is None or original_ns not in small_molecules or \
                original_ns not in proteins:
            term_match = gilda_res[1][0]["term"]
            return term_match["entry_name"], term_match["db"], term_match["id"]

        if original_ns in small_molecules:
            ns_set = small_molecules
        elif original_ns in proteins:
            ns_set = proteins
        else:
            raise ValueError(f"Unknown namespace {original_ns}")

        # Get the top ranked match withing the namespace class: go by the
        # index of the entity in ns_set. If tied, go by the gilda score
        res = []
        for term_match in gilda_res[1]:
            term = term_match["term"]
            term_name = term["entry_name"]
            term_ns = term["db"]
            term_id = term["ib"]
            rank = ns_set.index(term_ns) or len(ns_set)
            score = term_match["score"]
            res.append((rank, score, term_name, term_ns, term_id))

        # Sort by rank, then by score
        return sorted(res, key=lambda x: (x[0], x[1]))[0][2:]


def _get_node_info(
    cx, cx_id, sif_name_map=None, use_gilda: bool = True
) -> Optional[Tuple[str, str, str]]:
    # Shorthand for looking up ns-id pair from name
    def _nim(nom: str) -> Optional[Tuple[str, str]]:
        return sif_name_map.get(nom)

    # Shorthand for getting potential aliases for a node
    def _get_aliases(nid) -> List[str]:
        ad = cx.get_node_attribute(nid, "aliases")
        if ad:
            return ad.get("v", [])
        return []

    # Shorthand for getting the name of a node through bio_ontology.get_name
    def _ont_n(db_name: str, db_id: str, default: str) -> str:
        return bio_ontology.get_name(db_name, db_id) or default

    # Get the node attribute
    meta = cx.get_node(cx_id)
    if not meta:
        logger.warning(f"Node {cx_id} not found in the network")
        return None

    # Normalize the ns-id strings to be the same as the INDRA SIF dump
    name = meta.get("n")
    ns_id = meta.get("r")
    # print(f"id: {cx_id}; name: {name}, ns_id: {ns_id}")  # Debug
    # Skip ungrounded entries
    if ":" not in ns_id:
        return None

    # assume the format is "<db_name>:<id>"
    raw_ns, _id = ns_id.split(":")

    # Normalize the ns string to be the same as the INDRA SIF dump
    ns = _normalize(raw_ns)

    # If ns == uniprot/UP, try to get HGNC id from Sif name map
    if ns == "UP":

        # Get the HGNC id from the SIF name map
        hgnc_id = _nim(name) if sif_name_map else None

        # If we have protein mapping, return it
        if hgnc_id is not None and hgnc_id[0] == "HGNC":
            hns, hid = hgnc_id
            return name, hns, hid

        # Try to get the primary ID from the "grounded" ns/id, then map
        # it to the highest ranked ns/id
        up_prim_id = get_primary_id(_id)
        if up_prim_id:
            xrefs = bio_ontology.get_mappings("UP", up_prim_id)
            for xns, xid in sorted(xrefs, key=lambda x: _rank(default_ns_order, x[0])):
                if xns in default_ns_order:
                    return _ont_n(xns, xid, name), xns, xid

        # If we still don't have a primary ID, try to get it from the aliases
        up_aliases = _get_aliases(cx_id)
        for up_alias in up_aliases:
            up_id = up_alias.split(":")[1]
            up_prim_id = get_primary_id(up_id)
            if up_prim_id:
                xrefs = bio_ontology.get_mappings("UP", up_prim_id)
                for xns, xid in sorted(
                    xrefs, key=lambda x: _rank(default_ns_order, x[0])
                ):
                    if xns in default_ns_order:
                        return _ont_n(xns, xid, name), xns, xid

    # XREFS for non-uniprot namespaces
    else:
        # Check if entity is already valid, return normalized name with the ns-id
        if ns in default_ns_order:
            xrefs = bio_ontology.get_mappings(ns, _id)
            xrefs.append((ns, _id))
            for xns, xid in sorted(xrefs, key=lambda x: _rank(default_ns_order, x[0])):
                if xns in default_ns_order:
                    return _ont_n(xns, xid, name), xns, xid

        # Loop the xrefs and grab the first one that's in default_ns_order
        xrefs = bio_ontology.get_mappings(ns, _id)
        for xns, xid in sorted(xrefs, key=lambda x: _rank(default_ns_order, x[0])):
            if xns in default_ns_order:
                return _ont_n(xns, xid, name), xns, xid

        # If we still don't have a grounding, get it from the aliases
        aliases = _get_aliases(cx_id)
        for alias in aliases:
            ans, aid = alias.split(":")
            xrefs = bio_ontology.get_mappings(ans, aid)
            for xns, xid in sorted(xrefs, key=lambda x: _rank(default_ns_order, x[0])):
                if xns in default_ns_order:
                    return _ont_n(xns, xid, name), xns, xid

    # Last resort: try to get a name from gilda
    if use_gilda:
        grounding = _get_grounding(name, ns)
        if grounding:
            return grounding
        # Try to get a grounding without 'family' in the name. This works for
        # e.g. "Gi" vs "Gi family"
        if "family" in name.lower():
            name = name.replace("family", "").strip()
            grounding = _get_grounding(name, ns)
            # Only return if gilda found a family in FPLX
            if grounding and grounding[1] == "FPLX":
                return grounding


def get_node_mapping(
    cx, sif_name_map=None, use_gilda=True
) -> Dict[str, Tuple[str, str, str]]:
    """Get a mapping from node id to INDRA entity."""

    id_to_entity = {}
    nnodes = len(cx.nodes)
    bio_ontology.initialize()
    logger.info("Mapping nodes to entities")
    for node in tqdm(cx.nodes, total=nnodes):
        node_meta = cx.get_node(node)

        # Get name and ns-id
        entity = _get_node_info(cx, node, sif_name_map, use_gilda)

        # Add to the mapping
        if entity is None:
            entity = node_meta.get("n"), "TEXT", node_meta.get("r")

        id_to_entity[node] = entity

    return id_to_entity


def build_cx_sif(cx, node_id_to_entity) -> pd.DataFrame:
    # Loop the edges and add them to a list and then to a DataFrame
    s_names = []
    s_ns_list = []
    s_id_list = []
    t_names = []
    t_ns_list = []
    t_id_list = []
    int_list = []
    pmid_list = []

    logger.info("Adding edges to the CX SIF input")
    nedges = len(cx.edges)
    for e in tqdm(cx.edges, total=nedges):
        ed = cx.get_edge(e)
        s, interaction, t = ed["s"], ed["i"], ed["t"]

        s_entity = node_id_to_entity.get(s)
        t_entity = node_id_to_entity.get(t)

        if s_entity is None or t_entity is None:
            continue

        citations = cx.get_edge_attribute_value(e, "citation")
        assert isinstance(citations, list)

        s_name, s_ns, s_id = s_entity
        t_name, t_ns, t_id = t_entity
        s_names.append(s_name)
        s_ns_list.append(s_ns)
        s_id_list.append(s_id)
        t_names.append(t_name)
        t_ns_list.append(t_ns)
        t_id_list.append(t_id)
        int_list.append(interaction)
        pmid_list.append(citations)

        # If the interaction is a complex, add the reverse edge as well
        if interaction == "in-complex-with":
            s_names.append(t_name)
            s_ns_list.append(t_ns)
            s_id_list.append(t_id)
            t_names.append(s_name)
            t_ns_list.append(s_ns)
            t_id_list.append(s_id)
            int_list.append(interaction)
            pmid_list.append(citations)

    # Create a DataFrame
    logger.info("Creating the CX DataFrame")
    cx_sif = pd.DataFrame(
        {
            "agA_name": s_names,
            "agA_ns": s_ns_list,
            "agA_id": s_id_list,
            "agB_name": t_names,
            "agB_ns": t_ns_list,
            "agB_id": t_id_list,
            "interaction": int_list,
            "pmids": pmid_list,
        }
    )

    # Drop duplicates and reset the index (must set values as str to avoid
    # error on lists being unhashable)
    cx_sif = cx_sif.loc[cx_sif.astype(str).drop_duplicates().index]
    cx_sif.reset_index(inplace=True, drop=True)
    return cx_sif


def merge_dfs(sif, cx, merge_how="outer") -> pd.DataFrame:
    # Add a new columns to both of the data frames that maps statement
    # type/interaction to boolean indicating if the interaction is directed
    # direction of the interaction
    logger.info("Adding directed column to the DataFrames")
    sif["directed"] = sif.stmt_type != "Complex"
    cx["directed"] = cx.interaction != "in-complex-with"

    groupby_cols = [
        "agA_name",
        "agA_ns",
        "agA_id",
        "agB_name",
        "agB_ns",
        "agB_id",
        "directed",
    ]
    mergeon_cols = ["agA_ns", "agA_id", "agB_ns", "agB_id", "directed"]

    # Group the CX SIF by entity pair A B and directed
    logger.info("Grouping the CX SIF by entity pair")
    cx = (
        cx.groupby(groupby_cols)
        # Aggregate pmids to list of lists
        .aggregate(
            {"interaction": pd.Series.tolist, "pmids": pd.Series.tolist}
        ).reset_index(groupby_cols)
    )

    # Group the INDRA SIF by entity and stmt type
    logger.info("Grouping the INDRA SIF by entity pair")
    sif = (
        sif.groupby(groupby_cols)
        .aggregate(
            {
                "evidence_count": np.sum,
                "stmt_hash": pd.Series.tolist,
                "belief": pd.Series.tolist,
                "source_counts": pd.Series.tolist,
                "stmt_type": pd.Series.tolist,
            }
        )
        .reset_index(groupby_cols)
    )

    # Merge the two DataFrames on A-B pairs + interaction type
    logger.info("Merging the CX SIF with the INDRA SIF")
    merged_df = sif.merge(
        cx,
        on=mergeon_cols,
        how=merge_how,
        suffixes=("_sif", "_cx"),
        indicator=True,
    )
    return merged_df


def get_merged_df(
    sif_file: str,
    ncipid_file: str,
    out_dir: str,
    regenerate: bool = False,
    plot_venn: bool = False,
) -> pd.DataFrame:
    """Get the merged DataFrame from the SIF and CX files."""
    merged_file = Path(out_dir).joinpath("merged_df.pkl")
    if not regenerate and merged_file.exists():
        logger.info(
            f"Loading the merged DataFrame from local cache "
            f"{merged_file.absolute().as_posix()}"
        )
        merged_df = pd.read_pickle(merged_file)

        # Do venn diagram plotting if requested
        if plot_venn:
            venn_plots(merged_df, out_dir)
        return merged_df

    # Load the INDRA SIF dump
    logger.info(f"Loading the INDRA SIF from {sif_file}")
    sif_df = pd.read_pickle(sif_file)

    # Make a name to NS-ID mapping from the sif dump
    sif_ns_id_map = {
        n: (ns, _id)
        for ns, _id, n in set(zip(sif_df.agA_ns, sif_df.agA_id, sif_df.agA_name)).union(
            set(zip(sif_df.agB_ns, sif_df.agB_id, sif_df.agB_name))
        )
    }

    # Load the CX network
    logger.info(f"Loading the CX network from {ncipid_file}")
    nci_cx = create_nice_cx_from_file(ncipid_file)

    # Get node id to entity map
    node_id_to_entity = get_node_mapping(nci_cx, sif_ns_id_map)

    # Build the CX SIF file
    logger.info("Building CX SIF file")
    cx_sif = build_cx_sif(nci_cx, node_id_to_entity)

    # Merge the two data frames
    merged_df = merge_dfs(sif_df, cx_sif, merge_how='outer')
    merged_df.to_pickle(merged_file.absolute().as_posix())

    # Do venn diagram plotting if requested
    if plot_venn:
        venn_plots(merged_df, out_dir)

    return merged_df


def venn_plots(merged_df: pd.DataFrame, out_dir: str):
    """Plot the venn diagrams for the given merged_df"""
    # Gather index sets for venn plotting
    logger.info("Gathering index sets for venn plotting")
    t = tqdm(total=9)
    cx_tot_ix = set(
        merged_df.query(
            "_merge in ['right_only', 'both'] & agA_ns == 'HGNC' & agB_ns == 'HGNC'"
        ).index.to_list()
    )
    t.update()
    sif_merged: pd.DataFrame = merged_df.query(
        "_merge in ['left_only', 'both'] & agA_ns == 'HGNC' & agB_ns == 'HGNC'"
    )
    t.update()
    sif_tot_ix = set(sif_merged.index.to_list())
    t.update()

    # Belief cutoffs
    sif_b7_ix = set(
        sif_merged[
            sif_merged.belief.apply(lambda bl: _belief_filter(bl, 0.7))
        ].index.to_list()
    )
    t.update()

    sif_b9_ix = set(
        sif_merged[
            sif_merged.belief.apply(lambda bl: _belief_filter(bl, 0.9))
        ].index.to_list()
    )
    t.update()

    sif_b99_ix = set(
        sif_merged[
            sif_merged.belief.apply(lambda bl: _belief_filter(bl, 0.99))
        ].index.to_list()
    )
    t.update()

    # Has at least two readers in support (per hash, not combined e.g. not
    # {'sparser': 1} {'reach': 2} combined, but rather:
    # any(len(set(readers) & set(sc)) >= 2
    #     for sc in [{'sparser': 1}, {'reach': 2}])
    sif_min2_readers_ix = set(
        sif_merged[
            sif_merged.source_counts.apply(lambda scl: _reader_count_filter(scl, 2))
        ].index.to_list()
    )
    t.update()

    sif_min3_readers_ix = set(
        sif_merged[
            sif_merged.source_counts.apply(lambda scl: _reader_count_filter(scl, 3))
        ].index.to_list()
    )
    t.update()

    # Has DB support
    sif_has_db_ix = set(
        sif_merged[
            sif_merged.source_counts.apply(lambda scl: _db_count_filter(scl, 1))
        ].index.to_list()
    )
    t.update()
    t.close()

    # Plot the cx sif indices vs each of the sif indices in venn diagrams in
    # subplots
    logger.info("Plotting the venn diagrams")
    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(10, 10))
    axs = axs.flatten()
    for i, (ix, title) in enumerate(
        [
            (sif_tot_ix, "Unfiltered"),
            (sif_b7_ix, "Belief > 0.7"),
            (sif_b9_ix, "Belief > 0.9"),
            (sif_b99_ix, "Belief > 0.99"),
            (sif_min2_readers_ix, "2+ Readers"),
            (sif_min3_readers_ix, "3+ Readers"),
            (sif_has_db_ix, "Has DB Support"),
        ]
    ):
        ax = axs[i]
        venn2(
            subsets=(cx_tot_ix, ix),
            set_labels=("CX SIF", title),
            ax=ax,
        )
        ax.set_title(title)
    fig.tight_layout()
    outpath = os.path.join(out_dir, "venn_plots.png")
    fig.savefig(outpath)
    logger.info(f"Saved venn diagrams to {outpath}")


def get_nci_only_edges(merged_df: pd.DataFrame):
    """Get the missing edges from the merged_df"""
    # Get the interactions that are missing in the INDRA SIF
    nci_only_hgnc = merged_df[["agA_name_cx", "agB_name_cx", "interaction"]][
        (merged_df.agA_ns == "HGNC")
        & (merged_df.agB_ns == "HGNC")
        & (merged_df._merge == "right_only")
    ]

    # Return a list of tuples of the missing edges with their interactions
    return list(
        map(tuple, nci_only_hgnc[["agA_name_cx", "agB_name_cx", "interaction"]].values)
    )


def identify_missing_edges_in_cx_graph(
    cx_graph: NiceCXNetwork, edges
) -> List[Tuple[str, str]]:

    # Get node to name mapping
    node_to_name = {node: cx_graph.get_node(node)["n"] for node in cx_graph.nodes}

    nci_only_edges = set()

    # Check if the edges are in the graph
    for e_id in cx_graph.edges:
        # Get source and target of edge
        ed = cx_graph.get_edge(e_id)
        source, target = ed["s"], ed["t"]
        source_name = node_to_name[source]
        target_name = node_to_name[target]

        # Check if the edge is in the list of edges (check reverse edge as well)
        if (source_name, target_name) in edges:
            nci_only_edges.add((source_name, target_name))
        elif (target_name, source_name) in edges:
            nci_only_edges.add((target_name, source_name))

    return list(sorted(nci_only_edges, key=lambda x: x[0]))


def main(
    sif_file,
    ncipid_file,
    network_set_id,
    out_dir,
    regenerate_merged_df=False,
    plot_venn=False,
):
    """Compare NCI and INDRA Sif, get statements per NCI graph

    Parameters
    ----------
    sif_file :
        The INDRA SIF dump file location
    ncipid_file :
        The nci-pid CX dump file location
    out_dir :
        The output directory
    network_set_id :
        The NCI network set id
    regenerate_merged_df :
        Whether to save the merged DataFrame to a pickle file.
    plot_venn :
        Whether to plot Venn diagrams of the interactions in the merged
        DataFrame.
    """
    # Create the output directory if it doesn't exist
    Path(out_dir).exists() or Path(out_dir).mkdir(parents=True)

    # Get the pairs that are only in the NCI graphs
    nci_only_pairs = get_nci_only_edges(
        get_merged_df(
            sif_file,
            ncipid_file,
            out_dir,
            regenerate=regenerate_merged_df,
            plot_venn=plot_venn,
        )
    )
    client = Ndex2(
        **{(k if k != "server" else "host"): v for k, v in NDEX_ARGS.items()}
    )

    # Get the network ids for the NCI graphs
    nci_graph_ids = _get_networks_in_set(network_set_id, client)

    # For each graph, download and unzip the owl file, get the statements from
    # it and pickle those, save as a csv the nci only edges
    for graph_id in tqdm(nci_graph_ids):
        # Get the name and owl ftp url
        name, owl_url = _get_ndex_graph_info(graph_id, client)
        graph_dir = Path(out_dir).joinpath(name)
        graph_dir.exists() or graph_dir.mkdir()

        # Download and unzip the owl file
        try:
            owl_file_path = _download_extract_owl_file(
                owl_url, graph_dir.absolute().as_posix()
            )
        except URLError as err:
            logger.warning(f"Error downloading {owl_url} for {graph_id}: {err}")
            continue
        # Get the statements from the owl file and pickle them
        statements = _get_nci_statements(owl_file_path)
        if statements:
            pickle_path = graph_dir.joinpath("statements.pkl")
            with pickle_path.open("wb") as f:
                pickle.dump(statements, f)

        # Get the cx graph
        cx = _get_cx_graph_from_server(graph_id, graph_dir)

        # Get the nci only edges
        nci_only_edges = identify_missing_edges_in_cx_graph(cx, nci_only_pairs)

        # Write the nci only edges to a csv
        if nci_only_edges:
            nci_only_edges_path = graph_dir.joinpath("nci_only_edges.csv")
            with nci_only_edges_path.open("w") as f:
                writer = csv.writer(f, delimiter=",")
                writer.writerow(["source", "target"])
                writer.writerows(nci_only_edges)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge the INDRA SIF with the nci-pid CX SIF and find "
        "the sets of interactions in INDRA vs NCIPID"
    )
    parser.add_argument(
        "sif_file",
        help="The INDRA SIF dump file location",
    )
    parser.add_argument(
        "ncipid_file",
        help="The nci-pid CX dump file location",
    )
    parser.add_argument(
        "network_set_id",
        help="The network set id to grab cx and owl files from",
    )
    parser.add_argument(
        "--regenerate-merged-df",
        action="store_true",
        help="Whether to regenerate the merged DataFrame from scratch. If "
        "not specified, the merged DataFrame will be loaded from cache.",
    )
    parser.add_argument(
        "--plot-venn", action="store_true", help="Plot the venn diagrams"
    )
    parser.add_argument(
        "--out-dir",
        help="The output directory for the merged dataframe, the venn "
        "diagrams and the dataframe of graph-interaction comparisons",
        default=".",
    )
    args = parser.parse_args()

    # Run the main function if the dataframe does not already exist
    main(
        args.sif_file,
        args.ncipid_file,
        args.network_set_id,
        args.out_dir,
        regenerate_merged_df=args.regenerate_merged_df,
        plot_venn=args.plot_venn,
    )
