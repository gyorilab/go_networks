"""
Load the nci-pid CX network and convert it to a SIF file so be compared with
the INDRA SIF dump.
"""
import shutil
from contextlib import closing
from urllib import parse, request
import argparse
import logging
import os
import pickle
from collections import Counter
from typing import Tuple, Optional, List, Dict, Union

import ndex2.client
import numpy as np
import pandas as pd
from lxml import etree
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from tqdm import tqdm

from indra.databases import ndex_client
from indra.statements.agent import default_ns_order
from indra.ontology.bio import bio_ontology
from indra.preassembler.grounding_mapper.gilda import get_grounding
from indra.util.statement_presentation import reader_sources, db_sources

from protmapper.uniprot_client import get_primary_id

from ndex2 import create_nice_cx_from_file

logger = logging.getLogger(__name__)


NDEX_BASE_URL_V2 = "http://public.ndexbio.org/v2"
NDEX_BASE_URL = "http://public.ndexbio.org"
NETWORK_ENDPOINT = f"{NDEX_BASE_URL_V2}/network"
NETWORKSET_ENDPOINT = f"{NDEX_BASE_URL_V2}/networkset"
NCIPID_SET = "8a2d7ee9-1513-11e9-bb6a-0ac135e8bacf"
NDEX_FTP_BASE = (
    "ftp://ftp.ndexbio.org/NCI_PID_BIOPAX_2016-06-08-PC2v8-API/{pathway_name}.owl.gz"
)


def _get_default_ndex_args() -> Dict[str, str]:
    usr, pwd = ndex_client.get_default_ndex_cred(ndex_cred=None)
    return {"username": usr, "password": pwd, "server": NDEX_BASE_URL}


NDEX_ARGS = _get_default_ndex_args()


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


def _get_networks_in_set(network_set_id: str) -> List[str]:
    """Get the UUIDs of the networks in a network set."""
    nd = ndex2.client.Ndex2(
        **{(k if k != "server" else "host"): v for k, v in NDEX_ARGS.items()}
    )
    network_set = nd.get_networkset(network_set_id)
    return network_set["networks"]


def _get_ndex_graph_info(network_id: str) -> Tuple[str, str]:
    """Get name and ftp url for a network

    Parameters
    ----------
    network_id :
        The UUID of the network

    Returns
    -------
    :
        Tuple of name and ftp url
    """
    nd = ndex2.client.Ndex2(
        **{(k if k != "server" else "host"): v for k, v in NDEX_ARGS.items()}
    )
    info_dict = nd.get_network_summary(network_id)
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


def _download_owl_file(ftp_url: str, out_file: str):
    """Download the owl file for a network from NDEx. and return the name"""
    # Download the file at the url for the ftp server to the given file path
    with closing(request.urlopen(ftp_url)) as r:
        with open(out_file, 'wb') as f:
            shutil.copyfileobj(r, f)


def _get_grounding(name):
    """Get the grounding from the gilda"""
    gilda_res = get_grounding(name)
    if gilda_res and gilda_res[1]:  # gilda_res == ({...}, [...])
        # Get first match
        term_match = gilda_res[1][0]["term"]
        return term_match["entry_name"], term_match["db"], term_match["id"]


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

        # If we have a mapping, return it
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
        grounding = _get_grounding(name)
        if grounding:
            return grounding
        # Try to get a grounding without 'family' in the name. This works for
        # e.g. "Gi" vs "Gi family"
        if "family" in name.lower():
            name = name.replace("family", "").strip()
            grounding = _get_grounding(name)
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


def main(sif_file, ncipid_file, merge_how="outer"):
    """Convert the nci-pid CX network to a SIF file and merge with the INDRA SIF

    Parameters
    ----------
    sif_file :
        The INDRA SIF dump file location
    ncipid_file :
        The nci-pid CX dump file location
    merge_how :
        How to merge the INDRA SIF with the nci-pid CX SIF. Allowed values are
        "outer", "left", "right", and "inner". This is passed to the "how"
        parameter for pandas.DataFrame.merge(). The sif dump is "left" and
        the nci-pid CX SIF is "right": sif_df.merge(cx_df, how=merge_how).
    """
    if merge_how not in ["left", "right", "outer", "inner", "cross"]:
        raise ValueError(
            f"Invalid merge_how value {merge_how}. Allowed "
            f"values are 'left', 'right', 'outer', 'inner', 'cross'"
        )
    # Load the CX network
    nci_cx = create_nice_cx_from_file(ncipid_file)

    # Load the INDRA SIF dump
    with open(sif_file, "rb") as fh:
        sif_df: pd.DataFrame = pickle.load(fh)

    # Make a name to NS-ID mapping from the sif dump
    sif_ns_id_map = {
        n: (ns, _id)
        for ns, _id, n in set(zip(sif_df.agA_ns, sif_df.agA_id, sif_df.agA_name)).union(
            set(zip(sif_df.agB_ns, sif_df.agB_id, sif_df.agB_name))
        )
    }

    # Get node id to entity map
    node_id_to_entity = get_node_mapping(nci_cx, sif_ns_id_map)

    # Build the CX SIF file
    logger.info("Building CX SIF file")
    cx_sif = build_cx_sif(nci_cx, node_id_to_entity)

    # Merge the two data frames
    return merge_dfs(sif_df, cx_sif, merge_how=merge_how)


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


def get_missing_edges(merged_df: pd.DataFrame):
    """Get the missing edges from the merged_df"""
    # Get the interactions that are missing in the INDRA SIF
    nci_only_hgnc = merged_df[["agA_ns", "agA_id", "agB_ns", "agB_id", "interaction"]][
        (merged_df.agA_ns == "HGNC")
        & (merged_df.agB_ns == "HGNC")
        & (merged_df._merge == "right_only")
    ]

    # Get the count for each interaction type
    return sum(nci_only_hgnc.interaction.apply(Counter).values, Counter())


def identify_cx_graph_w_missing_edges(dir_path: str, edges):
    from pathlib import Path

    # Edges are hgnc symbols
    # Loop all cx files in the directory
    graph_files = []
    for cx_file in tqdm(
        Path(dir_path).glob("*.cx"), total=len(list(Path(dir_path).glob("*.cx")))
    ):
        if cx_file == "NCI_PID_Complete_Interactions.cx":
            print(f"Skipping {cx_file}")
            continue

        # Load the cx file
        cx_graph = create_nice_cx_from_file(cx_file)

        # Get node to name mapping
        node_to_name = {node: cx_graph.get_node(node)["n"] for node in cx_graph.nodes}

        # Check if the edges are in the graph
        for e_id in cx_graph.edges:
            # Get source and target of edge
            ed = cx_graph.get_edge(e_id)
            source, target = ed["s"], ed["t"]
            source_name = node_to_name[source]
            target_name = node_to_name[target]

            # Check if the edge is in the list of edges (check reverse edge
            # as well)
            if (source_name, target_name) in edges:
                graph_files.append((cx_file, (source_name, target_name)))
            elif (target_name, source_name) in edges:
                graph_files.append((cx_file, (target_name, source_name)))

    print(f"Found {len(graph_files)} missing edges")
    return graph_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge the INDRA SIF with the nci-pid CX SIF"
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
        "--merge-how",
        help="How to merge the INDRA SIF with the nci-pid CX SIF. This is "
        "passed to the 'how' parameter for pandas.DataFrame.merge(). The "
        "sif dump is 'left' and the nci-pid CX SIF is 'right': "
        "sif_df.merge(cx_df, how=merge_how)",
        default="outer",
    )
    parser.add_argument(
        "--out-dir",
        help="The output directory for the merged dataframe and plots",
        default=".",
    )
    args = parser.parse_args()

    # Run the main function if the dataframe does not already exist
    if not os.path.exists(os.path.join(args.out_dir, "merged_df.pkl")):
        df = main(
            args.sif_file,
            args.ncipid_file,
            merge_how=args.merge_how,
        )
        df.to_pickle(os.path.join(args.out_dir, "merged_df.pkl"))
    else:
        df = pd.read_pickle(os.path.join(args.out_dir, "merged_df.pkl"))

    # Plot the venn diagrams
    venn_plots(merged_df=df, out_dir=args.out_dir)
