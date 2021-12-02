"""
Load the nci-pid CX network and convert it to a SIF file so be compared with
the INDRA SIF dump.
"""
import argparse
import logging
import os
import pickle
from typing import Tuple, Optional, List, Dict, Union

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from tqdm import tqdm

from indra.statements.agent import default_ns_order
from indra.ontology.bio import bio_ontology
from indra.preassembler.grounding_mapper.gilda import get_grounding
from indra.util.statement_presentation import reader_sources, db_sources

from protmapper.uniprot_client import get_primary_id

from ndex2 import create_nice_cx_from_file

logger = logging.getLogger(__name__)


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


def _get_grounding(name):
    """Get the grounding from the gilda"""
    gilda_res = get_grounding(name)
    if gilda_res and gilda_res[1]:  # gilda_res == ({...}, [...])
        # Get first match
        term_match = gilda_res[1][0]["term"]
        return term_match["entry_name"], term_match["db"], term_match["id"]


def _get_node_info(
    cx, cx_id, sif_name_map, use_gilda
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
        hgnc_id = _nim(name)

        # If we have a mapping, return it
        if hgnc_id is not None:
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
    cx, sif_name_map, use_gilda=False
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

    # Drop duplicates and reset the index
    cx_sif.drop_duplicates(inplace=True)
    cx_sif.reset_index(inplace=True, drop=True)
    return cx_sif


def merge_dfs(sif, cx, merge_how="outer") -> pd.DataFrame:
    # Add a new columns to both of the data frames that maps statement
    # type/interaction to boolean indicating if the interaction is directed
    # direction of the interaction
    logger.info("Adding directed column to the DataFrames")
    sif["directed"] = sif.stmt_type != "Complex"
    cx["directed"] = cx.interaction != "in-complex-with"

    # Group the CX SIF by entity pair A B and directed
    logger.info("Grouping the CX SIF by entity pair A B and directed")
    cx = (
        cx.groupby(["agA_ns", "agA_id", "agB_ns", "agB_id", "directed"])
        # Aggregate pmids to list of lists
        .aggregate(
            {"interaction": pd.Series.tolist, "pmids": pd.Series.tolist}
        ).reset_index(["agA_ns", "agA_id", "agB_ns", "agB_id", "directed"])
    )

    # Group the INDRA SIF by entity and stmt type
    logger.info("Grouping the INDRA SIF by entity and stmt type")
    sif = (
        sif.groupby(["agA_ns", "agA_id", "agB_ns", "agB_id", "directed"])
        .aggregate(
            {
                "evidence_count": np.sum,
                "stmt_hash": pd.Series.tolist,
                "belief": pd.Series.tolist,
                "source_counts": pd.Series.tolist,
                "stmt_type": pd.Series.tolist,
            }
        )
        .reset_index(["agA_ns", "agA_id", "agB_ns", "agB_id", "directed"])
    )

    # Merge the two DataFrames on A-B pairs + interaction type
    logger.info("Merging the CX SIF with the INDRA SIF")
    merged_df = sif.merge(
        cx,
        on=["agA_ns", "agA_id", "agB_ns", "agB_id", "directed"],
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
        How to merge the INDRA SIF with the nci-pid CX SIF. This is passed to
        the "how" parameter for pandas.DataFrame.merge(). The sif dump is
        "left" and the nci-pid CX SIF is "right":
        sif_df.merge(cx_df, how=merge_how)
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
    cx_tot_ix = set(merged_df.query("_merge in ['right_only', 'both']").index.to_list())
    t.update()
    sif_merged: pd.DataFrame = merged_df.query("_merge in ['left_only', 'both']")
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
            sif_merged.source_counts.apply(lambda scl: _reader_count_filter(
                scl, 3))
        ].index.to_list()
    )
    t.update()

    # Has DB support
    sif_has_db_ix = set(
        sif_merged[
            sif_merged.source_counts.apply(lambda scl: _db_count_filter(scl, 1))
        ].index.to_list()
    )
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
    fig.savefig(os.path.join(out_dir, "venn_plots.png"))


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
