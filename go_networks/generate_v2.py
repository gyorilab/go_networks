"""
Generate GO Networks from a list of GO terms and the Sif dump.
"""
from collections import defaultdict
import logging
import pickle
import obonet
from itertools import combinations
from pathlib import Path
from typing import Optional, Dict, Set, Tuple, Union

import ndex2.client
from ndex2 import create_nice_cx_from_server, NiceCXNetwork
import networkx as nx
import pandas as pd
from tqdm import tqdm

from go_networks.util import (
    DIRECTED_TYPES,
    load_latest_sif,
)
from go_networks.network_assembly import GoNetworkAssembler
from indra_db.client.principal.curation import get_curations
from indra.databases import uniprot_client, hgnc_client
from indra.util.statement_presentation import reader_sources

# Derived types
PropAgg = Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]
PropDict = Dict[Tuple[str, str], PropAgg]
Go2Genes = Dict[str, Set[str]]
EvCountDict = Dict[str, Dict[str, int]]
NameEntityMap = Dict[str, Tuple[str, str]]

# Constants
HERE = Path(__file__).absolute().parent.parent
GO_ANNOTS_PATH = HERE.joinpath("goa_human.gaf").absolute().as_posix()
GO_OBO_PATH = HERE.joinpath("go.obo").absolute().as_posix()
PROPS_FILE = HERE.joinpath("props.pkl").absolute().as_posix()
NETWORKS_FILE = HERE.joinpath("networks.pkl").absolute().as_posix()
DEFAULT_NDEX_SERVER = "http://ndexbio.org"
TEST_GO_ID = None

min_gene_count = 5
max_gene_count = 200

logger = logging.getLogger(__name__)


def get_curation_set() -> Set[int]:
    def _is_incorrect(stmt_hash, evid_hash):
        # Evidence is incorrect if it was only curated as incorrect
        if (
            evid_hash in correct_stmt_evid[stmt_hash]["incorrect"]
            and evid_hash not in correct_stmt_evid[stmt_hash]["correct"]
        ):
            return True
        return False

    correct_tags = {"correct", "hypothesis", "act_vs_amt"}

    try:
        curations = get_curations()
    except Exception as e:
        logger.error(f"Could not get curations: {e}")
        return set()

    # If a statement only has incorrect curations, it is overall incorrect
    # If a statement has any correct curations, it is overall correct
    correct = {c["pa_hash"] for c in curations if c["tag"] in correct_tags}
    incorrect = {c["pa_hash"] for c in curations if c["pa_hash"] not in correct}

    correct_stmt_evid = {}
    for c in curations:
        pa_hash = c["pa_hash"]
        if pa_hash in correct:
            if pa_hash not in correct_stmt_evid:
                correct_stmt_evid[pa_hash] = defaultdict(set)
            if c["tag"] in correct_tags:
                correct_stmt_evid[pa_hash]["correct"].add(c["source_hash"])
            else:
                correct_stmt_evid[pa_hash]["incorrect"].add(c["source_hash"])

    hashes_out = set()
    for c in curations:
        pa_hash = c["pa_hash"]
        src_hash = c["source_hash"]
        if pa_hash not in incorrect:
            # Check evidences in detail
            if pa_hash in correct_stmt_evid:
                if _is_incorrect(pa_hash, src_hash):
                    hashes_out.add(pa_hash)
                else:
                    pass
            else:
                pass
        else:
            hashes_out.add(pa_hash)

    logger.info(f"Found {len(hashes_out)} hashes to filter out")

    return hashes_out


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


def quality_filter(sif_df: pd.DataFrame) -> pd.DataFrame:
    """Apply quality filters to the SIF dataframe

    Filters applied:
        - Remove all rows that are curated as incorrect
        - Remove all rows with only one evidence that is from a reader
        - Remove all rows where there is only one source, the source is
          sparser and the stmt type is Complex

    Parameters
    ----------
    sif_df : pd.DataFrame
        The SIF dataframe

    Returns
    -------
    pd.DataFrame
        The filtered SIF dataframe
    """
    # Filter out curated incorrect statements
    wrong_hashes = get_curation_set()
    if wrong_hashes:
        t = tqdm(desc="Quality filtering", total=3)
        sif_df = sif_df[~sif_df.stmt_hash.isin(wrong_hashes)]
        t.update()
    else:
        t = tqdm(desc="Quality filtering", total=2)

    # Filter out statements with only one evidence from a reader source
    reader_sources_set = set(reader_sources)
    sif_df = sif_df[
        ~(
            (sif_df.evidence_count == 1)
            & (
                sif_df.source_counts.apply(
                    lambda d: d is None
                    or (len(d) == 1 and bool((set(d) & reader_sources_set)))
                )
            )
        )
    ]
    t.update()

    # Remove all rows where the source is sparser and the stmt type is Complex
    sif_df = sif_df[
        ~(
            (sif_df.stmt_type == "Complex")
            & (
                sif_df.source_counts.apply(
                    lambda d: d is None or (set(d) == {"sparser"})
                )
            )
        )
    ]
    t.update()
    t.close()

    return sif_df


def generate_props(
    sif_file: Optional[str] = None,
    props_file: Optional[str] = None,
    apply_filters: bool = True,
    regenerate: bool = False,
) -> PropDict:
    """Generate properties per pair from the Sif dump

    For each pair of genes (A,B) (excluding self loops), generate the
    following properties:
        - "SOURCE => TARGET": aggregate number of evidences by statement
          type for A->B statements
        - "TARGET => SOURCE": aggregate number of evidences by statement
          type for B->A statements
        - "SOURCE - TARGET": aggregate number of evidences by statement type
          for A-B undirected statements

    Parameters
    ----------
    sif_file :
        The SIF dump file
    props_file :
        The file to cache the properties to
    apply_filters :
        Whether to apply quality filters to the SIF dataframe
    regenerate :
        Whether to regenerate the props. If a props file path is provided,
        the old props file will be overwritten.

    Returns
    -------
    :
        The properties dictionary
    """
    # If the props file exists, load it, unless we want to regenerate
    if not regenerate and props_file is not None and Path(props_file).is_file():
        logger.info(f"Loading property lookup from {props_file}")
        with Path(props_file).open(mode="rb") as fr:
            props_by_pair = pickle.load(fr)
    else:
        logger.info("Generating property lookup")

        # Load the latest INDRA SIF dump
        sif_df = get_sif(sif_file)

        # Filter to HGNC-only rows
        sif_df = filter_to_hgnc(sif_df)

        # Filter out self-loops
        sif_df = filter_self_loops(sif_df)

        # Run filters if enabled
        if apply_filters:
            logger.info("Applying quality filters")
            sif_df = quality_filter(sif_df)

        hashes_by_pair = defaultdict(set)
        props_by_hash = {}

        def get_direction(row, pair):
            directed = row.stmt_type in DIRECTED_TYPES
            direction = row.agA_name == pair[0]
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
                    "direction": get_direction(row, pair),
                }
        hashes_by_pair = dict(hashes_by_pair)

        def aggregate_props(props) -> PropAgg:
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
            props_by_pair[pair] = aggregate_props([props_by_hash[h] for h in hashes])

        # Write to file if provided
        if props_file:
            Path(props_file).absolute().parent.mkdir(exist_ok=True, parents=True)
            with Path(props_file).open(mode="wb") as fo:
                logger.info(f"Saving property lookup to {props_file}")
                pickle.dump(obj=props_by_pair, file=fo)

    return props_by_pair


def genes_by_go_id():
    """Load the gene/GO annotations as a pandas data frame."""
    go_dag = obonet.read_obo(GO_OBO_PATH)

    goa = pd.read_csv(
        GO_ANNOTS_PATH,
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
    # Filter out all "NOT" negative evidences
    goa["Qualifier"].fillna("", inplace=True)
    goa = goa[~goa["Qualifier"].str.startswith("NOT")]
    # We can filter to just GO terms in the ontology since
    # obsolete terms are not included in the GO DAG
    goa = goa[goa["GO_ID"].isin(go_dag)]

    genes_by_go_id = defaultdict(set)
    for go_id, up_id in zip(goa.GO_ID, goa.DB_ID):
        if go_dag.nodes[go_id]["namespace"] != "biological_process":
            continue
        hgnc_id = uniprot_client.get_hgnc_id(up_id)
        if hgnc_id:
            gene_name = hgnc_client.get_hgnc_name(hgnc_id)
            genes_by_go_id[go_id] = genes_by_go_id[go_id] | {gene_name}

    for go_id in set(genes_by_go_id):
        for child_go_id in nx.ancestors(go_dag, go_id):
            genes_by_go_id[go_id] |= genes_by_go_id[child_go_id]

    return genes_by_go_id


def build_networks(
    go2genes_map: Go2Genes,
    pair_props: PropDict,
) -> Dict[str, Dict[str, Union[NiceCXNetwork, float]]]:
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
    skipped = 0
    # Only pass the relevant parts of the pair_props dict
    for go_id, gene_set in tqdm(go2genes_map.items(), total=len(go2genes_map)):
        if TEST_GO_ID and go_id != TEST_GO_ID:
            continue

        def _key(g1, g2):
            return tuple(sorted([g1, g2]))

        # Get relevant pairs from pair_properties
        prop_dict = {
            _key(g1, g2): pair_props[_key(g1, g2)]
            for g1, g2 in combinations(gene_set, 2)
            if _key(g1, g2) in pair_props
        }

        if not prop_dict:
            # logger.info(f"No statements for ID {go_id}")
            skipped += 1
            continue

        gna = GoNetworkAssembler(
            identifier=go_id,
            entity_list=list(gene_set),
            pair_properties=prop_dict,
        )
        gna.assemble()
        networks[go_id] = {
            "network": gna.network,
            "max_score": max(gna.rel_scores),
            "min_score": min(gna.rel_scores),
        }

    logger.info(f"Skipped {skipped} networks without statements")
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
    return {
        go_id: genes
        for go_id, genes in go2genes_map.items()
        if min_gene_count <= len(genes) <= max_gene_count
    }


def generate(
    sif_file: Optional[str] = None,
    props_file: Optional[str] = None,
    apply_filters: bool = True,
    regenerate: bool = False,
) -> Dict[str, Dict[str, Union[NiceCXNetwork, float]]]:
    """Generate new GO networks from INDRA statements

    Parameters
    ----------
    sif_file :
        If provided, load sif dump from this file. Default: load from S3.
    props_file :
        If provided, load property lookup from this file. Default: generate
        from sif dump.
    apply_filters :
        If True, apply filters to the data before generating networks.
    regenerate :
        If True, regenerate the props from scratch and overwrite the existing
        props file, if it exists.
    """
    # Make genes by GO ID dict
    go2genes_map = genes_by_go_id()

    # Filter GO IDs
    go2genes_map = filter_go_ids(go2genes_map)

    # Generate properties
    sif_props = generate_props(
        sif_file=sif_file,
        props_file=props_file,
        apply_filters=apply_filters,
        regenerate=regenerate,
    )

    # Iterate by GO ID and for each list of genes, build a network
    return build_networks(go2genes_map, sif_props)


def format_and_upload_network(
    ncx: NiceCXNetwork,
    network_set_id: str,
    style_ncx: NiceCXNetwork,
    ndex_client: ndex2.client.Ndex2,
) -> Tuple[str, bool]:
    """Take a NiceCXNetwork and upload it to NDEx."""
    # Fixme: NiceCXNetwork.apply_style_from_network() does not exist in
    #  ndex2==2.0.1 but ==2.0.1 is the requirement for INDRA
    #  This method was added in 3.1.0:
    #  https://github.com/ndexbio/ndex2-client/issues/43
    ncx.apply_style_from_network(style_ncx)
    network_url = ncx.upload_to(client=ndex_client)
    network_id = network_url.split("/")[-1]
    failed = False
    try:
        ndex_client.make_network_public(network_id)
    except Exception as e:
        logger.warning(f"Failed to make network {network_id} public: {e}")
        failed = True

    try:
        ndex_client.add_networks_to_networkset(network_set_id, [network_id])
    except Exception as e:
        logger.warning(
            f"Failed to add network {network_id} to network set "
            f"{network_set_id}: {e}"
        )

    return network_id, failed


def _update_style_network(style_ncx: NiceCXNetwork, min_score: float, max_score: float):
    if NiceCXNetwork.CY_VISUAL_PROPERTIES in style_ncx.opaqueAspects:
        vis_prop = NiceCXNetwork.CY_VISUAL_PROPERTIES
        for style_dict in style_ncx.opaqueAspects[vis_prop]:
            if style_dict["properties_of"] == "edges:default":
                style_str = style_dict["mappings"]["EDGE_WIDTH"]["definition"]
                styles = style_str.split(",")
                min_ix = -1
                max_ix = -1
                for st in styles:
                    if "OV=0=" in st:
                        min_ix = styles.index(st)
                    if "OV=1=" in st:
                        max_ix = styles.index(st)
                    if min_ix != -1 and max_ix != -1:
                        break
                styles[min_ix] = f"OV=0={min_score}"
                styles[max_ix] = f"OV=1={max_score}"
                style_dict["mappings"]["EDGE_WIDTH"]["definition"] = ",".join(styles)
                break

        assert f"OV=0={min_score}" in style_dict["mappings"]["EDGE_WIDTH"]["definition"]
        assert f"OV=1={max_score}" in style_dict["mappings"]["EDGE_WIDTH"]["definition"]


def main(
    network_set: str,
    style_network: str,
    regenerate: bool,
    local_sif: Optional[str] = None,
    test_go_term: Optional[str] = None,
    ndex_server_style: str = DEFAULT_NDEX_SERVER,
):
    global TEST_GO_ID
    logger.info(f"Using network set id {network_set}")
    logger.info(f"Using network style {style_network}")
    if test_go_term:
        logger.info(f"Testing GO term {test_go_term}")
        TEST_GO_ID = test_go_term

    logger.info(f"Using ndex server {ndex_server_style} for style")

    networks = generate(
        sif_file=local_sif,
        props_file=PROPS_FILE,
        apply_filters=True,
        regenerate=regenerate,
    )

    # Only cache the networks if we're not testing
    if not test_go_term:
        with open(NETWORKS_FILE, "wb") as f:
            logger.info(f"Writing networks to {NETWORKS_FILE}")
            pickle.dump(networks, f)

    style_ncx = create_nice_cx_from_server(server=ndex_server_style, uuid=style_network)

    from indra.databases import ndex_client

    username, password = ndex_client.get_default_ndex_cred(ndex_cred=None)
    ndex_args = {
        "server": "http://public.ndexbio.org",
        "username": username,
        "password": password,
    }
    ndex_web_client = ndex2.client.Ndex2(
        **{(k if k != "server" else "host"): v for k, v in ndex_args.items()}
    )

    failed_to_set_public = []
    for go_id, network_dict in tqdm(sorted(networks.items(), key=lambda x: x[0])):
        network = network_dict["network"]
        min_score = network_dict["min_score"]
        max_score = network_dict["max_score"]

        # Update style network
        _update_style_network(style_ncx, min_score=min_score, max_score=max_score)

        network_id, failed = format_and_upload_network(
            network, network_set, style_ncx, ndex_client=ndex_web_client
        )
        if failed:
            failed_to_set_public.append(network_id)

        if TEST_GO_ID:
            print(f"Testing network uuid {network_id}")

    if failed_to_set_public:
        logger.warning(f"Failed to set {len(failed_to_set_public)} networks public")
        logger.info("Retrying to set networks public...")

        # Retry to set public
        for failed_uuid in failed_to_set_public:
            try:
                ndex_web_client.make_network_public(network_id=failed_uuid)
            except Exception as e:
                logger.warning(f"Failed to set network {failed_uuid} public again: {e}")
