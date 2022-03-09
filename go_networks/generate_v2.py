"""
Generate GO Networks from a list of GO terms and the Sif dump.
"""
import argparse
from collections import defaultdict
import logging
import pickle
import obonet
from itertools import combinations
from pathlib import Path
from typing import Optional, Dict, Set, Tuple, List

import ndex2.client
from ndex2 import create_nice_cx_from_server
import networkx as nx
import pandas as pd
from tqdm import tqdm

from go_networks.util import (
    DIRECTED_TYPES,
    UNDIRECTED_TYPES,
    load_latest_sif,
)
from go_networks.network_assembly import GoNetworkAssembler
from indra_db.client.principal.curation import get_curations
from indra.databases import uniprot_client, hgnc_client
from indra.util.statement_presentation import reader_sources

# Derived types
Go2Genes = Dict[str, Set[str]]
EvCountDict = Dict[str, Dict[str, int]]
NameEntityMap = Dict[str, Tuple[str, str]]

# Constants
HERE = Path(__file__).absolute().parent.parent
GO_ANNOTS_PATH = HERE.joinpath("goa_human.gaf").absolute().as_posix()
GO_OBO_PATH = HERE.joinpath("go.obo").absolute().as_posix()
PROPS_FILE = HERE.joinpath("props.pkl").absolute().as_posix()
NETWORKS_FILE = HERE.joinpath("networks.pkl").absolute().as_posix()

min_gene_count = 5
max_gene_count = 200

logger = logging.getLogger(__name__)


def get_curation_set() -> Set[int]:
    correct_tags = {"correct", "hypothesis", "activity_amount"}
    try:
        curations = get_curations()
        # curations == {'pa_hash': 123456, 'tag': '<grounding tag>'}
        # Only keep the hashes for the curations that are not correct
        wrong_hashes = {
            c['pa_hash'] for c in curations if c['tag'] not in correct_tags
        }
        logger.info(f"Found {len(wrong_hashes)} hashes to filter out")
        return wrong_hashes
    except Exception as e:
        logger.error(f"Could not get curations: {e}")


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
        ~((sif_df.evidence_count == 1) &
          (sif_df.source_counts.apply(
            lambda d: d is None or
            (len(d) == 1 and bool((set(d) & reader_sources_set))))
          ))
    ]
    t.update()

    # Remove all rows where the source is sparser and the stmt type is Complex
    sif_df = sif_df[~(
            (sif_df.stmt_type == "Complex") &
            (sif_df.source_counts.apply(lambda d: d is None or
                                        (set(d) == {"sparser"})))
    )]
    t.update()
    t.close()

    return sif_df


def generate_props(
    sif_file: str, props_file: Optional[str] = None, apply_filters: bool = True
) -> Dict[str, List[Dict[str, int]]]:
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
    skipped = 0
    # Only pass the relevant parts of the pair_props dict
    for go_id, gene_set in tqdm(go2genes_map.items(), total=len(go2genes_map)):
        def _key(g1, g2):
            return tuple(sorted([g1, g2]))
        # Get relevant pairs from pair_properties
        prop_dict = {_key(g1, g2): pair_props[_key(g1, g2)]
                     for g1, g2 in combinations(gene_set, 2)
                     if _key(g1, g2) in pair_props}

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
        networks[go_id] = gna.network

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
    return {go_id: genes for go_id, genes in go2genes_map.items()
            if min_gene_count <= len(genes) <= max_gene_count}


def generate(sif_file: Optional[str] = None,
             props_file: Optional[str] = None,
             apply_filters: bool = True):
    """Generate new GO networks from INDRA statements

    Parameters
    ----------
    sif_file :
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
    sif_props = generate_props(sif_file, props_file, apply_filters=apply_filters)

    # Iterate by GO ID and for each list of genes, build a network
    return build_networks(go2genes_map, sif_props)


def format_and_upload_network(ncx, network_set_id, style_ncx,
                              **ndex_args):
    """Take a NiceCXNetwork and upload it to NDEx."""
    ncx.apply_style_from_network(style_ncx)
    network_url = ncx.upload_to(**ndex_args)
    network_id = network_url.split('/')[-1]
    nd = ndex2.client.Ndex2(**{(k if k != 'server' else 'host'): v
                               for k, v in ndex_args.items()})
    try:
        nd.make_network_public(network_id)
    except Exception as e:
        logger.warning(f"Failed to make network {network_id} public: {e}")

    try:
        nd.add_networks_to_networkset(network_set_id, [network_id])
    except Exception as e:
        logger.warning(f"Failed to add network {network_id} to network set "
                       f"{network_set_id}: {e}")

    return network_id


if __name__ == "__main__":
    # Setup argument parser
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--local-sif',
                        help='Local SIF dump file to load. If not provided, '
                             'the latest SIF dump will be fetched from S3.')
    parser.add_argument('--network-set',
                        help='Network set ID to add the new networks to.',
                        default='bdba6a7a-488a-11ec-b3be-0ac135e8bacf')
    parser.add_argument('--style-network',
                        help='Network ID of the style network',
                        default='4c2006cd-9fef-11ec-b3be-0ac135e8bacf')
    parser.add_argument('--regenerate', action='store_true',
                        help='Regenerate the networks and re-cache them')
    args = parser.parse_args()
    logger.info(f"Using network set id {args.network_set}")
    logger.info(f"Using network style {args.style_network}")

    if Path(NETWORKS_FILE).exists() and not args.regenerate:
        with open(NETWORKS_FILE, 'rb') as fh:
            networks = pickle.load(fh)
    else:
        if args.regenerate:
            logger.info("Regenerating props and networks")
            props_file = None
        else:
            props_file = PROPS_FILE
        networks = generate(args.local_sif, props_file, apply_filters=True)
        with open(NETWORKS_FILE, 'wb') as f:
            pickle.dump(networks, f)

    style_ncx = create_nice_cx_from_server(
        server='http://ndexbio.org',
        uuid=args.style_network)

    from indra.databases import ndex_client
    username, password = ndex_client.get_default_ndex_cred(ndex_cred=None)
    ndex_args = {'server': 'http://public.ndexbio.org',
                 'username': username,
                 'password': password}
    for go_id, network in tqdm(sorted(networks.items(), key=lambda x: x[0])):
        network_id = format_and_upload_network(network, args.network_set,
                                               style_ncx, **ndex_args)
