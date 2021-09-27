import os
import copy
import tqdm
import pickle
import obonet
import logging
import networkx
import itertools
import pandas as pd
from collections import defaultdict
import ndex2.client
from ndex2 import create_nice_cx_from_server
from indra.util import batch_iter
from indra.statements import Complex
from indra.sources import indra_db_rest
import indra.tools.assemble_corpus as ac
from indra.databases import uniprot_client, ndex_client
from indra.assemblers.cx import NiceCxAssembler
from indra.preassembler import Preassembler
from indra.ontology.bio import bio_ontology
from indra.preassembler.custom_preassembly import agents_stmt_type_matches
from indra_db.client.principal.curation import get_curations

logger = logging.getLogger('go_networks')
logging.getLogger('indra.sources.indra_db_rest.util').setLevel(logging.WARNING)

HERE = os.path.dirname(os.path.abspath(__file__))
GO_OBO_PATH = os.path.join(HERE, os.pardir, 'go.obo')
GO_ANNOTS_PATH = os.path.join(HERE, os.pardir, 'goa_human.gaf')
INDRA_SIF_PICKLE = '/Users/ben/data/indradb/2021-08-09/sif.pkl'

go_dag = obonet.read_obo(GO_OBO_PATH)

STMTS_BY_HASH_CACHE = {}


def make_genes_by_go_id(path):
    """Load the gene/GO annotations as a pandas data frame."""
    goa = pd.read_csv(path, sep='\t',
                      skiprows=31, dtype=str,
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

    genes_by_go_id = defaultdict(set)
    for idx, (go_id, up_id) in enumerate(zip(goa.GO_ID, goa.DB_ID)):
        gene_name = uniprot_client.get_gene_name(up_id)
        if gene_name:
            genes_by_go_id[go_id].add(gene_name)
    return genes_by_go_id


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    df = df[(df.agA_ns == 'HGNC') & (df.agB_ns == 'HGNC')]
    logger.info('Filtered to %d rows with genes only' % len(df))
    return df


def get_hashes_by_gene_pair(df):
    hashes_by_gene_pair = defaultdict(set)
    for a, b, hs in zip(df.agA_name, df.agB_name, df.stmt_hash):
        hashes_by_gene_pair[(a, b)].add(hs)
    return hashes_by_gene_pair


def filter_out_medscan(stmts):
    logger.info('Filtering out medscan evidence on %d statements' % len(stmts))
    new_stmts = []
    for stmt in stmts:
        new_evidence = [e for e in stmt.evidence if e.source_api != 'medscan']
        if not new_evidence:
            continue
        stmt.evidence = new_evidence
        if not stmt.evidence:
            continue
        new_stmts.append(stmt)
    logger.info('%d statements after filter' % len(new_stmts))
    return new_stmts


def download_statements(hashes):
    """Download the INDRA Statements corresponding to a set of hashes.
    """
    stmts_by_hash = {h: STMTS_BY_HASH_CACHE[h] for h in hashes
                     if h in STMTS_BY_HASH_CACHE}
    hashes_not_in_cache = {h for h in hashes if h not in STMTS_BY_HASH_CACHE}
    if not hashes_not_in_cache:
        return stmts_by_hash
    for group in tqdm.tqdm(batch_iter(hashes_not_in_cache, 1000),
                           total=int(len(hashes_not_in_cache)/1000)):
        idbp = indra_db_rest.get_statements_by_hash(list(group), ev_limit=10)
        for stmt in idbp.statements:
            stmts_by_hash[stmt.get_hash()] = stmt
            STMTS_BY_HASH_CACHE[stmt.get_hash()] = stmt
    return stmts_by_hash


def expand_complex(stmt):
    """Replace a Complex statement with binary ones."""
    stmts = []
    added = set()
    for m1, m2 in itertools.combinations(stmt.members, 2):
        keys = (m1.entity_matches_key(), m2.entity_matches_key())
        if keys in added:
            continue
        if len(set(keys)) == 1:
            continue
        ordered = sorted([m1, m2], key=lambda x: x.entity_matches_key())
        c = Complex(ordered, evidence=copy.deepcopy(stmt.evidence))
        stmts.append(c)
        added.add(keys)
    return stmts


def get_genes_for_go_id(go_id):
    """Return genes that are annotated with a given go ID."""
    gene_names = genes_by_go_id[go_id]
    for child_go_id in networkx.ancestors(go_dag, go_id):
        gene_names |= genes_by_go_id[child_go_id]
    gene_names = sorted(gene_names)
    return gene_names


def get_go_ids():
    """Get a list of all GO IDs."""
    go_ids = [n for n in go_dag.nodes
              if go_dag.nodes[n]['namespace'] == 'biological_process']
    return go_ids


def format_and_upload_network(ncx, network_set_id, style_ncx,
                              **ndex_args):
    """Take a NiceCXNetwork and upload it to NDEx."""
    ncx.apply_style_from_network(style_ncx)
    network_url = ncx.upload_to(**ndex_args)
    network_id = network_url.split('/')[-1]
    nd = ndex2.client.Ndex2(**{(k if k != 'server' else 'host'): v
                               for k, v in ndex_args.items()})
    nd.make_network_public(network_id)
    nd.add_networks_to_networkset(network_set_id, [network_id])
    return network_id


def get_statement_hashes(go_id):
    """For a given GO ID, find the statement hashes."""
    genes = get_genes_for_go_id(go_id)
    logger.debug('%d genes for %s' % (len(genes), go_id))

    metadata = {}
    hashes = set()
    metadata['num_genes'] = len(genes)
    if len(genes) < min_gene_count:
        metadata['error'] = 'TOO_FEW_GENES'
    elif len(genes) > max_gene_count:
        metadata['error'] = 'TOO_MANY_GENES'
    else:
        for ga, gb in itertools.product(genes, genes):
            hashes |= hashes_by_gene_pair.get((ga, gb), set())
        if not hashes:
            metadata['error'] = 'NO_STATEMENTS'
    return hashes, metadata


def assemble_network_stmts(stmts, genes):
    metadata = {}
    metadata['num_all_raw_stmts'] = len(stmts)
    stmts = filter_out_medscan(stmts)
    metadata['num_filtered_raw_stmts'] = len(stmts)
    all_stmts = []
    for stmt in stmts:
        if isinstance(stmt, Complex):
            all_stmts += expand_complex(stmt)
        else:
            all_stmts.append(stmt)
    # This is to make sure that expanded complexes don't add nodes that
    # shouldn't be in the scope of the network
    all_stmts = ac.filter_gene_list(all_stmts, genes, policy='all')
    pa = Preassembler(bio_ontology, stmts=all_stmts,
                      matches_fun=agents_stmt_type_matches)
    stmts = pa.combine_duplicates()
    stmts = ac.filter_by_curation(stmts, curations=db_curations)
    metadata['num_assembled_stmts'] = len(stmts)
    return stmts, metadata


def make_cx_network(stmts, go_id):
    network_name = '%s (%s)' % (go_id, go_dag.nodes[go_id]['name'])
    logger.info('===============================')
    network_attributes = {
        'networkType': 'pathway',
        'GO ID': go_id,
        'GO hierarchy': 'biological process',
        'Prov:wasGeneratedBy': 'INDRA',
        'Organism': 'Homo sapiens (Human)',
        'Description': go_dag.nodes[go_id]['name'],
        'Methods': 'This network was assembled automatically by INDRA ('
                   'http://indra.bio) by processing all available '
                   'biomedical literature with multiple machine reading '
                   'systems, and integrating curated pathway '
                   'databases. The network represents '
                   'mechanistic interactions between genes/proteins that '
                   'are associated with this GO process.',
    }
    # Assemble CX network
    ca = NiceCxAssembler(stmts, network_name)
    ncx = ca.make_model(self_loops=False,
                        network_attributes=network_attributes)
    return ncx


if __name__ == '__main__':
    min_gene_count = 5
    max_gene_count = 200
    network_set_id = 'd4b51854-0b2b-11ec-b666-0ac135e8bacf'
    style_network_id = '058c452f-b0d6-11ea-a4d3-0660b7976219'
    style_ncx = create_nice_cx_from_server(
        server='http://test.ndexbio.org',
        uuid=style_network_id)
    username, password = ndex_client.get_default_ndex_cred(ndex_cred=None)
    ndex_args = {'server': 'http://public.ndexbio.org',
                 'username': username,
                 'password': password}

    indra_df = load_indra_df(INDRA_SIF_PICKLE)
    hashes_by_gene_pair = get_hashes_by_gene_pair(indra_df)
    genes_by_go_id = make_genes_by_go_id(GO_ANNOTS_PATH)
    go_ids = get_go_ids()

    # Stage 1. Get all hashes for each GO ID
    metadata = {}
    network_hashes = {}
    for go_id in tqdm.tqdm(go_ids):
        network_hashes[go_id], metadata[go_id] = get_statement_hashes(go_id)
        metadata[go_id]['go_name'] = go_dag.nodes[go_id]['name']

    with open('network_hashes.pkl', 'wb') as fh:
        pickle.dump(network_hashes, fh)

    with open('metadata.pkl', 'wb') as fh:
        pickle.dump(metadata, fh)

    # Stage 2. Download all statements by hash
    if not os.path.exists('stmts_by_hash.pkl'):
        all_hashes = set.union(*network_hashes.values())
        stmts_by_hash = download_statements(all_hashes)
        with open('stmts_by_hash.pkl', 'wb') as fh:
            pickle.dump(stmts_by_hash, fh)
    else:
        with open('stmts_by_hash.pkl', 'rb') as fh:
            stmts_by_hash = pickle.load(fh)

    # Stage 3. Assemble statements for each GO ID
    db_curations = get_curations()
    networks = {}
    for go_id in go_ids:
        logger.info('Assembling network for %s' % go_id)
        if metadata[go_id].get('error'):
            continue
        stmts = [stmts_by_hash.get(h) for h in network_hashes[go_id]]
        stmts = {stmt for stmt in stmts if stmt}
        logger.info('Got %d statements by hash' % len(stmts))
        metadata[go_id]['new_raw_stmts_obtained'] = len(stmts)
        if not stmts:
            metadata[go_id]['error'] = 'NO_RAW_STMTS'
            continue
        stmts, md = assemble_network_stmts(stmts, genes_by_go_id[go_id])
        metadata[go_id].update(md)
        if not stmts:
            metadata[go_id]['error'] = 'NO_ASSEMBLED_STMTS'
            continue
        ncx = make_cx_network(stmts, go_id)
        metadata[go_id]['num_ncx_nodes'] = len(ncx.nodes)
        if not ncx.nodes:
            metadata[go_id]['error'] = 'NO_NETWORK_NODES'
            logger.info('Skipping: no nodes in network.')
            continue
        metadata[go_id]['included'] = True
        networks[go_id] = ncx

    with open('networks.pkl', 'wb') as fh:
        pickle.dump(networks, fh)

    with open('metadata.pkl', 'wb') as fh:
        pickle.dump(metadata, fh)

    # Stage 4. Upload networks
    for go_ids, ncx in networks.items():
        network_id = format_and_upload_network(ncx, network_set_id,
                                               style_ncx, **ndex_args)
        logger.info('Uploaded network with ID: %s' % network_id)
