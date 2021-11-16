"""
Batch read text_refs.tsv and statements.tsv to align PMC IDs with evidence
counts
"""
from pathlib import Path
import json
import logging
from indra.statements import Statement
import codecs

logger = logging.getLogger(__name__)


def generate_reading_pmc_mapping(text_refs_tsv):
    tsv_path = Path(text_refs_tsv) if isinstance(text_refs_tsv, str) else text_refs_tsv

    id_names = ['TRID', 'PMID', 'PMCID', 'DOI', 'PII', 'URL', 'MANUSCRIPT_ID']
    pmc_ix = id_names.index('PMCID') + 1
    pmc_map = {}

    # Open tsv file and save mapping to pmc_map_file
    with tsv_path.open('r') as fi:
        line = fi.readline()
        while line:
            ll = line.split('\t')

            # columns are [<reading id>, *id_names]
            reading_id = ll[0]
            pmc_id = ll[pmc_ix]

            if pmc_id != r'\N':
                pmc_map[reading_id] = pmc_id

            line = fi.readline()

    return pmc_map


def load_statement_tsv_to_dict(stmts_tsv, ignore_ungrounded=True,
                               ignore_sources=None):
    stmts_path = Path(stmts_tsv) if isinstance(stmts_tsv, str) else stmts_tsv
    reading_counts = {}

    # Open tsv file and save statements to evidence_counts_tsv
    with stmts_path.open('r') as fi:
        line = fi.readline()
        while line:
            ll = line.strip().split('\t')

            # columns are [<raw id>, <reading id>, <hash>, <stmt_json>]
            reading_id = ll[1]
            if reading_id == r'\N':
                line = fi.readline()
                continue

            try:
                raw_json = ll[3]
                esc_json = codecs.escape_decode(
                    codecs.escape_decode(raw_json)[0].decode()
                )[0].decode()
                stmt = Statement._from_json(json.loads(esc_json))
            except json.JSONDecodeError as e:
                print(f'Could not parse JSON for {reading_id}: {e}')
                line = fi.readline()
                continue

            # Check grounded if requested
            if ignore_ungrounded and \
                    any((a is None or a.get_grounding() == (None, None))
                        for a in stmt.agent_list()):
                line = fi.readline()
                continue

            # Ignore sources if given
            if ignore_sources and stmt.evidence[0].source_api in ignore_sources:
                line = fi.readline()
                continue

            try:
                reading_counts[reading_id] += 1
            except KeyError:
                reading_counts[reading_id] = len(stmt.evidence)

            line = fi.readline()

    return reading_counts
