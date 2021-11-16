"""
Batch read text_refs.tsv and statements.tsv to align PMC IDs with evidence
counts
"""
from pathlib import Path
import json
import logging

logger = logging.getLogger(__name__)


def generate_reading_pmc_mapping(text_refs_tsv, pmc_map_file) -> Path:
    tsv_path = Path(text_refs_tsv) if isinstance(text_refs_tsv, str) else text_refs_tsv
    map_path = Path(pmc_map_file) if isinstance(pmc_map_file, str) else pmc_map_file

    id_names = ['TRID', 'PMID', 'PMCID', 'DOI', 'PII', 'URL', 'MANUSCRIPT_ID']
    pmc_ix = id_names.index('PMCID') + 1

    # Open tsv file and save mapping to pmc_map_file
    with tsv_path.open('r') as fi, map_path.open('w') as fo:
        line = fi.readline()
        while line:
            ll = fi.readline().split('\t')

            # columns are [<reading id>, *id_names]
            reading_id = ll[0]
            pmc_id = ll[pmc_ix]

            if pmc_id != r'\N':
                fo.write(f'{reading_id}\t{pmc_id}\n')

            line = fi.readline()

    return map_path
