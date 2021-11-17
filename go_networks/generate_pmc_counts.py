"""
Batch read text_refs.tsv and statements.tsv to align PMC IDs with evidence
counts and output the results to a new file.
"""
from pathlib import Path
import json
import logging
from indra.statements import Statement
from tqdm import tqdm
import codecs

logger = logging.getLogger(__name__)


def generate_reading_pmc_mapping(text_refs_tsv):
    tsv_path = Path(text_refs_tsv) if isinstance(text_refs_tsv, str) else text_refs_tsv

    id_names = ['TRID', 'PMID', 'PMCID', 'DOI', 'PII', 'URL', 'MANUSCRIPT_ID']
    pmc_ix = id_names.index('PMCID') + 1
    pmc_map = {}

    # Open tsv file and save mapping to pmc_map_file
    logger.info('Generating reading id -> pmc map from text_refs.tsv')
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


def statement_tsv_to_reading_counts(stmts_tsv, ignore_ungrounded=True,
                                    ignore_sources=None):
    stmts_path = Path(stmts_tsv) if isinstance(stmts_tsv, str) else stmts_tsv
    reading_counts = {}

    logger.info('Generating reading id -> ev count dict')
    with stmts_path.open('r') as fi:
        line = fi.readline()
        while line:
            ll = line.strip().split('\t')

            # columns are [<raw id>, <reading id>, <hash>, <stmt_json>]
            reading_id = ll[1]
            if reading_id == r'\N':
                line = fi.readline()
                continue

            raw_json = ll[3]
            esc_json = codecs.escape_decode(
                codecs.escape_decode(raw_json)[0].decode()
            )[0].decode()
            stmt = Statement._from_json(json.loads(esc_json))

            # Check grounded if requested
            if ignore_ungrounded and \
                    any((a is None or a.get_grounding() == (None, None))
                        for a in stmt.agent_list()):
                line = fi.readline()
                continue

            # Ignore sources if given
            if ignore_sources:
                count = 0
                for ev in stmt.evidence:
                    if ev.source_api not in ignore_sources:
                        count += 1
                if count == 0:
                    line = fi.readline()
                    continue

            else:
                count = len(stmt.evidence)

            try:
                reading_counts[reading_id] += count
            except KeyError:
                reading_counts[reading_id] = count

            line = fi.readline()

    return reading_counts


def main(text_refs_tsv, stmts_tsv, pmc_map_file, ignore_ungrounded=True,
         ignore_sources=None):
    pmc_map = generate_reading_pmc_mapping(text_refs_tsv)
    reading_counts = statement_tsv_to_reading_counts(stmts_tsv,
                                                     ignore_ungrounded,
                                                     ignore_sources)

    # Aggregate counts per PMC ID: PMC ID -> reading_id is a one to many
    # mapping
    logger.info('Aggregating counts per PMC ID')
    pmc_counts = {}
    for reading_id, count in tqdm(reading_counts.items(), total=len(reading_counts)):
        pmc_id = pmc_map.get(reading_id)
        if pmc_id is not None:
            try:
                pmc_counts[pmc_id] += count
            except KeyError:
                pmc_counts[pmc_id] = count

    # Write to pmc_map_file
    logger.info('Writing pmc_map to %s' % pmc_map_file)
    with open(pmc_map_file, 'w') as fo:
        for pmc_id, count in tqdm(pmc_counts.items(), total=len(pmc_counts)):
            try:
                fo.write(f'{pmc_id}\t{count}\n')
            except KeyError:
                continue


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('text_refs_tsv', help='Path to text_refs.tsv')
    parser.add_argument('stmts_tsv', help='Path to statements.tsv')
    parser.add_argument('pmc_map_file', help='Path to output file')
    parser.add_argument('--ignore-ungrounded', action='store_true',
                        help='Ignore ungrounded statements')
    parser.add_argument('--ignore-sources', nargs='+',
                        help='Ignore statements from these sources')
    args = parser.parse_args()
    main(args.text_refs_tsv, args.stmts_tsv, args.pmc_map_file,
         ignore_ungrounded=args.ignore_ungrounded,
         ignore_sources=args.ignore_sources)
