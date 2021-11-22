"""Get PMC XML from the DB and extract metadata

# DATA COLUMNS USEFUL FOR BIOFACTOID

* journal
* article
* email
* corresponding_author (as Boolean)
* year
* indra_statement_count

# XPATH STATEMENTS

## Authors

.//contrib[@contrib-type="author"]

## Corresponding Author

.//contrib[@contrib-type="author" and @corresp="yes"]

## Author Sub-Properties (Per Author)

For what I have done in the past, I have just taken the first entries

### Last Names

.//surname

### First Name

.//given-names

### Email

.//email
"""
import argparse
import codecs
import logging
from pathlib import Path

from lxml import etree
from tqdm import tqdm

from indra_db_lite.construction import query_to_csv
from gzip import decompress

logger = logging.getLogger(__name__)


def text_ref_xml_to_csv(out_path):
    """Get XML content keyed by text ref id from the DB and save to CSV

    The text_content table has the following columns (copied from
    indra_db/schemas/principal_schema.py):

    - **id** ``integer PRIMARY KEY``: The auto-generated primary key of
      the table. These are elsewhere called Text Content IDs, or TCIDs.
    - **text_ref_id** ``integer NOT NULL``: A foreign-key constrained
      reference to the appropriate entry in the :func:`text_ref <text_ref>`
      table.
    - **source** ``varchar(250) NOT NULL``: The name of the source, e.g.
      "pubmed" or "pmc_oa". The list of content names can be found in the
      class attributes in content managers.
    - **format** ``varchar(250) NOT NULL``: The file format of the
      content, e.g. "XML" or "TEXT".
    - **text_type** ``varchar(250) NOT NULL``: The type of the text, e.g.
      "abstract" of "fulltext".
    - **preprint** ``boolean``: Indicate whether the content is from
      a preprint.
    - **license** [``varchar``]: Record the license that applies to the
      content.
    - **content** ``bytea``: The raw compressed bytes of the content.

    **Metadata Columns**

    - **insert_data** ``timestamp without time zone``: The date the record
      was added.
    - **last_updated** ``timestamp without time zone``: The most recent
      time the record was edited.



    Parameters
    ----------
    out_path :
        Path to the output CSV file.
    """
    query = """
    SELECT
        text_ref_id, encode(content, 'hex')
    FROM 
        text_content
    WHERE
        format = 'xml'
        AND
        text_ref_id IN (
            SELECT
                id
            FROM
                text_ref
            WHERE
                pmcid IS NOT NULL
        )
    """
    query_to_csv(query, out_path)


def text_ref_id_pmc_id_dump(out_path):
    """Get text_ref_id-pmcid pairs from the DB

    The text_ref table has the following columns (copied from
    indra_db/schemas/principal_schema.py):

    - **id** ``integer PRIMARY KEY``: The primary key of the TextRef
      entry. Elsewhere this is often referred to as a "text ref ID" or
      "trid" for short.
    - **pmid** ``varchar(20)``: The identifier from pubmed.
    - **pmcid** ``varchar(20)``: The identifier from PubMed Central (e.g.
      "PMC12345")
    - **doi** ``varchar(100)``: The ideally universal identifier.
    - **pii** ``varchar(250)``: The identifier used by Springer.
    - **url** ``varchar UNIQUE``: For sources found exclusively online
      (e.g. wikipedia) use their URL.
    - **manuscript_id** ``varchar(100) UNIQUE``: The ID assigned documents
      given to PMC author manuscripts.

    **Metadata Columns**

    In addition we also track some basic metadata about the entry and
    updates to the data in the table.

    - **create_date** ``timestamp without time zone``: The date the record
      was added.
    - **last_updated** ``timestamp without time zone``: The most recent
      time the record was edited.
    - **pub_year** ``integer``: The year the article was published, based
      on the first report we find (in order of PubMed, PMC, then PMC
      Manuscripts).

    Parameters
    ----------
    out_path : str
        Path to the output CSV file.
    """
    query = """
    SELECT
        id, pmcid
    FROM 
        text_ref
    WHERE
        pmcid IS NOT NULL
    """
    query_to_csv(query, out_path)


def hex_bin_to_str(raw_hex_bin: str) -> str:
    """Convert hex-encoded raw string to string

    Parameters
    ----------
    raw_hex_bin :
        Hex-encoded bytes as a plain string.

    Returns
    -------
    str
        String
    """
    decode_hex = codecs.getdecoder("hex_codec")
    # It's a plain text string containing hex-encoded bytes, so it's first
    # two characters are escaping the hex-encoding: '\\x1f8b0808......'
    hex_str = decode_hex(raw_hex_bin[2:])[0]
    return decompress(hex_str).decode()


def extract_info_from_pmc_xml(xml_str: str) -> dict:
    """Extract metadata from PMC XML

    * journal
    * article
    * email (from corresponding author, or any email if no corresponding
      author email exists)
    * corresponding_author (as Boolean)
    * year

    Parameters
    ----------
    xml_str :
        PMC XML as a string

    Returns
    -------
    :
    """
    tree = etree.fromstring(xml_str)
    corr_author_query = ".//contrib[@contrib-type='author' and @corresp='yes']"

    def _get_email(root):
        # Get the email of the corresponding author or any email if the
        # corresponding author doesn't have an email
        corr_author = root.xpath(corr_author_query)
        if corr_author and corr_author[0].xpath(".//email"):
            return corr_author[0].xpath(".//email")[0].text
        else:
            any_email = root.xpath(".//email")
            return (any_email[0].text or None) if any_email else None

    def _get_pub_year(root):
        pub_years = [y.text for y in
                     root.xpath(".//pub-date/year")]
        if pub_years:
            # Get earliest publication year
            return min(pub_years)

    # Corresponding author
    corr_auth = tree.xpath(corr_author_query)

    # Get email
    email = _get_email(tree)

    # Journal name
    journal = (tree.xpath(".//journal-title")[0].text or None) if \
        tree.xpath(".//journal-title") else None

    # Article title
    article_title = (tree.xpath(".//article-title")[0].text or None) if \
        tree.xpath(".//article-title") else None

    # Year
    year = _get_pub_year(tree)

    return {
        'journal': journal,
        'article': article_title,
        'email': email,
        'corresponding_author': bool(corr_auth),
        'year': year,
    }


def _read_text_ref_id_pmc_csv(path: str) -> dict:
    with open(path, 'r') as f:
        # text ref ID -> PMC is a many-to-one mapping
        rid_pmc_map = {}
        line = f.readline()
        while line:
            # Assumes the columns are <text ref ID>,<PMC ID>
            rid, pmc = line.strip().split(',')
            rid_pmc_map[rid] = pmc
            line = f.readline()

    return rid_pmc_map


def _read_trid_xml_csv(path: str) -> dict:
    with open(path, 'r') as f:
        trid_xml_map = {}
        line = f.readline()
        while line:
            trid, raw_xml = line.strip().split(',')

            # Convert hex-encoded raw string to string
            xml_str = hex_bin_to_str(raw_xml)

            # Extract metadata from PMC XML
            trid_xml_map[trid] = extract_info_from_pmc_xml(xml_str)

            line = f.readline()

    return trid_xml_map


def main(pmc_reading_id_path: str,
         reading_xml_path: str,
         pmc_count_path: str,
         out_path: str):
    # Get the PMC count
    with open(pmc_count_path, 'r') as f:
        pmc_counts = {}
        line = f.readline()
        while line:
            pmc, count = line.strip().split('\t')
            pmc_counts[pmc] = int(count)
            line = f.readline()

    # Get the reading ID -> PMC mapping
    if not Path(pmc_reading_id_path).exists():
        text_ref_id_pmc_id_dump(pmc_reading_id_path)
    trid_pmc_map = _read_text_ref_id_pmc_csv(pmc_reading_id_path)

    # Get the reading id, XML mapping
    if not Path(reading_xml_path).exists():
        text_ref_xml_to_csv(reading_xml_path)
    rid_xml_info_map = _read_trid_xml_csv(reading_xml_path)

    # Loop the XML info and map the reading id to the PMC id, then write the
    # info to the output TSV with the columns:
    # PMC ID, Journal, Article Title, Corresponding Author, Year, Evidence Count
    processed_ids = set()
    with open(out_path, 'w') as fo:
        for trid, xml_info in tqdm(rid_xml_info_map.items(),
                                   total=len(rid_xml_info_map)):
            pmc = trid_pmc_map.get(trid)

            if pmc is None:
                continue

            assert pmc not in processed_ids

            # Get the evidence count
            count = pmc_counts.get(pmc)
            if count is None:
                continue

            # Write to the output file:
            # * PMC ID
            # * journal
            # * article
            # * email
            # * corresponding_author (as Boolean)
            # * year
            # * indra_statement_count
            fo.write(f'{pmc}\t'
                     f'{xml_info["journal"]}\t'
                     f'{xml_info["article"]}\t'
                     f'{xml_info["email"]}\t'
                     f'{xml_info["corresponding_author"]}\t'
                     f'{xml_info["year"]}\t'
                     f'{count}\n')

            processed_ids.add(pmc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pmc_reading_id_path',
                        help='Path to the PMC reading ID CSV file')
    parser.add_argument('--reading_xml_path',
                        help='Path to the reading id XML CSV file')
    parser.add_argument('--pmc_count_path',
                        help='Path to the PMC count TSV file')
    parser.add_argument('--out_path',
                        help='Path to the output TSV file')
    args = parser.parse_args()

    assert args.pmc_reading_id_path.endswith('.csv'),\
        'PMC reading ID CSV file must be a CSV file'
    assert args.reading_xml_path.endswith('.csv'), \
        'Reading XML CSV file must be a CSV file'
    assert args.pmc_count_path.endswith('.tsv'), \
        'PMC count file must be a TSV file'
    assert args.out_path.endswith('.tsv'), \
        'The output file must be a TSV file'

    main(args.pmc_reading_id_path,
         args.reading_xml_path,
         args.pmc_count_path,
         args.out_path)
