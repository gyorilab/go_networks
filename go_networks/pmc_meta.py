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
import csv
import argparse
import codecs
import logging
from pathlib import Path
from typing import Optional

from lxml import etree
from tqdm import tqdm

from indra_db_lite.construction import query_to_csv
from gzip import decompress

from indra.literature.pmc_client import _remove_elements_by_tag, \
    _replace_unwanted_elements_with_their_captions, _select_from_top_level

logger = logging.getLogger(__name__)


def buf_count_newlines_gen(fname: str) -> int:
    # Source https://stackoverflow.com/a/68385697/10478812
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b:
                break
            yield b

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count


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
    start_ix = 2 if raw_hex_bin.startswith("\\") else 0
    hex_str = decode_hex(raw_hex_bin[start_ix:])[0]
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
    tree = etree.fromstring(xml_str.encode('utf-8'))
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
        pub_years = [y.text for y in root.xpath(".//pub-date/year")]
        if pub_years:
            # Get earliest publication year
            return min(pub_years)

    def _get_title(root):
        # Remove namespaces if any exist
        if root.tag.startswith('{'):
            for element in root.getiterator():
                # The following code will throw a ValueError for some
                # exceptional tags such as comments and processing instructions.
                # It's safe to just leave these tag names unchanged.
                try:
                    element.tag = etree.QName(element).localname
                except ValueError:
                    continue
            etree.cleanup_namespaces(root)
        # Strip out latex
        _remove_elements_by_tag(root, 'tex-math')
        # Strip out all content in unwanted elements except the captions
        _replace_unwanted_elements_with_their_captions(root)
        # First process front element. Titles alt-titles and abstracts
        # are pulled from here.
        front_elements = _select_from_top_level(root, 'front')
        title_xpath = './article-meta/title-group/article-title'
        for front_element in front_elements:
            for element in front_element.xpath(title_xpath):
                return ' '.join(element.itertext())

    # Corresponding author
    corr_auth = tree.xpath(corr_author_query)

    # Get email
    email = _get_email(tree)

    # Journal name
    journal = (
        (tree.xpath(".//journal-title")[0].text or None)
        if tree.xpath(".//journal-title")
        else None
    )

    # Article title
    try:
        article_title = _get_title(tree)
    except ValueError:
        article_title = None

    # Year
    year = _get_pub_year(tree)

    return {
        "journal": journal,
        "article": article_title,
        "email": email,
        "corresponding_author": bool(corr_auth),
        "year": year,
    }


def _read_text_ref_id_pmc_csv(path: str) -> dict:
    logger.info(f"Reading text ref ID - PMC CSV: {path}")
    with open(path, "r") as f:
        # text ref ID -> PMC is a many-to-one mapping???
        trid_pmc_map = {}
        line = f.readline()
        # Check if the first line is a header
        if line.startswith(("trid,", "text_ref_id,")):
            line = f.readline()
        while line:
            # Assumes the columns are <text ref ID>,<PMC ID>
            trid, pmc = line.strip().split(",")
            trid_pmc_map[trid] = pmc
            line = f.readline()

    return trid_pmc_map


def _read_trid_xml_csv(path: str) -> dict:
    logger.info(f"Reading TRID - XML CSV: {path}")
    with open(path, "r") as f:
        trid_xml_map = {}
        line = f.readline()
        # Check if the first line is the header
        if line.startswith(("text_ref_id,", "trid,")):
            line = f.readline()
        while line:
            trid, raw_xml = line.strip().split(",")

            # Convert hex-encoded raw string to string
            xml_str = hex_bin_to_str(raw_xml)

            # Extract metadata from PMC XML
            trid_xml_map[trid] = extract_info_from_pmc_xml(xml_str)

            line = f.readline()

    return trid_xml_map


def main(
        pmc_reading_id_path: str,
        reading_xml_path: str,
        pmc_count_path: str,
        out_path: str,
        xml_lines: Optional[int] = None,
):
    # Get the PMC count
    with open(pmc_count_path, "r") as f:
        logger.info(f"Reading PMC count file: {pmc_count_path}")
        pmc_counts = {}
        line = f.readline()
        while line:
            pmc, count = line.strip().split("\t")
            pmc_counts[pmc] = int(count)
            line = f.readline()

    # Get the reading ID -> PMC mapping
    if not Path(pmc_reading_id_path).exists():
        text_ref_id_pmc_id_dump(pmc_reading_id_path)
    trid_pmc_map = _read_text_ref_id_pmc_csv(pmc_reading_id_path)

    # Get the reading id, XML mapping
    if not Path(reading_xml_path).exists():
        text_ref_xml_to_csv(reading_xml_path)

    # Loop the XML info and map the reading id to the PMC id, then write the
    # info to the output TSV with the columns:
    # PMC ID, Journal, Article Title, Corresponding Author, Year, Evidence Count
    processed_ids = set()
    missing_pmc_mapping = 0
    missing_counts = 0
    duplicate_pmc_mapping = 0
    if xml_lines is None:
        xml_lines = buf_count_newlines_gen(reading_xml_path)

    logger.info(f"Processing {xml_lines} lines from XML {reading_xml_path}.")
    t = tqdm(total=xml_lines)
    with open(reading_xml_path, "r") as fi, \
            open(out_path, "w", newline='') as fo, \
            open('failed_xml.csv', 'w') as f_failed:
        # Get csv writer
        writer = csv.writer(fo, delimiter="\t")

        # Add header to output file
        writer.writerow(["pmc_id", "journal", "article_title", "email",
                         "corresponding_author", "year", "evidence_count"])
        line = fi.readline()
        read_lines = 1
        while line:
            # Update progress bar
            t.update()

            # Get content
            trid, raw_xml = line.strip().split(",")

            # Get PMC ID
            pmc = trid_pmc_map.get(trid)

            # Skip if no PMC ID
            if pmc is None:
                missing_pmc_mapping += 1
                line = fi.readline()
                continue

            # Skip if already processed
            elif pmc in processed_ids:
                duplicate_pmc_mapping += 1
                line = fi.readline()
                continue

            # Get the evidence count
            count = pmc_counts.get(pmc, 0)

            # Count the no counts, but don't skip
            if not count:
                missing_counts += 1

            try:
                # Convert hex-encoded raw string to string and extract metadata
                # from PMC XML
                xml_str = hex_bin_to_str(raw_xml)
                xml_info = extract_info_from_pmc_xml(xml_str)
            except ValueError:
                logger.warning(
                    f"Failed to parse XML for PMC {pmc}; TRID {trid}")
                f_failed.write(f"{pmc},{trid}\n")
                line = fi.readline()
                continue

            # Write to the output file:
            # * PMC ID
            # * journal
            # * article
            # * email
            # * corresponding_author (as Boolean)
            # * year
            # * indra_statement_count
            writer.writerow([
                pmc,
                xml_info["journal"],
                xml_info["article"],
                xml_info["email"],
                xml_info["corresponding_author"],
                xml_info["year"],
                count
            ])

            processed_ids.add(pmc)

            if read_lines > xml_lines:
                logger.info(f"Read {read_lines} lines")
                break

            line = fi.readline()
            read_lines += 1

    t.close()

    if missing_pmc_mapping:
        logger.info(f"{missing_pmc_mapping} missing PMC mapping(s)")

    if missing_counts:
        logger.info(f"{missing_counts} missing PMC count(s)")

    if duplicate_pmc_mapping:
        logger.info(f"{duplicate_pmc_mapping} duplicate PMC mapping(s)")

    logger.info(f"Wrote PMC IDs that failed to process to "
                f"{Path('./failed_xml.csv').absolute().as_posix()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pmc_reading_id_path", help="Path to the PMC reading ID CSV file"
    )
    parser.add_argument(
        "--reading_xml_path", help="Path to the reading id XML CSV file"
    )
    parser.add_argument("--pmc_count_path",
                        help="Path to the PMC count TSV file")
    parser.add_argument("--out_path", help="Path to the output TSV file")
    parser.add_argument(
        "--xml_lines", type=int,
        help="Number of lines to read from the XML file. If "
             "not specified, all lines will be read.",
    )
    args = parser.parse_args()

    assert args.pmc_reading_id_path.endswith(
        ".csv"
    ), "PMC reading ID CSV file must be a CSV file"
    assert args.reading_xml_path.endswith(
        ".csv"
    ), "Reading XML CSV file must be a CSV file"
    assert args.pmc_count_path.endswith(
        ".tsv"), "PMC count file must be a TSV file"
    assert args.out_path.endswith(".tsv"), "The output file must be a TSV file"

    main(
        args.pmc_reading_id_path,
        args.reading_xml_path,
        args.pmc_count_path,
        args.out_path,
        args.xml_lines,
    )
