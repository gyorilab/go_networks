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
import codecs
import logging

from lxml import etree

from indra_db_lite.construction import query_to_csv
from gzip import decompress

logger = logging.getLogger(__name__)


def pmc_xml_to_csv(out_path):
    """Get PMC XML from the DB and extract metadata

    Parameters
    ----------
    out_path : str
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
    """Get text_ref_id and pmcid from the DB

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
    * email
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

    # Authors
    def _extract_authors(root):
        authors = root.xpath(".//contrib[@contrib-type='author']")
        author_dict = {}

        # Loop authors and extract info per author
        for author in authors:
            ad = {}
            # Last name
            ad['last_name'] = author.xpath(".//surname")[0].text
            # First name
            ad['first_name'] = author.xpath(".//given-names")[0].text
            # Email
            ad['email'] = author.xpath(".//email")[0].text
            # Corresponding author
            ad['corresponding_author'] = author.attrib.get('corresp', 'no') == 'yes'
            author_dict[author.attrib.get('id')] = ad
        return author_dict

    # Author info
    author_dict = _extract_authors(tree)

    # Get corresponding author email from the author info dict
    for k, v in author_dict.items():
        if v['corresponding_author']:
            corresponding_author_email = v['email']
            break
    else:
        corresponding_author_email = None

    # Journal name
    journal = tree.xpath(".//journal-title")[0].text

    # Article title
    # FixMe: sometimes the title and the abstract are extracted to the
    #  article-title tag
    article_title = tree.xpath(".//article-title")[0].text

    # Year
    year = tree.xpath(".//pub-date[@pub-type='epub']")[0].xpath(".//year")[0].text

    return {
        'journal': journal,
        'article_title': article_title,
        'corresponding_author_email': corresponding_author_email,
        'year': year,
        'corresponding_author': corresponding_author_email is not None,
    }
