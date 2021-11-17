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
import logging
from indra_db_lite.construction import query_to_csv

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
