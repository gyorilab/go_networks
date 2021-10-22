import logging
from typing import List, Dict, Optional

from ndex2 import NiceCXNetwork

from indra.ontology.bio import bio_ontology
from indra.preassembler import Preassembler
from indra.preassembler.custom_preassembly import agents_stmt_type_matches
from indra.statements import Statement
from indra.assemblers.cx import NiceCxAssembler
import indra.tools.assemble_corpus as ac
from go_networks.data_models import PairProperty
from indra_db.client import get_curations


logger = logging.getLogger(__file__)
try:
    CURATIONS = get_curations()
except Exception as exc:
    logger.warning("Error trying to get curations")
    logger.warning(exc)
    CURATIONS = []


class CxNoNodes(Exception):
    """Raised for no nodes in graph"""


class GoNetworkAssembler:
    namespace = "GO"

    def __init__(
        self,
        identifier: str,  # go ID
        identifier_description: str,  # Description of GO ID
        entity_list: List[str],  # Associated genes
        pair_properties: Dict[str, PairProperty],  # Lookup for these genes
    ):

        self.identifier = identifier
        self.identifier_description = identifier_description
        self.entity_list = entity_list
        self._pp = pair_properties
        self.model: Optional[NiceCXNetwork] = None  # Save the model here
        self.no_nodes: bool = False
        self.network_attributes = {
            "networkType": "pathway",
            "GO ID": identifier,
            "GO hierarchy": "biological process",
            "Prov:wasGeneratedBy": "INDRA",
            "Organism": "Homo sapiens (Human)",
            "Description": identifier_description,
            "Methods": "This network was assembled automatically by INDRA "
            "(http://indra.bio) by processing all available biomedical "
            "literature with multiple machine reading systems, and "
            "integrating curated pathway databases. The network represents "
            "mechanistic interactions between genes/proteins that are "
            "associated with this GO process.",
        }
        self.assembled_stmts = []

    def _add_node(self):
        # Add db refs
        pass

    def _add_edge(self):
        # Add nodes of edge
        # Add meta data from pair properties
        # Get and add sentence from a top ranked statements of the edge
        # Figure out most specific edge label
        # - interacts with (undirected)
        #   - affects (for directed statements)
        #     - up/downregulates (for signed-like stmt types)
        #       - stmt type (for specific statements types)
        pass

    def _stmt_assembly(self):
        stmts = []
        # Get all statements from pair properties
        for pair_prop in self._pp.values():
            stmts.extend(pair_prop.statements.directed.values())
            for cplx_stmts in pair_prop.statements.undirected.values():
                stmts.extend(cplx_stmts)

        stmts = filter_out_medscan(stmts)

        stmts = ac.filter_gene_list(stmts, self.entity_list, policy="all")
        pa = Preassembler(
            bio_ontology, stmts=stmts, matches_fun=agents_stmt_type_matches
        )
        stmts = pa.combine_duplicates()
        if CURATIONS:
            logger.info("Filtering by curations")
            stmts = ac.filter_by_curation(stmts, curations=CURATIONS)
        self.assembled_stmts = stmts

    def assemble(self):
        """
        Iterate by GO ID and for each list of genes, build a network:
            - Make a node for each gene with metadata representing db_refs
            - Add all the edges between the list of genes with metadata coming
              from the pre-calculated data in (4), i.e. the pair properties
                - When generating links to different statement types, use query
                  links instead of hash-based links e.g., instead of
                  https://db.indra.bio/statements/from_hash/-26724098956735262?format=html,
                  link to
                  https://db.indra.bio/statements/from_agents?subject=HBP1&object=CDKN2A&type=IncreaseAmount&format=html
            - Should we include any sentences?
                - Maybe choose the statement involving A and B with the highest
                  evidence count / belief score and choose one example sentence
                  to add to the edge?
            - How to generate the edge label (the general pattern is
              "A (interaction type) B")?
                - Make these as specific as possible, the least specific being
                  "A (interacts with) B"

        Returns
        -------

        """
        if self.model is not None:
            print("Model is already assembled!")
            return

        if not self.assembled_stmts:
            logger.info("Running preassembly")
            self._stmt_assembly()

        # pre-assemble
        logger.info("Assembling CX model")
        ca = NiceCxAssembler(
            self.assembled_stmts, f"{self.identifier} ({self.identifier_description})"
        )
        ncx = ca.make_model(
            self_loops=False, network_attributes=self.network_attributes
        )
        self.model = ncx
        if not ncx.nodes:
            self.no_nodes = True
            raise CxNoNodes(
                f"No nodes in CxNetwork for {self.namespace}:" f"{self.identifier}"
            )

        # Apply Layout after creation - see format_and_upload_network()
        pass


def filter_out_medscan(stmts: List[Statement]) -> List[Statement]:
    logger.info(f"Filtering out medscan evidence on {len(stmts)} statements")
    filtered_stmts = []
    for stmt in stmts:
        new_evidence = [e for e in stmt.evidence if e.source_api != "medscan"]
        if not new_evidence:
            continue
        stmt.evidence = new_evidence
        if not stmt.evidence:
            continue
        filtered_stmts.append(stmt)
    logger.info("%d statements after filter" % len(filtered_stmts))
    return filtered_stmts
