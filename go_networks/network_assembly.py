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


# Derived types
DbRefs = Dict[str, str]


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
        # Filled out during cx assembly
        self.entity_lookup: Dict[str, DbRefs] = {}
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
        self.assembled_stmts: List[Statement] = []

    def _add_edge(self):
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

        # Run entity lookup now that statements are assembled
        self._set_entity_lookup()

    def _cx_assembly(self):
        # assemble cx network
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
                f"No nodes in CxNetwork for {self.namespace}:{self.identifier}"
            )

    def _add_metadata(self):
        # Loop edges and add edge and node meta data
        nodes_set = set()
        for e in self.model.edges:
            # Set node meta data
            source = self.model.edges[e]["s"]
            if source not in nodes_set:
                self._add_node_metadata(source)
                nodes_set.add(source)

            target = self.model.edges[e]["t"]
            if target not in nodes_set:
                self._add_node_metadata(target)
                nodes_set.add(target)

            # Set edge metadata
            self._add_edge_metadata(e, source, target)

    def _add_node_metadata(self, node):
        node_name = self.model.nodes[node]["n"]
        db_refs = self.entity_lookup.get(node_name)
        if db_refs is None:
            logger.warning(f"No DB refs for node {node_name}")
            return
        db_refs_list = [_format_entities(ns, _id) for ns, _id in db_refs.items()]
        self.model.set_node_attribute(node, "DB refs", db_refs_list)

    def _add_edge_metadata(self, edge, source, target):
        # ToDo:
        #  Add meta data from pair properties
        #  Get and add sentence from a top ranked statements of the edge
        #  Figure out most specific edge label
        #  - interacts with (undirected)
        #    - affects (for directed statements)
        #      - up/downregulates (for signed-like stmt types)
        #        - stmt type (for specific statements types)
        source_name = (self.model.get_node(source) or {}).get("n")
        target_name = (self.model.get_node(target) or {}).get("n")
        if source_name is None or target_name is None:
            logger.warning(f"Could not get node ids for edge {edge}")
            return
        pair = f"{source_name}|{target_name}"
        rev_pair = f"{target_name}|{source_name}"
        # The edges in the model are from the statements contained in the
        # associated PairProperties, so should not get KeyError here
        try:
            for key, value in (
                self._pp[pair].dict(exclude={"a", "b", "statements"}).items()
            ):
                self.model.set_edge_attribute(edge, key, value)
        # Can happen with Complex
        except KeyError:
            try:
                for key, value in (
                    self._pp[rev_pair].dict(exclude={"a", "b", "statements"}).items()
                ):
                    self.model.set_edge_attribute(edge, key, value)
            except KeyError:
                logger.warning(f"Could not set edge attributes for {pair} "
                               f"or {rev_pair}")
                return

    def assemble(self):
        """Assemble cx network

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
        """
        if self.no_nodes:
            print("Model contains no nodes")
            return

        if self.model is not None:
            print("Model is already assembled!")
            return

        if not self.assembled_stmts:
            logger.info("Running statement preassembly")
            self._stmt_assembly()

        # Assemble cx graph
        self._cx_assembly()

        # Set meta-data
        self._add_metadata()

        # Apply Layout after creation - see format_and_upload_network()
        pass

    def _set_entity_lookup(self):
        # If already set, return
        if self.entity_lookup:
            return

        # Loop agents of statements and create an agent lookup
        for stmt in self.assembled_stmts:
            for ag in stmt.agent_list():
                if ag.name not in self.entity_lookup:
                    self.entity_lookup[ag.name] = ag.db_refs
                elif set(self.entity_lookup[ag.name]).symmetric_difference_update(
                    ag.db_refs
                ):
                    self.entity_lookup[ag.name].update(ag.db_refs)


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


def _format_entities(db_name: str, db_id: str) -> str:
    if ':' in db_id:
        return db_id
    else:
        return f"{db_name}:{db_id}"