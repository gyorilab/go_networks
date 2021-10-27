import json
import logging
from typing import List, Dict

from ndex2 import NiceCXNetwork

from indra.ontology.bio import bio_ontology
from indra.databases import identifiers
from indra.ontology.standardize import standardize_db_refs
from indra.databases import hgnc_client


logger = logging.getLogger(__name__)

context_prefixes = {
    'hgnc': 'https://identifiers.org/hgnc:',
    'hgnc.symbol': 'https://identifiers.org/hgnc.symbol:',
    'pubmed': 'https://identifiers.org/pubmed:',
    'uniprot': 'http://identifiers.org/uniprot:',
    'ncbigene': 'http://identifiers.org/ncbigene:',
    'mesh': 'https://identifiers.org/mesh:',
}


def get_aliases(gene_name):
    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
    db_refs = standardize_db_refs({'HGNC': hgnc_id})
    curies = [
        f'{identifiers.get_identifiers_ns(db_ns)}:{db_id}'
        for db_ns, db_id in db_refs.items()
    ]
    return curies


def edge_attribute_from_ev_counts(source, target, ev_counts, directed):
    parts = []
    for stmt_type, cnt in sorted(ev_counts.items(), lambda x: x[1], reverse=True):
        if directed:
            url = f'https://db.indra.bio/statements/from_agents?' \
                f'subject={source}&object={target}&type={stmt_type}&format=html'
        else:
            url = f'https://db.indra.bio/statements/from_agents?' \
                f'agent0={source}&agent1={target}&type={stmt_type}&format=html'
        part = f'{stmt_type}(<a href="{url}" target="_blank">{cnt}</a>)'
        parts.append(part)
    return parts


class GoNetworkAssembler:
    namespace = "GO"

    def __init__(
        self,
        identifier: str,  # go ID
        entity_list: List[str],  # Associated genes
        pair_properties: Dict[str, List[Dict[str, int]]],  # Lookup for these genes
    ):
        self.identifier = identifier
        self.entity_list = entity_list
        self.pair_properties = pair_properties
        self.go_name = bio_ontology.get_name('GO', identifier)
        self.network_attributes = {
            "networkType": "pathway",
            "GO ID": identifier,
            "GO hierarchy": "biological process",
            "Prov:wasGeneratedBy": "INDRA",
            "Organism": "Homo sapiens (Human)",
            "Description": self.go_name,
            "Methods": "This network was assembled automatically by INDRA "
            "(http://indra.bio) by processing all available biomedical "
            "literature with multiple machine reading systems, and "
            "integrating curated pathway databases. The network represents "
            "mechanistic interactions between genes/proteins that are "
            "associated with this GO process.",
        }

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
        name = f"{self.identifier} ({self.go_name})"
        self.network = NiceCXNetwork()
        self.network.set_network_attribute('name', name)
        self.network.set_network_attribute('@context',
            json.dumps(context_prefixes))
        for k, v in self.network_attributes.items():
            self.network.set_network_attribute(k, v, 'string')

        node_keys = {}
        for gene in self.entity_list:
            node_id = self.network.create_node(gene,
                node_represents=f'hgnc.symbol:{gene}')
            node_keys[gene] = node_id
            aliases = get_aliases(gene)
            self.network.add_node_attribute(property_of=node_id,
                                            name='aliases',
                                            values=aliases,
                                            type='list_of_string')
            self.network.add_node_attribute(property_of=node_id,
                                            name='type',
                                            values='protein',
                                            type='string')
        for (source, target), (forward, reverse, undirected) \
                in self.pair_properties.items():
            edge_id = self.network.create_edge(node_keys[source],
                                               node_keys[target],
                                               'interacts with')
            self.network.add_edge_attribute(edge_id, 'directed',
                                            True if forward else False,
                                            'boolean')
            self.network.add_edge_attribute(edge_id, 'reverse directed',
                                            True if reverse else False,
                                            'boolean')
            if forward:
                self.network.add_edge_attribute(edge_id, 'SOURCE => TARGET',
                    edge_attribute_from_ev_counts(source, target, forward,
                                                 True),
                    'list_of_string')
            if reverse:
                self.network.add_edge_attribute(edge_id, 'TARGET => SOURCE',
                    edge_attribute_from_ev_counts(target, source, reverse,
                                                  True),
                    'list_of_string')
            if undirected:
                self.network.add_edge_attribute(edge_id, 'SOURCE - TARGET',
                    edge_attribute_from_ev_counts(source, target, undirected,
                                                  False),
                    'list_of_string')