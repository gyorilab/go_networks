import json
import logging
from typing import List, Dict

import networkx as nx
from networkx.drawing import layout

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
    for stmt_type, cnt in sorted(ev_counts.items(), key=lambda x: x[1], reverse=True):
        if directed:
            url = f'https://db.indra.bio/statements/from_agents?' \
                f'subject={source}&object={target}&type={stmt_type}&format=html'
        else:
            url = f'https://db.indra.bio/statements/from_agents?' \
                f'agent0={source}&agent1={target}&type={stmt_type}&format=html'
        part = f'{stmt_type}(<a href="{url}" target="_blank">View {cnt} statements</a>)'
        parts.append(part)
    return parts


def _untangle_layout(g: nx.Graph,
                     pos: Dict[str, List[float]]):
    # Find the nodes that are disconnected from the rest of the graph
    disconnected_nodes = []
    for node in pos:
        if nx.degree(g, node) == 0:
            disconnected_nodes.append(node)

    # If there are no disconnected nodes, return
    if not disconnected_nodes:
        return

    # Find min and max x and y values as well as the distance between them
    y_max = max(pos.values(), key=lambda x: x[1])[1]
    y_min = min(pos.values(), key=lambda x: x[1])[1]
    y_dist = y_max - y_min

    x_max = max(pos.values(), key=lambda x: x[0])[0]
    x_min = min(pos.values(), key=lambda x: x[0])[0]
    x_dist = x_max - x_min

    # Sort the nodes lexicographically by name: 0-9, A-Z, a-z
    disconnected_nodes = sorted(sorted(disconnected_nodes), key=str.upper)

    # Move the disconnected nodes to below the graph at 10 % of the
    # y-distance, then set the x position linearly from the min to max with
    # no more than 10 nodes per row with 10 % of the x-distance between rows
    for n, node in enumerate(disconnected_nodes):
        pos[node][0] = x_min + ((n % 10) + 0.5) * (x_dist / 10)
        pos[node][1] = y_min - 0.1 * y_dist * (1 + ((n//10) % 10))


def _get_cx_layout(network: NiceCXNetwork, scale_factor: float = 500,
                   untangle: bool = True) -> Dict:
    # Convert to networkx graph
    g = network.to_networkx()

    # Get layout
    pos = layout.kamada_kawai_layout(g, scale=scale_factor)

    # Untangle nodes
    if untangle:
        _untangle_layout(g, pos)

    # Convert to cx layout: need to invert the y-axis
    cx_pos = {node: [x, -y] for node, (x, y) in pos.items()}
    return cx_pos


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

        # Get layout by name now that the network is built
        node_layout_by_name = _get_cx_layout(self.network)

        # Loop the entities again to add coordinates for each node
        layout_aspect = []
        for gene in self.entity_list:
            x, y = node_layout_by_name[gene]
            node_id = node_keys[gene]
            layout_aspect.append({'node': node_id, 'x': x, 'y': y})

        self.network.add_opaque_aspect('cartesianLayout', layout_aspect)
