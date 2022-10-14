"""
This file contains the GoNetworkAssembler class, which is used to assemble
an Ndex network from a set of statements related to a given GO term.
"""
import json
import logging
from itertools import combinations
from math import dist, pi, cos, sin
from random import random
from typing import List, Dict, Tuple

import networkx as nx
from networkx.drawing import layout
from numpy import log

from ndex2 import NiceCXNetwork

from indra.assemblers.english import EnglishAssembler
from indra.statements import get_statement_by_name, Agent
from indra.ontology.bio import bio_ontology
from indra.databases import identifiers
from indra.ontology.standardize import standardize_db_refs
from indra.databases import hgnc_client

Coord = Tuple[float, float]
NodeCoords = Dict[str, List[float]]


logger = logging.getLogger(__name__)

context_prefixes = {
    "hgnc": "https://identifiers.org/hgnc:",
    "hgnc.symbol": "https://identifiers.org/hgnc.symbol:",
    "pubmed": "https://identifiers.org/pubmed:",
    "uniprot": "http://identifiers.org/uniprot:",
    "ncbigene": "http://identifiers.org/ncbigene:",
    "mesh": "https://identifiers.org/mesh:",
}


def get_aliases(gene_name):
    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
    db_refs = standardize_db_refs({"HGNC": hgnc_id})
    curies = [
        f"{identifiers.get_identifiers_ns(db_ns)}:{db_id}"
        for db_ns, db_id in db_refs.items()
    ]
    return curies


def _get_english_from_stmt_type(stmt_type: str, source: str, target: str) -> str:
    # Get Statement class from type
    stmt_cls = get_statement_by_name(stmt_type)
    # Get mock statement
    if stmt_cls.__name__ == "Complex":
        stmt = stmt_cls([Agent(source), Agent(target)])
    else:
        stmt = stmt_cls(Agent(source), Agent(target))
    # Get English
    return EnglishAssembler([stmt]).make_model()


def edge_attribute_from_ev_counts(source, target, ev_counts, directed) -> List[str]:
    parts = []
    for stmt_type, cnt in sorted(ev_counts.items(), key=lambda x: x[1], reverse=True):
        english = _get_english_from_stmt_type(stmt_type, source, target)

        # Strip out trailing period
        english = english[:-1] if english[-1] == "." else english
        if directed:
            url = (
                f"https://db.indra.bio/statements/from_agents?"
                f"subject={source}&object={target}&type="
                f"{stmt_type}&format=html&expand_all=true"
            )
        else:
            url = (
                f"https://db.indra.bio/statements/from_agents?"
                f"agent0={source}&agent1={target}&type="
                f"{stmt_type}&format=html&expand_all=true"
            )
        part = f'{english} (<a href="{url}" target="_blank">{cnt}</a>)'
        parts.append(part)
    return parts


def _get_position_on_circle(center: Coord, radius: float):
    # Get a random angle (in radians)
    rad = random() * 2 * pi
    x = radius * cos(rad) + center[0]
    y = radius * sin(rad) + center[1]

    return x, y


def _is_too_close(coord: Coord, node_positions: NodeCoords, min_dist: float):
    # Check if coords are within min_dist of any of the coordinates in the
    # keyed positions
    for other_coord in node_positions.values():
        if dist(coord, other_coord) < min_dist:
            return True

    return False


def _move_nodes_apart(
    pos: NodeCoords, min_dist: float, x_mm: Coord, y_mm: Coord
):
    # Any nodes that are on top of each other should be moved apart, without
    # creating any new overlap

    def _is_within(c: Coord):
        return x_mm[0] < c[0] < x_mm[1] and y_mm[0] < c[1] < y_mm[1]

    # First calculate the euclidian distance between all nodes
    node_pairs = [
        sorted(c, key=lambda x: x[0]) for c in combinations(pos.keys(), 2)
    ]
    distances = {(n1, n2): dist(pos[n1], pos[n2]) for n1, n2 in node_pairs}

    # Loop the pairs and update the position if it's too close
    for n1, n2 in node_pairs:
        if distances[(n1, n2)] < min_dist:
            # Move n1 to a random position a radius == min_dist away from n2
            xy1_new = _get_position_on_circle(pos[n2], min_dist)
            new_min = min_dist

            # While the new position is still within the min-max of all
            # coordinates and there still is overlap, try a new distance and
            # angle
            while _is_within(xy1_new) and _is_too_close(xy1_new, pos, new_min):
                new_min += 0.5*min_dist
                xy1_new = _get_position_on_circle(pos[n2], new_min)

            # Replace the current position of n1
            pos[n1] = list(xy1_new)


def _untangle_layout(g: nx.Graph, pos: NodeCoords):
    # Find the nodes that are disconnected from the rest of the graph and
    # find and fix overlapping nodes
    disconnected_nodes = []
    for node in pos:
        if nx.degree(g, node) == 0:
            disconnected_nodes.append(node)

    # Find min and max x and y values as well as the distance between them
    y_max = max(pos.values(), key=lambda x: x[1])[1]
    y_min = min(pos.values(), key=lambda x: x[1])[1]
    y_dist = y_max - y_min

    x_max = max(pos.values(), key=lambda x: x[0])[0]
    x_min = min(pos.values(), key=lambda x: x[0])[0]
    x_dist = x_max - x_min

    # If there are disconnected nodes, move them to the bottom
    if disconnected_nodes:
        # Sort the nodes lexicographically by name: 0-9, A-Z, a-z
        disconnected_nodes = sorted(sorted(disconnected_nodes), key=str.upper)

        # Move the disconnected nodes to below the graph at 10 % of the
        # y-distance, then set the x position linearly from the min to max with
        # no more than 10 nodes per row with 10 % of the x-distance between rows
        for n, node in enumerate(disconnected_nodes):
            pos[node][0] = x_min + ((n % 10) + 0.5) * (x_dist / 10)
            pos[node][1] = y_min - 0.1 * y_dist * (1 + ((n // 10) % 10))

    # Now check if any nodes are still on top of each other and move them a
    # distance apart if possible
    _move_nodes_apart(pos, min_dist=0.05 * x_dist,
                      x_mm=(x_min, x_max), y_mm=(y_min, y_max))


def get_cx_layout(
    network: NiceCXNetwork, scale_factor: float = 500, untangle: bool = True
) -> Dict:
    """Get a NiceCXNetwork compatible set of coordinates for the nodes

    Returns a dict of coordinates keyed by the node
    """
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


def _get_symb(forward: bool, reverse: bool) -> str:
    if forward and reverse:
        return "<=>"
    elif forward and not reverse:
        return "=>"
    elif not forward and reverse:
        return "<="
    else:
        return "-"


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
        self.go_name = bio_ontology.get_name("GO", identifier)
        self.network_attributes = {
            "networkType": "pathway",
            "GO ID": identifier,
            "GO hierarchy": "biological process",
            "GO term": self.go_name,
            "Prov:wasGeneratedBy": "INDRA",
            "organism": "Human, 9606, H.sapiens",
            "description": "This network was assembled automatically by INDRA "
            "(http://indra.bio) by processing all available biomedical "
            "literature with multiple machine reading systems, and "
            "integrating curated pathway databases. The network represents "
            "mechanistic interactions between genes/proteins that are "
            "associated with this GO process.",
        }
        self.rel_scores = []

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
        self.network.set_network_attribute("name", name)
        self.network.set_network_attribute("@context", json.dumps(context_prefixes))
        for k, v in self.network_attributes.items():
            self.network.set_network_attribute(k, v, "string")

        node_keys = {}
        for gene in self.entity_list:
            node_id = self.network.create_node(
                gene, node_represents=f"hgnc.symbol:{gene}"
            )
            node_keys[gene] = node_id
            aliases = get_aliases(gene)
            self.network.add_node_attribute(
                property_of=node_id,
                name="aliases",
                values=aliases,
                type="list_of_string",
            )
            self.network.add_node_attribute(
                property_of=node_id, name="type", values="protein", type="string"
            )
        for (source, target), (
            forward,
            reverse,
            undirected,
        ) in self.pair_properties.items():
            interaction_symbol = _get_symb(bool(forward), bool(reverse))
            edge_id = self.network.create_edge(
                node_keys[source], node_keys[target], interaction_symbol
            )
            self.network.add_edge_attribute(
                edge_id, "__directed", True if forward else False, "boolean"
            )
            self.network.add_edge_attribute(
                edge_id, "__reverse_directed", True if reverse else False, "boolean"
            )
            # Get lists of english assembled statements with linkouts
            total_ev_count = sum(
                sum(d.values()) for d in (forward, reverse, undirected)
            )
            # Add __evidence_count
            self.network.add_edge_attribute(
                edge_id, "__evidence_count", total_ev_count, "integer"
            )
            # Add __relationship_score as attribute == ln(1+total_ev_count)
            rel_score = 0.3 + log(total_ev_count)
            self.rel_scores.append(rel_score)
            self.network.add_edge_attribute(
                edge_id, "__relationship_score", rel_score, "double"
            )

            all_statements = []
            if forward:
                all_statements.extend(
                    edge_attribute_from_ev_counts(source, target, forward, True)
                )
            if reverse:
                all_statements.extend(
                    edge_attribute_from_ev_counts(target, source, reverse, True)
                )
            if undirected:
                all_statements.extend(
                    edge_attribute_from_ev_counts(source, target, undirected, False)
                )

            # Format the statements to an unordered list:
            html_list = "<ul>"
            for english in all_statements:
                html_list += f"<li>{english}</li>"
            html_list += "</ul>"

            # Add a linkout to all statements involving the pair
            if forward and reverse:
                agent_query = f"agent0={source}&agent1={target}"
            elif forward and not reverse:
                agent_query = f"subject={source}&object={target}"
            elif not forward and reverse:
                agent_query = f"subject={target}&object={source}"
            else:  # not forward and not reverse:
                # this implies there are only undirected statements and
                # subject/object would also include directed statements
                agent_query = f"agent0={source}&agent1={target}"
            url = (
                "https://db.indra.bio/statements/from_agents?"
                f"{agent_query}&format=html&expand_all=true"
            )
            # Edge attribute value should be:
            # View All Evidence (123)
            #   - Statement 1
            #   - Statement 2
            # etc.
            self.network.add_edge_attribute(
                edge_id,
                f"Relationship",
                f'<a href="{url}" target="_blank">View all evidences '
                f"({total_ev_count})</a> {html_list}",
                "string",
            )

        # Get layout by name now that the network is built
        node_layout_by_name = _get_cx_layout(self.network)

        # Loop the entities again to add coordinates for each node
        layout_aspect = []
        for gene in self.entity_list:
            x, y = node_layout_by_name[gene]
            node_id = node_keys[gene]
            layout_aspect.append({"node": node_id, "x": x, "y": y})

        self.network.add_opaque_aspect("cartesianLayout", layout_aspect)
