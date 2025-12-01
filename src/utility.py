import argparse
import re, os
from pyvis.network import Network
#import matplotlib.pyplot as plt
import networkx as nx
from networkx.exception import NetworkXNoCycle

##########################################################
##########################################################
### Elementary functions 
##########################################################

def reverse_complement(dna_sequence):
    """
    Return the opposite sign as a string.

    :param sign: A string, either '+' or '-'.
    :return: '-' if input is '+', '+' if input is '-'.
    :raises ValueError: If the input is not '+' or '-'.
    """
    dna_sequence = dna_sequence.strip().upper()

    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    try:
        return ''.join(complement[nuc] for nuc in reversed(dna_sequence))
    except KeyError as e:
        raise ValueError(f"Invalid nucleotide in sequence: {e.args[0]}")


##########################################################
##########################################################
### De Bruijn graph 
##########################################################

def reverse_sign(sign):
    """
    Return the opposite sign as a string.

    :param sign: A string, either '+' or '-'.
    :return: '-' if input is '+', '+' if input is '-'.
    :raises ValueError: If the input is not '+' or '-'.
    """
    if sign == "+":
        return "-"
    elif sign == "-":
        return "+"
    else:
        raise ValueError(f"Invalid sign: {sign}")

def add_node_with_metadata(G, node_id, sequence="", verbose=False):
    """
    Add a node with label and optional title to the graph.

    :param G: The De Bruijn graph (NetworkX MultiDiGraph).
    :param node_id: Node identifier.
    :param sequence: Optional sequence string to show as tooltip.
    :param verbose: Whether to print debug info.
    """
    if not G.has_node(node_id):
        G.add_node(node_id, label=node_id, title=sequence)
        if verbose:
            print(f"[INFO] Node '{node_id}' created.")
    elif verbose:
        print(f"[SKIP] Node '{node_id}' already exists.")

def edge_with_sign_exists(G, u, v, sign_begin, sign_end):
    """
    Check if an edge from u to v exists in the graph with matching signs.

    :param G: The De Bruijn graph (NetworkX MultiDiGraph).
    :param u: Source node.
    :param v: Target node.
    :param sign_begin: Expected sign at the start of the edge (e.g., '+' or '-').
    :param sign_end: Expected sign at the end of the edge.
    :return: True if such an edge exists, False otherwise.
    """
    if not G.has_edge(u, v):
        return False
    # --- Multigraphs ------------------------------------------------------
    # edge_data is a dict keyed by edge keys → attribute dicts
    for _, edge_attrs in G.get_edge_data(u, v).items(): # Multigraph formulation
        if (
            edge_attrs.get('sign_begin') == sign_begin and
            edge_attrs.get('sign_end') == sign_end
        ):
            return True
    return False

def parse_de_bruijn_link(link):
    """
    Validate and parse a De Bruijn link string.

    :param link: A string of the format 'L:<sign1>:<target>:<sign2>'.
    :return: Tuple containing (sign1, target, sign2).
    :raises ValueError: If the link is malformed.
    """
    if not link:
        raise ValueError("Empty link string")

    parts = link.strip().split(":")
    if len(parts) != 4 or parts[0] != "L":
        raise ValueError(f"Invalid link format: {link}")

    _, sign1, target, sign2 = parts
    if sign1 not in {"+", "-"} or sign2 not in {"+", "-"}:
        raise ValueError(f"Invalid signs in link: {link}")

    return sign1, target, sign2

def add_de_bruijn_edge(G, v1, v2, s1, s2, verbose=False):
    """
    Add a pair of directed edges with opposite signs between two nodes.

    Adds two directed edges to the De Bruijn graph G:
    - One from v1 to v2 with signs (s1, s2).
    - One from v2 to v1 with reversed signs.
    
    :param G: The De Bruijn graph (NetworkX MultiDiGraph).
    :param v1: The source node of the first edge.
    :param v2: The target node of the first edge.
    :param s1: The sign at the start of the v1 → v2 edge ('+' or '-').
    :param s2: The sign at the end of the v1 → v2 edge ('+' or '-').
    :param verbose: If True, print info and warnings.
    """
    # Optionally validate signs
    if s1 not in {"+", "-"} or s2 not in {"+", "-"}:
        raise ValueError(f"Invalid signs: s1={s1}, s2={s2}. Must be '+' or '-'.")
    
    if not edge_with_sign_exists(G, v1, v2, s1, s2):
        G.add_edge(v1, v2, sign_begin=s1, sign_end=s2, label=s1+">"+s2)
        G.add_edge(v2, v1, sign_begin=reverse_sign(s2), sign_end=reverse_sign(s1), label=reverse_sign(s2)+">"+reverse_sign(s1))
        if verbose:
            print(f"Added edges between {v1} and {v2} with signs {s1}, {s2}")



def add_de_bruijn_edges(G, v, links, verbose=False):
    """
    Parse and add multiple De Bruijn edges from a list of link strings.

    :param G: The De Bruijn graph (NetworkX MultiDiGraph).
    :param v: The source node.
    :param links: A list of strings, each formatted as 'L:<sign1>:<target>:<sign2>'.
    :param verbose: If True, print info and warnings.
    """
    for link in links:
        try:
            sign1, target, sign2 = parse_de_bruijn_link(link)
            add_de_bruijn_edge(G, v, target, sign1, sign2, verbose)

        except ValueError as e:
            print(f"Skipping invalid link '{link}': {e}")

def generate_de_bruijn_graph_from_fa(path, verbose=False):
    """
    Generate a directed De Bruijn MultiDiGraph from a sequence-based file.

    Parses a file (e.g., GFA-like) where each node is defined by a line starting with '>',
    followed by a sequence line. Edges are described using 'L:' formatted links.

    :param path: Path to the input file.
    :param verbose: If True, prints progress and debug information.
    :return: A NetworkX MultiDiGraph with nodes and signed De Bruijn edges.
    """
    with open(path) as f:
        lines = f.read().splitlines()

    G = nx.MultiDiGraph()
    entries = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith(">"):
            parts = line.split()
            l_elements = [p for p in parts if p.startswith("L:")]

            match = re.search(r"LN:[^\d]*(\d+)", line)
            if not match:
                if verbose:
                    print(f"[WARN] 'LN' not found in header line: {line}")
                i += 1
                continue

            node_val = parts[0][1:]  # Remove leading '>'
            sequence = lines[i + 1].strip() if i + 1 < len(lines) else ""

            G.add_node(node_val, label=node_val, title=sequence)
            if verbose:
                print(f"[INFO] Node {node_val} created")

            entries.append({
                "val": node_val,
                "L_elements": l_elements,
                "sequence": sequence
            })

            i += 2
        else:
            i += 1

    for entry in entries:
        add_de_bruijn_edges(G, entry["val"], entry["L_elements"], verbose=verbose)

    return G

##########################################################
##########################################################
### General function for graphs
##########################################################
def save_graph_to_html(graph, filename, width='2000px', height='2000px', directed=True, notebook=False):
    """
    Save a NetworkX graph as an interactive HTML visualization using pyvis.

    :param G: A NetworkX graph (e.g., nx.Graph, nx.DiGraph, nx.MultiDiGraph).
    :param filename: Base filename (without extension) to save the HTML file.
    :param width: Width of the visualization (default '2000px').
    :param height: Height of the visualization (default '2000px').
    :param directed: Whether the graph is directed (default True).
    :param notebook: Whether to enable Jupyter notebook mode (default False).
    """
    nt = Network(directed=directed, notebook=notebook)
    #nt = Network(width, height, directed=directed, notebook=notebook)

    # Populate pyvis network from NetworkX graph (note: multiedges not supported)
    nt.from_nx(graph)

    # Optional: show physics options in the browser
    nt.show_buttons(filter_=['physics'])

    # Write to HTML file
    nt.write_html(f"{filename}.html")

def remove_isolated_nodes(G, verbose=False):
    """
    Remove all isolated nodes (degree 0).
    
    :param G: A NetworkX graph (e.g., nx.Graph, nx.DiGraph, nx.MultiDiGraph).
    :param verbose: If True, prints the number and list of removed nodes.
    :return: None. The graph is modified in place.
    """
    isolated_nodes = [n for n in G.nodes() if G.degree(n) == 0]
    G.remove_nodes_from(isolated_nodes)
    if verbose:
        print(f"[INFO] Removed {len(isolated_nodes)} isolated node(s): {isolated_nodes}")

def convert_de_bruijn_to_digraph(G, verbose=False):
    """
    Convert a De Bruijn MultiDiGraph into a simplified directed graph with strand orientation.

    For each node in the graph, this function adds both the original node and its reverse complement.
    Then, for each edge, it uses the 'sign_begin' and 'sign_end' attributes to determine whether
    to use the normal or complement form of the source and destination nodes.

    Nodes that remain isolated in the final graph are removed.

    :param G: The De Bruijn graph (NetworkX MultiDiGraph).
    :param verbose: If True, prints debug messages.
    :return: A directional graph (NetworkX DiGraph).
    """
    DG = nx.DiGraph()
    DG.add_nodes_from(G.nodes(data=True))  # Copy original nodes with attributes

    # Add complementary nodes
    for node in G.nodes():
        val = G.nodes[node].get("val", str(node))
        seq = G.nodes[node].get("title", "")
        comp_node = "c" + val
        comp_seq = reverse_complement(seq)
        DG.add_node(comp_node, label=comp_node, title=comp_seq)
        if verbose:
            print(f"[INFO] Complement node '{comp_node}' added.")

    # Add edges based on sign attributes
    for u, v, attrs in G.edges(data=True):
        sign_begin = attrs.get("sign_begin")
        sign_end = attrs.get("sign_end")

        if sign_begin not in ('+', '-') or sign_end not in ('+', '-'):
            if verbose:
                print(f"[WARN] Skipping edge ({u} → {v}) due to missing/invalid signs.")
            continue

        # Determine actual source and destination
        from_node = u if sign_begin == '+' else "c" + str(u)
        to_node = v if sign_end == '+' else "c" + str(v)

        DG.add_edge(from_node, to_node)
        if verbose:
            print(f"[INFO] Edge added: {from_node} → {to_node}")

    # Remove isolated nodes
    remove_isolated_nodes(DG, verbose)

    return DG

def find_and_print_cycles(G, verbose=True):
    """
    Finds and prints all simple cycles in the directed graph G.

    :param G: A NetworkX directed graph (e.g., nx.DiGraph, nx.MultiDiGraph).
    :param verbose: If True, print info and warnings.
    """
    cycles = list(nx.simple_cycles(G))
    if not cycles:
        if verbose:
            print("No cycles found.")
        return
    if verbose:
        for i, cycle in enumerate(cycles, 1):
            print(f"Cycle {i} found: {cycle}")

def print_connected_components_info(G):
    """
    Prints the number of weakly and strongly connected components of a directed graph G.
    
    :param G: A NetworkX directed graph (e.g., nx.DiGraph, nx.MultiDiGraph).
    """

    num_weak = nx.number_weakly_connected_components(G)
    num_strong = nx.number_strongly_connected_components(G)
    print(f"Number of weakly connected components: {num_weak}")
    print(f"Number of strongly connected components: {num_strong}")

def find_biggest_component(G):
    """
    Finds the largest weakly connected component in a directed graph.

    :param G: A NetworkX directed graph (`nx.DiGraph` or similar).
    :return: A subgraph (`nx.DiGraph`) corresponding to the largest weakly connected component of G.

    This function computes all weakly connected components in the graph `G`,
    selects the largest one, and returns it as a new subgraph.
    """
    list_connected_compossant = sorted(
        nx.weakly_connected_components(G), key=len, reverse=True
    )  # Generate a sorted list of weakly connected components, largest first.

    largest_component = list_connected_compossant[0]
    SG = nx.subgraph(G, largest_component).copy()
    return SG

def k_hop_neighborhood(G, node, k=1, include_self=True):
    """
    Returns the subgraph containing all nodes within k hops of a given node.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    node : node in G
        The central vertex.
    k : int, optional
        Number of layers/hops to include (default is 1).
    include_self : bool, optional
        Whether to include the central node itself (default is True).

    Returns
    -------
    H : networkx.Graph
        The induced subgraph of the k-hop neighborhood.
    """
    # Use BFS to get all nodes within k hops
    nodes_within_k = nx.single_source_shortest_path_length(G, node, cutoff=k).keys()
    
    if not include_self:
        nodes_within_k = [n for n in nodes_within_k if n != node]
    
    H = G.subgraph(nodes_within_k).copy()
    return H


##########################################################
##########################################################
### Main
##########################################################
       
def main():
    parser = argparse.ArgumentParser(description="Demo of De Bruijn graph construction and its conversion to a simplified directed graph.")
    parser.add_argument("input_path", help="Path to the input file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--info", action="store_true", help="Enable verbose output")
    args = parser.parse_args()
    path = args.input_path
    verbose = args.verbose
    filename = os.path.splitext(os.path.basename(path))[0]

    if verbose:
        print(f"[INFO] Generating De Bruijn graph from: {path}")
    G = generate_de_bruijn_graph_from_fa(path,verbose)
    save_graph_to_html(G,filename+"_de-Bruijn")
    remove_isolated_nodes(G,verbose)
    if args.info:
        print("Information about the de Bruijn graph:")
        #find_and_print_cycles(G)
        print_connected_components_info(G)
    if verbose:
        print("[INFO] Converting De Bruijn graph to simplified directed graph")
    DG = convert_de_bruijn_to_digraph(G,verbose)
    if args.info:
        print("Information about the directed graph:")
        #find_and_print_cycles(DG)
        print_connected_components_info(DG)
    save_graph_to_html(DG,filename+"_directed")
    
if __name__ == "__main__":
    main()