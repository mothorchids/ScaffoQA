"""
graph_path_problem.py

Defines the GraphPathProblem class to solve pathfinding tasks on directed graphs.
The graph is derived from a typed De Bruijn graph, and a DNA sequence can be reconstructed by traversing a path within the graph.

Author : Jui-Ting Lu
"""
import os
import matplotlib.pyplot as plt
import networkx as nx
from bidict import bidict

from utility import *

class GraphPathProblem:
    """
    A class to solve pathfinding tasks on directed graphs.
    """

    def __init__(self, graph=nx.DiGraph, k_mer=0):
        """
        Initialize the graph path problem.

        :param graph: A directed graph (NetworkX DiGraph).
        :param k_mer: The length of substrings representing the overlap between two nodes.
        :param edges_register: A bidirectional dictionary (bidict) that maps the edge notation (i,j) in the graph 
            to its corresponding variable index k in the QUBO matrix, and vice versa.
        """
        self.graph = graph
        self.k_mer = k_mer
        self.edges_register = bidict()
    
    def update_edges_register(self):
        er = bidict()
        for i, edge in enumerate(self.graph.edges()):
            er[edge] = i
        self.edges_register = er

    def draw_graph(self, colors=None, pos=None):
        colors = ["green" for _ in self.graph.nodes()]
        pos = nx.spring_layout(self.graph)
        default_axes = plt.axes(frameon=True)
        nx.draw_networkx(self.graph, node_color=colors, node_size=600, alpha=0.8, ax=default_axes, pos=pos)
        edge_labels = nx.get_edge_attributes(self.graph, "weight")
        nx.draw_networkx_edge_labels(self.graph, pos=pos, edge_labels=edge_labels)

    def combine_DNA(self, s:str, v:list) -> str:
        """
        Combines the DNA sequence of a node with an existing DNA string.

        Appends the DNA fragment from the given node to the existing sequence, taking
        into account the k-mer overlap length. If the sequences do not overlap
        correctly, an error message is printed.

        :param s: Current assembled DNA sequence. Can be empty on the first call.
        :param v: Current node whose DNA fragment will be appended.
        :return: New DNA sequence obtained by combining ``s`` with the fragment at node ``v``.
        """
        sequence_of_v = self.graph.nodes[v]['title']
        if s:
            if s[len(s)-self.k_mer:len(s)] != sequence_of_v[:self.k_mer]:
                print(f"Assembly error at node {v}: sequence not matching")
            return s[:len(s)-self.k_mer] + sequence_of_v
        else: # when s is empty
            return sequence_of_v

    def recontruct_DNA(self, path_on_graph:list)->str:
        """
        Reconstruct the DNA sequence from a path of nodes in a graph.

        :param path_on_graph: A list of nodes representing a path through the graph.

        :return: The reconstructed DNA sequence.
        """
        dna_sequence = ""
        for vertex in path_on_graph:
            dna_sequence = self.combine_DNA(dna_sequence, vertex)
        return dna_sequence
    
    def show_highlight_edges(self, highlight_edges):
        """
        Create and visualize a graph, highlighting specific edges.

        :param highlight_edges: List of edges to highlight, each represented as a pair of nodes.
        """
        colors = ["r" for node in self.graph.nodes()]
        pos = nx.spring_layout(self.graph)
        default_axes = plt.axes(frameon=True)
        nx.draw_networkx(self.graph, node_color=colors, node_size=600, alpha=0.8, ax=default_axes, pos=pos)
        edge_labels = nx.get_edge_attributes(self.graph, "weight")
        # Draw all edges with default style
        nx.draw_networkx_edges(self.graph, pos, edgelist=self.graph.edges(), edge_color='gray')
        # Draw highlighted edges with different style
        nx.draw_networkx_edges(self.graph, pos, edgelist=highlight_edges, edge_color='green', width=2)

    def show_solution(self,solution):
        """
        Display or process the given solution.

        :param solution: Usually a bit string representing the solution.
        :type solution: str
        """
        highlight_edges = []
        for i in range(len(solution)):
            if int(solution[i]) == 1:
                highlight_edges.append(self.edges_register.inverse[i])
        self.show_highlight_edges(highlight_edges)

    def find_path_always_first_neighbor(self, starting_point):
        """
        Finds a path in a graph starting from a given node, always following the first unvisited successor.

        The traversal stops when a node has no successors or when continuing would create a cycle.

        :param G: NetworkX directed graph.
        :param starting_point: Node from which the path traversal begins.
        :return: List of nodes representing the traversed path.
        """
        visited = set([starting_point])
        path = [starting_point]
        cat = starting_point  # current location

        while True:
            successors = list(self.graph.successors(cat))
            if not successors:
                break  # no outgoing edges

            a_neighbor = successors[0]  # just take the first one (arbitrarily)
            if a_neighbor in visited:
                break  # avoid cycles
            path.append(a_neighbor)
            visited.add(a_neighbor)
            cat = a_neighbor

        return path

    def find_longest_path_dfs(self, starting_point):
        """
        Finds the longest simple path starting from a given node using DFS.

        The path length is defined as the number of nodes in the path.

        :param starting_point: Node from which the search begins.
        :return: List of nodes representing the longest path without repetitions.
        """
        def dfs(current, path, visited):
            nonlocal longest_path

            visited.add(current)
            path.append(current)

            extended = False
            for neighbor in self.graph.successors(current):
                if neighbor not in visited:
                    dfs(neighbor, path, visited)
                    extended = True

            if not extended and len(path) > len(longest_path):
                longest_path = path.copy()

            # Backtrack
            path.pop()
            visited.remove(current)

        longest_path = []
        dfs(starting_point, [], set())
        return longest_path
    
    def find_longest_dna_path_dfs_from(self, starting_point):
        """
        Finds the longest simple path starting from a given node using DFS.

        The path length is measured as the total length of the DNA sequences associated 
        with the nodes in the path.

        .. todo:: This implementation is not optimal.

        :param G: NetworkX directed graph.
        :param starting_point: Node from which the search begins.
        :return: List of nodes representing the longest path without repetitions.
        """
        def dfs(current, path, visited):
            nonlocal longest_path
            nonlocal longest_length

            visited.add(current)
            path.append(current)

            extended = False
            for neighbor in self.graph.successors(current):
                if neighbor not in visited:
                    dfs(neighbor, path, visited)
                    extended = True

            if not extended and len(self.recontruct_DNA(self.graph,path,self.k_mer)) > longest_length:
                # update
                longest_path = path.copy()
                longest_length = len(self.recontruct_DNA(self.graph,path,self.k_mer))

            # Backtrack
            path.pop()
            visited.remove(current)

        longest_length = 0
        longest_path = []
        dfs(starting_point, [], set())
        return longest_path
    
    def find_longest_dna_path_dfs_to(self, finish_point):
        """
        Finds the longest simple path starting from a given node using DFS.

        The path length is measured as the total length of the DNA sequences associated 
        with the nodes in the path.

        .. todo:: This implementation is not optimal.

        :param G: NetworkX directed graph.
        :param starting_point: Node from which the search begins.
        :return: List of nodes representing the longest path without repetitions.
        """
        def dfs(current, path, visited):
            nonlocal longest_path
            nonlocal longest_length

            visited.add(current)
            path.append(current)

            extended = False
            for neighbor in self.graph.predecessors(current):
                if neighbor not in visited:
                    dfs(neighbor, path, visited)
                    extended = True

            if not extended and len(self.recontruct_DNA(self.graph,path,self.k_mer)) > longest_length:
                # update
                longest_path = path.copy()
                longest_length = len(self.recontruct_DNA(self.graph,path,self.k_mer))

            # Backtrack
            path.pop()
            visited.remove(current)

        longest_length = 0
        longest_path = []
        dfs(finish_point, [], set())
        return longest_path
    

##########################################################
##########################################################
### MAIN
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

    G = generate_de_bruijn_graph_from_fa(path,verbose)
    save_graph_to_html(G,filename+"_de-Bruijn")
    remove_isolated_nodes(G,verbose)
    DG = convert_de_bruijn_to_digraph(G,verbose)

    SG = find_biggest_component(DG)
  
if __name__ == "__main__":
    main()