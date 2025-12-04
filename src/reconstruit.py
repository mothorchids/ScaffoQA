import networkx as nx
from utility import *
import matplotlib.pyplot as plt
import numpy as np
import pulp
from ortools.sat.python import cp_model


def find_biggest_component(G):
    """
    @brief Finds the largest weakly connected component in a directed graph.

    @param G A NetworkX directed graph (`nx.DiGraph` or similar).

    @return A subgraph (`nx.DiGraph`) corresponding to the largest weakly connected component of G.

    This function computes all weakly connected components in the graph `G`,
    selects the largest one, and returns it as a new subgraph.
    """
    list_connected_compossant = sorted(
        nx.weakly_connected_components(G), key=len, reverse=True
    )  # Generate a sorted list of weakly connected components, largest first.

    largest_component = list_connected_compossant[0]
    SG = nx.subgraph(G, largest_component).copy()
    return SG

def combine_DNA(G, s, v, k):
    """
    @brief Combines the DNA sequence of a node with an existing DNA string.

    @param G A NetworkX graph where each node has a 'title' attribute representing a DNA fragment.
    @param s The current assembled DNA sequence (string). Can be empty on the first call.
    @param v The current node whose DNA fragment will be added to the sequence.
    @param k The overlap length (k-mer size).

    @return A new DNA sequence string formed by combining `s` and the fragment at node `v`.

    If the sequences do not overlap properly on `k` characters, an error message is printed.
    """
    sequence_of_v = G.nodes[v]['title']
    if s:
        if s[len(s)-k:len(s)] != sequence_of_v[:k]:
            print(f"Assembly error at node {v}: sequence not matching")
        return s[:len(s)-k] + sequence_of_v
    else: # when s is empty
        return sequence_of_v


def recontruct_DNA(G, path_on_graph, k):
    """
    @brief Reconstructs the DNA sequence from a path of nodes in a graph.

    @param G A NetworkX graph where each node has a 'title' attribute representing a DNA fragment.
    @param path_on_graph A list of nodes representing a path through the graph.
    @param k The overlap length (k-mer size).

    @return A string representing the reconstructed DNA sequence.
    """
    dna_sequence = ""
    for vertex in path_on_graph:
        dna_sequence = combine_DNA(G, dna_sequence, vertex, k)
    return dna_sequence

def find_path_always_first_neighbor(G, starting_point):
    """
    @brief Finds a path in a graph starting from a node, always following the first unvisited successor.

    @param G A NetworkX directed graph.
    @param starting_point The starting node for the path.

    @return A list of nodes representing the traversed path.

    The function stops when it reaches a node with no successors or when a cycle would occur.
    """
    visited = set([starting_point])
    path = [starting_point]
    cat = starting_point  # current location

    while True:
        successors = list(G.successors(cat))
        if not successors:
            break  # no outgoing edges

        a_neighbor = successors[0]  # just take the first one (arbitrarily)
        if a_neighbor in visited:
            break  # avoid cycles
        path.append(a_neighbor)
        visited.add(a_neighbor)
        cat = a_neighbor

    return path

def find_longest_path_dfs(G, starting_point):
    """
    @brief Finds the longest simple path starting from a given node using DFS.
    The length is the number of the nodes in the path. 

    @param G A NetworkX directed graph.
    @param starting_point The node to start the search from.

    @return A list representing the longest path (no repeated nodes).
    """
    def dfs(current, path, visited):
        nonlocal longest_path

        visited.add(current)
        path.append(current)

        extended = False
        for neighbor in G.successors(current):
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


def find_longest_dna_path_dfs(G, starting_point,k):
    """
    @brief Finds the longest simple path starting from a given node using DFS.
    The length is the length of the DNA sequence related to the path. 
    
    @todo this implementation is not optimal. 

    @param G A NetworkX directed graph.
    @param starting_point The node to start the search from.

    @return A list representing the longest path (no repeated nodes).
    """
    def dfs(current, path, visited):
        nonlocal longest_path
        nonlocal longest_length

        visited.add(current)
        path.append(current)

        extended = False
        for neighbor in G.successors(current):
            if neighbor not in visited:
                dfs(neighbor, path, visited)
                extended = True

        if not extended and len(recontruct_DNA(G,path,k)) > longest_length:
            # update
            longest_path = path.copy()
            longest_length = len(recontruct_DNA(G,path,k))

        # Backtrack
        path.pop()
        visited.remove(current)

    longest_length = 0
    longest_path = []
    dfs(starting_point, [], set())
    return longest_path



        