import networkx as nx
from utility import *
from reconstruit import *
import matplotlib.pyplot as plt
import numpy as np
import pulp
from ortools.sat.python import cp_model
import matplotlib.colors as mcolors

def generate_colors(n):
    """Génère n couleurs distinctes en format hexadécimal"""
    cmap = plt.colormaps.get_cmap('hsv')  # on charge juste la palette
    colors = [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, n)]
    return colors

def registre_nodes(G):
    """
    Creation of registre thats maps nodes of G onto liste [1,2,3,...]
    """ 
    registre_nodes={}
    num=0
    for node in G.nodes():
        registre_nodes[node]=num
        num+=1

    return registre_nodes

def registre_edges(G):
    """
    Creation of registre thats maps edges of G onto liste [1,2,3,...]
    """
    registre_edges={}
    num=0
    for edge in G.edges():
        registre_edges[edge]=num
        num+=1

    return registre_edges

def registre_pred_node(G):
    """
    Creation of registre that contains the predecessors of nodes
    """

    registre_predecessors={}
    for node in list( G.nodes() ):
        registre_predecessors[node]=list( G.predecessors(node) )

    return registre_predecessors

def registre_suc_node(G):
    """
    Creation of another registre that contains the successors of nodes
    """
    registre_successors={}
    for node in list( G.nodes() ):
        registre_successors[node]=list( G.successors(node) )

    return registre_successors

def check_path(G):
    """
    This finction chech that if a graph G is a road.
    G is a road iff it is connected graph and a path from 'start_node' to 'end_node' is exactly G
    @param G : a graph
    @pram return boolean 
    """
    node_start=[node for node,d in  G.in_degree() if d==0] 
    node_end=[node for node,d in G.out_degree() if d==0]
    if len(node_start ) != 1 :
        return False
    elif len(node_end) !=1: 
        return False

    else:
        nodes_without_start_end= [node for node in G.nodes if node != node_start[0] and node != node_end[0]]
        for node in nodes_without_start_end:
            if G.in_degree(node) !=1  or G.out_degree(node) != 1:
                return False
                
        return True      


def recontruct_path(G):
    """
    This function take a graph G (which has to be a path) and return an ordered list of nodes of  G .

    """
    if check_path(G):
        list_nodes=[]
        node_start=[node for node,d in  G.in_degree() if d==0] 
        if node_start !=[]: # G is not empty
            while node_start != []:
                list_nodes.append(node_start[0] )
                node_start=list( G.successors(node_start[0]) ) # a list


        return list_nodes
    else:
        print("G isn't a path")
         

