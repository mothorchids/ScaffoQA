import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pulp
from ortools.sat.python import cp_model
import matplotlib.colors as mcolors
from utility import *
from reconstruit import *
from tools_base import *

def bridge_graph(G):
    """
    This function determines the node of one incoming edge and one outgoing edge
    @parm G: a directe graphe

    @return a liste of node of degree 2: one incoming and one outgoing
    """
    node_deg_2=[]
    for node in list(G.nodes()):
        # One incoming  node and one outgoing node
        if G.in_degree[node]==1 & G.out_degree[node]==1:
            node_deg_2.append(node)

    return node_deg_2

def decomposition(G):
    """
    Function decompose into the subgraphs so that the nodes between two subgraphs are bridges
    Recall a node is a bridge is it has exactly one entrance edge and one exit edge

    @param G: a graph which is expected to be decomposed 

    @return a list of subgraphs with additional attribute on nodes
     The attributes indicate whether or node is a start or end of a subgraph

    """
    G_copy=G.copy()
    node_deg_2= bridge_graph(G)
    
    node_pre= []
    node_suc= []

    for node in node_deg_2:
        node_pre+= list( G.predecessors(node) ) # Always exist
        node_suc+= list( G.successors(node) ) 
        
        
    G_copy.remove_nodes_from(node_deg_2) # remove all bridge node

    list_connected_compossant = sorted(
        nx.weakly_connected_components(G_copy), key=len, reverse=True
    )  # Generate a sorted list of weakly connected components, largest first.
    
    #Creation of subgraph
    
    list_subgraph=[]
    k=0
    for set_nodes in list_connected_compossant:
        sub_graph= nx.subgraph(G, set_nodes).copy()
            
        for k,node in enumerate(node_deg_2):
            if sub_graph.has_node(node_pre[k]):
                #set_nodes.append(node)
                sub_graph.add_node(node)
                sub_graph.add_edge(node_pre[k], node)
                sub_graph.nodes[node]["nature"]="end_node"
                    
            

            if sub_graph.has_node(node_suc[k]) :
                #set_nodes.append(node)
                sub_graph.add_node(node)
                sub_graph.add_edge(node,node_suc[k])
                sub_graph.nodes[node]["nature"]="start_node"


            
            
        
        list_subgraph.append( sub_graph)        
    
    
    return  list_subgraph


def sharing_nodes(node_bridge, list_subgraph):
    """
    This function determines the subgraph of G sharing a specific  node

    @param G a grpah
    @return a registre where each term contains a liste of sharing subgraph of a specific node
    registre[node]=list
    
    
    """

    registre_bridge={}
    for node in node_bridge:
        list_sharing_sb_graphe=[]
        
        for subgraph in list_subgraph:
            if subgraph.has_node(node):
                list_sharing_sb_graphe.append(subgraph)
                
                
        registre_bridge[node]= list_sharing_sb_graphe 
                
            
    return  registre_bridge


def sharing_node_subgraph(node_bridge,list_subgraph):
    """
    Fonction which determines all bridges shared by two subgraph
    @param G : a graph
    @param return registre_bridge_int : is a registre containing a liste of bridge node
    @param return registres_subgraphe_int : is a registre containing a list of subgraph
    registre_bridge_int[k] contains all bridge node shared by the two subgraphs registres_subgraphe_int[k]
    """
    
    #node_bridge=bridge_graph(G) # list
    
    #liste_subgraph= decomposition(G) #list
    
    registre_bridge=sharing_nodes(node_bridge,list_subgraph) #registre

    node_bridge_copy=node_bridge.copy()
    registre_bridge_int={}
    registres_subgraphe_int={}
    k=0

    while node_bridge_copy !=[]:
        node=node_bridge_copy[0]

        liste_subgraphs= registre_bridge[node] # list of subgraph sharing node
        registres_subgraphe_int[k]= registre_bridge[node]
        registre_bridge_int[k]=[]
        (registre_bridge_int[k]).append(node)
        
        node_bridge_copy.remove(node)

        # Looking for other node sharing liste_subgraphs[0] and liste_subgraphs[1]
        for node in node_bridge_copy:
            if (liste_subgraphs[0]).has_node(node) & (liste_subgraphs[1]).has_node(node):
                (registre_bridge_int[k]).append(node)
                node_bridge_copy.remove(node)
                
        k+=1    
        

    return registre_bridge_int, registres_subgraphe_int


def sharing_subgraph(node_bridge,list_subgraph):
    registre_subgraph={}
    for subgraphe in list_subgraph:
        registre_subgraph[subgraphe]=[]
        
        for node in node_bridge:
            if subgraphe.has_node(node):
                (registre_subgraph[subgraphe]).append(node)


    return registre_subgraph


def opti_or_tool(G,start_node="",end_node="",alpha=20,beta=20,gamma=50,delta=10, limit_time=None):
    """
    This function uses the Google Or-Tools to determine the longest  path that visites each node at least once
    @param G is the graph
    @param start_node is the start node, 
    @param end_node is the end node
    @param alpha is the coefficient of the ongoing penalty 
    @param belta is the coefficient of the incoming penalty
    @ param gamma is the coefficient of the flux penalty
    @param delta is the coefficient of the total cost

    @return the longest path on the graph G and the value of the Hamiltonian if finded, otherwise, return None, None  
    
    """
    
    # ---------------------
    # Some useful tools
    # ---------------------
    edges_G=list( G.edges() )
    nodes_G=list( G.nodes() )
    
    registre_edges_G= registre_edges(G)
    registre_successors_G=registre_suc_node(G)
    registre_predecessors_G= registre_pred_node(G)

    n = len(edges_G)  # taille du problème

    # Paramètres

    #alpha = 20   # pondération des pénalités de ligne
    #beta = 20   # pondération des pénalités de colonne
    #gamma=10
    #delta=10

    model = cp_model.CpModel()

    # ---------------------
    # Variables binaires
    # ---------------------

    x = [model.NewBoolVar(f"x{i}") for i in range(n)]

    # ---------------------
    # Contraintes pour fixer certaines variables
    #----------------------
    if start_node !="": # start_node must be of degree 1
        successor_start_node=list( G.successors(start_node) )
        i_start=registre_edges_G[(start_node, successor_start_node[0] )]
        model.Add(x[i_start]==1) # x[edge_start] = 1
        

    if end_node !="": #end_node must be of degree 1
        predecessor_end_node=list(G.predecessors(end_node))
        i_end= registre_edges_G[(predecessor_end_node[0], end_node )]
        model.Add(x[i_end]==1)#x[edge_end]=1


    # ---------------------
    # Fonction objectif
    # ---------------------
    objective_terms = []

     #Terme lineaire
    for i in range(n):
            objective_terms.append(  -delta* x[i])

    for node in nodes_G:
    
        sucessor=registre_successors_G[node]
        for node_successor in sucessor:
             i=registre_edges_G[(node, node_successor )]
             objective_terms.append(  (gamma-alpha)* x[i])

        predecessor= registre_predecessors_G[node]
        for node_predecessor in  predecessor:
            j= registre_edges_G[(node_predecessor, node )]
            objective_terms.append(  (gamma-beta)* x[j])



    # Terme quadratique
    for node in nodes_G:
    
        sucessor=registre_successors_G[node]
        # Pénalités de lignes: 2(alpha+gamma)*x[i,j]*x[i,k] for j<k
        for node_successor_j in sucessor:
            for node_successor_k in sucessor:
                if node_successor_j != node_successor_k:
                
                    k= registre_edges_G[(node, node_successor_k )]
                    j=registre_edges_G[(node, node_successor_j )]
                    Y_kj= model.NewBoolVar(f"Y_{k}_{j}")
                
                    model.AddMultiplicationEquality(Y_kj, [x[k], x[j]])
                    objective_terms.append(  (alpha+gamma) * Y_kj)


        predecessor= registre_predecessors_G[node]
        # Pénalités de colonnes : 2*(beta+gamma)* x[i,j]*x[k,j] for i<k
        for node_predecessor_i in predecessor:
            for node_predecessor_k in predecessor:
                if node_predecessor_i != node_predecessor_k:
                
                    i=registre_edges_G[(node_predecessor_i,node)]
                    k=registre_edges_G[(node_predecessor_k,node)]
                    y_ik = model.NewBoolVar(f"y_{i}_{k}")
                
                
                    model.AddMultiplicationEquality(y_ik, [x[i], x[k]])
                    objective_terms.append(  (beta+gamma) * y_ik)

        # Pénalité de flux corrigée : -2*gamma* x[k,j]*x[i,k]
        for node_successor_k in sucessor:
            for node_predecessor_j in predecessor:
            
                k=registre_edges_G[(node, node_successor_k )]
                j=registre_edges_G[(node_predecessor_j,node)]
                z_jk = model.NewBoolVar(f"z_{j}_{k}")
            
            
                model.AddMultiplicationEquality(z_jk, [x[j], x[k]])
                objective_terms.append( -2* (gamma) * z_jk)



    model.Minimize(sum(objective_terms))

    # ---------------------
    # Résolution
    # ---------------------
    solver = cp_model.CpSolver()
    
    if limit_time !=None:
        solver.parameters.max_time_in_seconds = 10
        
    status = solver.Solve(model)

    # Récupération de la solution
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return [solver.Value(var) for var in x], solver.ObjectiveValue()
    else:
        return None, None
            
        
        