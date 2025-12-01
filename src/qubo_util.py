import numpy as np
from graph_path_problem import GraphPathProblem

def symetrize(A):
    return A + A.T - np.diag(A.diagonal())

def get_qubo_main(prob:GraphPathProblem):
    n = len(prob.edges_register)
    Q = np.zeros((n, n))
    for i in prob.graph.nodes():
        for j in prob.graph.successors(i):
            diag_ind = prob.edges_register[(i,j)]
            Q[diag_ind , diag_ind] = prob.graph.edges[i, j].get("weight", 1)
    return Q

##########################################################
##########################################################
### Elementary terms
##########################################################

# quad terms
def add_quad_in(prob:GraphPathProblem,Q,i,val):
    """
    Computes the QUBO matrix for quadratique term of a input node i
    """
    for j in prob.graph.predecessors(i):
        for k in prob.graph.predecessors(i):
            ind_r = prob.edges_register[(j,i)]
            ind_c = prob.edges_register[(k,i)]
            Q[ind_r, ind_c] += val

def add_quad_out(prob:GraphPathProblem,Q,i,val):
    for j in prob.graph.successors(i):
        for k in prob.graph.successors(i):
            ind_r = prob.edges_register[(i,j)]
            ind_c = prob.edges_register[(i,k)]
            Q[ind_r, ind_c] += val

def add_quad_mix(prob:GraphPathProblem,Q,i,val):
    for j in prob.graph.successors(i):
        for k in prob.graph.predecessors(i):
            ind_r = prob.edges_register[(i,j)]
            ind_c = prob.edges_register[(k,i)]
            Q[ind_r, ind_c] += val

# lineaire termes 
def add_lin_in(prob:GraphPathProblem,Q,j,val):
    """
    Add lineaire incoming term from node j
    """
    for i in prob.graph.predecessors(j):
        diag_ind = prob.edges_register[(i,j)]
        Q[diag_ind , diag_ind] += val

def add_lin_out(prob:GraphPathProblem,Q,i,val):
    """
    Add lineaire outgoing term from node i
    """
    for j in prob.graph.successors(i):
        diag_ind = prob.edges_register[(i,j)]
        Q[diag_ind , diag_ind] += val

##########################################################
##########################################################
### Penalty Matrices
##########################################################

def get_qubo_out(prob:GraphPathProblem):
    n = len(prob.edges_register)
    Q = np.zeros((n, n))
    for i in prob.graph.nodes():
        add_lin_out(prob,Q,i,-2)
        add_quad_out(prob,Q,i,1)
    return Q

def get_qubo_in(prob:GraphPathProblem):
    n = len(prob.edges_register)
    Q = np.zeros((n, n))
    for j in prob.graph.nodes():
        add_lin_in(prob,Q,j,-2)
        add_quad_in(prob,Q,j,1)
    return Q

# flow
def get_qubo_flow(prob:GraphPathProblem):
    n = len(prob.edges_register)
    Q = np.zeros((n, n))
    for i in prob.graph.nodes():
        add_quad_out(prob,Q,i,1)
        add_quad_mix(prob,Q,i,-2)
        add_quad_in(prob,Q,i,1)
    return Q

def get_qubo_flow_start(prob:GraphPathProblem,start):
    n = len(prob.edges_register)
    Q = np.zeros((n, n))
    add_lin_in(prob,Q,start,-2)
    add_lin_out(prob,Q,start,2)
    return Q

def get_qubo_flow_finish(prob:GraphPathProblem,finish):
    n = len(prob.edges_register)
    Q = np.zeros((n, n))
    add_lin_in(prob,Q,finish,2)
    add_lin_out(prob,Q,finish,-2)
    return Q

##########################################################
##########################################################
### QUBO 
##########################################################

def get_qubo(prob:GraphPathProblem, delta=1, alpha=1, beta=1, gamma=1, start=None, finish=None):
    Q_main = get_qubo_main(prob)
    Q_out = get_qubo_out(prob)
    Q_in = get_qubo_in(prob)
    Q_flow = get_qubo_flow(prob)
    if start is not None:
        Q_flow += get_qubo_flow_start(prob,start)
    if finish is not None:
        Q_flow += get_qubo_flow_finish(prob,finish)
    Q = delta*Q_main + alpha*Q_out + beta*Q_in + gamma*Q_flow
    return symetrize(Q)