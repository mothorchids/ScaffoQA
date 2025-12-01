import numpy as np
# custom package
from graph_path_problem import GraphPathProblem
from qubo_util import get_qubo
from utility import *


def main():
    parser = argparse.ArgumentParser(description="QUBO Matrix Computing")
    parser.add_argument("input_path", help="Path to the input file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--info", action="store_true", help="Enable verbose output")
    parser.add_argument("-k","--kmer", type=int, required=True, help="K-mer size (default: 0)")
    args = parser.parse_args()
    path = args.input_path
    verbose = args.verbose
    k = args.kmer
    #filename = os.path.splitext(os.path.basename(path))[0]

    G = generate_de_bruijn_graph_from_fa(path,verbose)
    remove_isolated_nodes(G,verbose)
    DG = convert_de_bruijn_to_digraph(G,verbose)
    SG = find_biggest_component(DG)

    prob = GraphPathProblem(SG,k)
    prob.update_edges_register()

    Q = get_qubo(prob, delta=10, alpha=1, beta=1, gamma=5, start=None, finish=None)
    np.save(f"Q_graph_{k}.npy", Q)
    print(f"Output a QUBO matrix of dimension = {len(prob.edges_register)}")
  
if __name__ == "__main__":
    main()