"""
Reads a QUBO matrix and computes a solution using the QAOA algorithm 
from the myQLM library.

The script takes as input a QUBO matrix representing the quadratic 
unconstrained binary optimization problem and returns the solution 
provided by the QAOA quantum algorithm.

Dependencies:
- myQLM 

Usage:
- Provide the QUBO matrix as input.
- The script invokes QAOA to find an approximate solution.
- Outputs the solution vector corresponding to the minimum QUBO energy.
"""

import numpy as np
import argparse
#import re, os
# from qat.plugins import ScipyMinimizePlugin
from qat.opt import QUBO
from qat.qpus import get_default_qpu
from qat.plugins import ScipyMinimizePlugin
# from qat.core import Job
# from itertools import product
# import ast

# Parameters
#p_max = 9  # Set to desired maximum p
#grid_res = 10  # Coarse grid resolution for p = 1
#qpu = get_default_qpu()
# Problem setup
def main():
    parser = argparse.ArgumentParser(description="Find solution with QAOA")
    parser.add_argument("input_path", help="Path to the input file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("-p","--piter", type=int, required=True, help="number of iteration in QAOA")
    parser.add_argument("-s","--nbshot", type=int, required=True, help="number of shots")
    args = parser.parse_args()
    path_mat = args.input_path
    p = args.piter
    nbshot = args.nbshot

    print("Load Q")
    q_mat = np.load(path_mat)
    problem_QUBO = QUBO(Q=q_mat, offset_q=0)
    print("Problem QUBO formulated")
    observable = problem_QUBO.get_observable("terms")
    print(observable)
    print("=======================================")

    qpu = get_default_qpu()
    stack = ScipyMinimizePlugin(method="COBYLA",
                                tol=1e-5, 
                                options={"maxiter": 200}) | qpu
    # We can directly call the to_job method of the Problem class to pack an Ansatz and 
    # the cost observable in a single abstract Job
    job = problem_QUBO.to_job("qaoa", p) # Here 3 is the depth of the Ansatz
    result = stack.submit(job)
    print("Final energy:", result.value)
    print("=======================================")

    #Binding the variables:
    sol_job = job(**eval(result.meta_data["parameter_map"]))

    #Rerunning in 'SAMPLE' mode to get the most probable states:
    sampling_job = sol_job.circuit.to_job()

    sol_res = qpu.submit(sampling_job)
    print("Most probable states | as bitstrings | with probability ")
    for sample in sol_res:
        if sample.probability > 0.05:
            print(f"{sample.state} | {sample.state.value[0]} | {100 * sample.probability:.2f}%")

    print("=======================================")        
    # Emulating a reasonnable setup:
    sampling_job = sol_job.circuit.to_job(nbshots=nbshot)
    sol_res = qpu.submit(sampling_job)
    # Picking the most probable cut
    best_path = max([(s.state.value[0], s.probability) for s in sol_res], key=lambda s: s[1])[0]
    print(f"Best_path: {best_path}")


if __name__ == "__main__":
    main()