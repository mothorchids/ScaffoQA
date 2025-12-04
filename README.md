# Scaffolding Quantum Approach (ScaffoQA)

This repository provides an implementation of an advanced approach to **Metagenomic de Novo Assembly** by reformulating the problem as a **Quadratic Unconstrained Binary Optimization (QUBO)** task. This formulation makes the assembly process compatible with both quantum algorithms, such as the **Quantum Approximate Optimization Algorithm (QAOA)**, and classical solvers, including Google’s OR-Tools.  The project contributes a specialized pipeline that adapts quantum optimization techniques to the complexities of real-world metagenomic data.

This work was developed during the [**CEMERACS**](https://cemracs2025.math.cnrs.fr/) summer school under the supervision of *Jérôme Gomar* (L'Oréal). The focus is on transforming the complex assembly graph structure into a quantum-ready format and delivering a complete pipeline, spanning the initial mathematical formulation through to preliminary evaluation on given data. To address scalability issues inherent to quantum computation, we also introduce a graph decomposition method designed to significantly reduce the size of the optimization problem before encoding it into a quantum system. This decomposition step improves the feasibility of applying QAOA and similar algorithms on current, near-term quantum devices.

## Documentation 

Location of the documentation: ```documentation```.

## Examples of usage

```
python graph_to_QUBO.py {path_to_data} -k {kmer}
```
Notes:
- kmer is an integer.

This commande will generate the QUBO matrix and save it in a ```.npy``` file.

```
python find_solution.py {path_QUBO_matrix.npy} -p 5 -s 1024
```
Notes:
- p is the number of layes in QAOA algorithm.
- s is the number of experiments.

## Structure
```
├── README.md
├── example.ipynb
├── decomposition_genome_k_121.ipynb
├── longest_path_genom.ipynb
├── documentation
│   ├── index.html
│   └── utility.html
├── src
│   ├── find_solution.py
│   ├── graph_path_problem.py
│   ├── graph_to_QUBO.py
│   ├── qubo_util.py
│   ├── reconstruit.py  
│   ├── tools.py
│   ├── tools_base.py
│   └── utility.py
├── input
│   └── spneumoniae_k_121.unitigs.fa
└──  output
   ├── graph_genom_longest_path.html
   └── graph_genome_de-Bruijn.html

```
