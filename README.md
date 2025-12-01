# Scaffolding Quantum Approach (ScaffoQA)

Source code for CEMERACS2025 research project is in ```src/```.

## Examples 
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

## Documentation 
Location of the documentation: ```documentation```.