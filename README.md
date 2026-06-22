Author: Ezra Baup

_Master's thesis, Spring 2026, Computational Science and Engineering (MA, EPFL)._

This repository contains all code associated to the numerical experiments done during my Master's thesis. It is organized in two folders, corresponding to the two main parts of my thesis:
- `RJDBASE_opti`: refinement of solutions from the RJD-BASE framework through optimization
- `RJEA`: refinement of solutions from the randomized joint eigenvalues method

## RJDBASE_opti
The main scripts are:
- `main.ipynb`: notebook to run experiment
- `helpers_opti`: helper functions for optimization implementation
- `helpers_output`: helper functions to format and export output
- `helpers_setup`: helper functions for input setup 
- `rjd_base`: implementation of original RJD-BASE algorithm

## RJEA
- `main_rje`: main experiment
- `main_rootfinding`: experiment about rootfinding application
- `analysis_nb_trials`: analysis method vs number of trials
- `analysis_n_ex1`: analysis method vs n (size of matrices)
- `test_greedymatch`: auxiliary test to compare our greedy matchinf with Matlab's built-in (in the appendix of the thesis)

The other files are functions and auxiliaries. In particular, `config` does the setup of the experiments, `newRayleighOpt` redirects to the right solver depending on the configuration, and `new[method_name]` is the implementation of our RJEA-based algorithm _method_name_.
The folder `rootfinding` contains auxiliaries & data for the rootfinding applications, see its README for  specific references.

## References
He, H., Plestenjak, B. (2024), RandomJointEig (v1.0.0), doi: 10.5281/zenodo.15144637, https://github.com/borplestenjak/RandomJointEig

H. He, D. Kressner, and B. Plestenjak. "Randomized methods for computing joint eigenvalues, with applications to multiparameter eigenvalue problems and root finding". Numerical Algorithms, 100(3):861-892, 2024.

H. He, A. Pados, and D. Kressner. "RJD-BASE: Multi-modal spectral clustering via randomized joint diagonalization", 2025.

