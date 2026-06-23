Author: Ezra Baup

_Master's thesis, Spring 2026, Computational Science and Engineering (MA, EPFL)._

This repository contains all code associated to the numerical experiments done during my Master's thesis. It is organized in two folders, corresponding to the two main parts of my thesis:
- `RJDBASE_opti`: refinement of solutions from the RJD-BASE framework through optimization
- `RJEA`: refinement of solutions from the randomized joint eigenvalues method

## RJDBASE_opti
Code files are:
- `main.ipynb`: notebook to run experiment
- `helpers_opti`: helper functions for optimization implementation
- `helpers_output`: helper functions to format and export output
- `helpers_setup`: helper functions for input setup 
- `rjd_base`: implementation of original RJD-BASE algorithm

The folder `caltech-101` contains data used for the Caltech example. The original data can be found at https://data.caltech.edu/records/mzrjq-6wc02, the preprocessed mat file was given by A. Pados [3].

To run the experiment, simply change the variable 'example' in the third cell to choose betwee, the synthetic graph 1, the weighted SBM (synthetic graph 2) or the real-world example from the Caltech dataset. Then run the notebook. The variable 'nb_repet' controls how many trials we do of each method, so that the results are then average over these trials. The results can be visualized in the second to last cell, and are exported to a LaTeX table. They represent the average of each measure with 95% confidence interval. The last cell also allows to save them as pickle files separately (mean, std, confidence interval).


## RJEA
The main scripts are:
- `main_rje`: main experiment
- `main_rootfinding`: experiment about rootfinding application
- `analysis_nb_trials`: analysis method vs number of trials
- `analysis_n_ex1`: analysis method vs n (size of matrices)
- `test_greedymatch`: auxiliary test to compare our greedy matchinf with Matlab's built-in (in the appendix of the thesis)

The other files are functions and auxiliaries. In particular, `config` does the setup of the experiments, `newRayleighOpt` redirects to the right solver depending on the configuration, and `new[method_name]` is the implementation of our RJEA-based algorithm _method_name_.
The folder `RJEA_He_Plestenjak` is a modified version of [1] for the original method.
The folder `rootfinding` contains auxiliaries & data for the rootfinding applications, see its README for  specific references.

To run the main experiments, one needs to determine the setting with a 'user' structure. The default setting is otherwise: example 1.1 with N=1000 runs and RQ1 version of the best residual (k=2) method. To modify the setting, run in the Matlab console:
```
user.field_name = var;
```
where field_name and var are summarized in the table below.

|   field_name  |    meaning    | values possible |
| ------------- | ------------- | --------------- |
| example  | family of input matrices  | between 0 and 3
| method  | algorithm used to solve the joint eigenvalues problem  | 'ref', 'best_res', 'mean', 'smaller_subspace'
| rq_type | Rayleigh quotient variant used by the solver | 'one-sided', 'two-sided', 'both' 
| nb_trials | number of trials of original method used by the new algorithms (k) | between 2 and size of the matrices, typically not too high
|Nsamples | number of runs done by the method (to obtain distribution of results) | integer value
| download_output | boolean that controls if we want to save the images and tables created by the script | true/false

About the 'example' field:
- 0.1 corresponds to Example 5.1 of [2]
- 1.1-1.3 corresponds to our Example 4.1 with n=10,50,100 respectively.
- 1.5 corresponds to our Example 4.2
- 2.1 and 2.2 correspond to Example 5.3 of [2], which is our Example 4.3
- (4.1 and 4.2 are past examples that were not analyzed further due to the small size of the matrices.)

Then you can run the script in the Matlab console as usual. For example for the main experiment:
```
main_rje
```
For this script, the console will print the state of the run, and output the results as 2 tables and 2 images.

## References
[1] He, H., Plestenjak, B. (2024), RandomJointEig (v1.0.0), doi: 10.5281/zenodo.15144637, https://github.com/borplestenjak/RandomJointEig

[2] H. He, D. Kressner, and B. Plestenjak. "Randomized methods for computing joint eigenvalues, with applications to multiparameter eigenvalue problems and root finding". Numerical Algorithms, 100(3):861-892, 2024.

[3] H. He, A. Pados, and D. Kressner. "RJD-BASE: Multi-modal spectral clustering via randomized joint diagonalization", 2025.
