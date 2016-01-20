##### Repository for *mutsel model benchmark paper of TBD name*, by SJS\* and COW. 

##### Contents of respository:

- `data/` contains all simulated alignments, and the 512-taxon balanced tree used during simulation.

- `results/` contains all inference results.
    - [`swmutsel/`](./results/swmutsel/) contains all inference results with swMutSel for a variety of penalizations, indicated in file name.
    - [`phylobayes/`](./results/phylobayes/) contains all inference results with PhyloBayes
    - [`dnds_coeffs_jsd/`](./results/dnds_coeffs_jsd/) contains all *dN/dS*, selection coefficients, and JSD (Jensen-Shannon distance) values calculated from mutation-selection model inferences.
    - `summarized_results/` contains summary csv files for data in [`dnds_coeffs_jsd/`](./results/dnds_coeffs_jsd/). These files were created with [`process_results.R`](../postprocessing/process_results.R).

- `simulation` contains all code used for simulating sequences, as well as simulating parameters for use in sequence simulation. 
    - [`ramsey2011_yeast_alignments`](./simulation/ramsey2011_yeast_alignments) contains all sequence alignments from [Ramsey et al. (2011)](http://www.genetics.org/cgi/pmidlookup?view=long&pmid=21467571). 
    - [`derive_simulation_parameters.py`](./simulation/derive_simulation_parameters.py) derives parameters for sequence simulation, using alignments in [`ramsey2011_yeast_alignments`](./simulation/ramsey2011_yeast_alignments). Resulting parameters (including true *dN/dS*, amino-acid fitness, codon frequencies, and selection coefficients) are saved in [`true_simulation_parameters`](./simulation/true_simulation_parameters).
    - [`simulate_alignments.py`](./simulation/simulate_alignments.py) simulates a sequence alignment, specifically on UT's PhyloCluster. This script is run in conjuction with submission scripts `simulate_alignments.qsub` and `submit_simulations.sh`.
    
- `inference` contains all code used for mutation-selection model inference. All scripts named `*.sh` and `*.qsub` are used for submitting jobs to UT's Phylocluster, and all `*.py` scripts conduct and process inferences. 

- `postprocessing` contains all code used to process, analyze, and plot data (mostly `R` scripts). All generated plots are also in this directory.

- `universal_functions.py` is simply a python module containing various functions used throughout the repository.


\*This repository is maintained by SJS. Please file any questions/comments in [Issues](https://github.com/sjspielman/mutsel_benchmark/issues/).
