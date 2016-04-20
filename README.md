### Repository for *Extensively-parameterized mutation–selection models reliably capture site-specific selective constraint*, by SJS\* and COW. 

### Contents of repository:

- `data/` contains all simulated alignments, and the 512-taxon balanced trees (with branch lengths of either 0.5 or 0.01) used during simulation.
    - Alignments with parameters from yeast alignments are named with this format: ``<dataset_name>_del<regime>_bl<0.5/0.01>.phy``, where ``<regime>`` is either ``weak`` or ``strong,`` and ``bl`` indicates the branch lengths of the tree used for simulation.
    - Alignments with parameters from deep mutational scanning data are named with this format: ``<dataset_name>_bl<0.5/0.01>.phy``, where ``bl`` indicates the branch lengths of the tree used for simulation.

- `results/` contains all inference results.
    - [`swmutsel/`](./results/swmutsel/) contains all inference results with swMutSel for a variety of penalizations, indicated in file name.
    - [`phylobayes/`](./results/phylobayes/) contains all inference results with PhyloBayes
    - [`dnds_coeffs_jsd/`](./results/dnds_coeffs_jsd/) contains all *dN/dS*, selection coefficients, and JSD (Jensen-Shannon distance) values calculated from mutation-selection model inferences. These files were created with [`calculate_dnds_sc_jsd.py`](./postprocessing/calculate_dnds_sc_jsd.py).
    - [`summarized_results/`](./results/summarized_results) contains summary .csv files for data in [`dnds_coeffs_jsd/`](./results/dnds_coeffs_jsd/). These files were created with [`process_results.R`](./postprocessing/process_results.R).

- `simulation/` contains all code used for simulating sequences, as well as simulating parameters for use in sequence simulation. 
    - [`ramsey2011_yeast_alignments`](./simulation/ramsey2011_yeast_alignments) contains all sequence alignments from [Ramsey et al. (2011)](http://www.genetics.org/cgi/pmidlookup?view=long&pmid=21467571). 
    - [`derive_yeast_simulation_parameters.py`](./simulation/derive_yeast_simulation_parameters.py) derives parameters for sequence simulation, using alignments in [`ramsey2011_yeast_alignments`](./simulation/ramsey2011_yeast_alignments). Resulting simulation parameters (including true *dN/dS*, amino-acid fitness, codon frequencies, and selection coefficients) are saved in [`true_simulation_parameters`](./simulation/true_simulation_parameters).
    - [`derive_dms_simulation_parameters.py`](./simulation/derive_dms_simulation_parameters.py) derives parameters for sequence simulation, using data from deep mutational scanning experiments for influenza NP and HA. The input propensity parameters and resulting simulation parameters (including true *dN/dS* and codon frequencies) are saved in [`true_simulation_parameters`](./simulation/true_simulation_parameters). Here, fitness parameters and selection coefficients (in terms of their meaning as MutSel model paramters!) are not saved because they cannot be analytically calculated, due to the nature of this data. The raw preferences should be reasonably close to the true fitness parameters, and they are found in ``<HA/NP>_preferences.txt``.
    - [`simulate_alignments.py`](./simulation/simulate_alignments.py) simulates a sequence alignment, specifically on UT's PhyloCluster. This script is run in conjuction with submission scripts `simulate_alignments.qsub` and `submit_simulations.sh`.
    
- `inference/` contains all code used for mutation-selection model inference. All scripts named `*.sh` and `*.qsub` are used for submitting jobs to UT's Phylocluster, and all `*.py` scripts conduct and process inferences. 

- `postprocessing/` contains all code used to process, analyze, and plot data (mostly `R` scripts). All generated plots are also in this directory. **Note**: R code requires the packages cowplot, dplyr, tidyr, readr, lme4, multcomp, and lmerTest.
    - [`calculate_dnds_sc_jsd.py`](./postprocessing/calculate_dnds_sc_jsd.py) calculates *dN/dS*, selection coefficient distribitions, and JSD for inferences in `results/`.
    - [`process_results.R`](./postprocessing/process_results.R) processes inference results in [`dnds_coeffs_jsd/`](./results/dnds_coeffs_jsd/) to create csv files in [`summarized_results/`](./results/summarized_results).
    - [`plot_figures.R`](./postprocessing/plot_figures.R) makes all the figures in the manuscript. Figures are saved in the subdirectory [`figures/`](./postprocessing/figures/)
    - [`run_statistics.R`](./postprocessing/run_statistics.R) performs some statistical analyses on the results. All R output is commented out in this script itself.

- `universal_functions.py` is just a python module containing various functions used throughout the repository.


**\*This repository is maintained by SJS. Please file any questions/comments in [Issues](https://github.com/sjspielman/mutsel_benchmark/issues/), or contact me at stephanie.spielman@gmail.com.**
