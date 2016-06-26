### Repository for *Extensively-parameterized mutation–selection models reliably capture site-specific selective constraint*, by SJS\* and COW.

### Contents of repository:

- `data/` contains all simulated alignments, and the 512-taxon balanced trees (with branch lengths of either 0.5 or 0.01) used during simulation. Alignments are named with this format: ``<dataset_name>_bl<0.5/0.01>.phy``, where ``bl`` indicates the branch lengths of the tree used for simulation..

- `simulation/` contains all code used for simulating sequences, as well as simulating parameters for use in sequence simulation.
    - [`ramsey2011_alignments`](./simulation/ramsey2011_alignments) contains all sequence alignments from [Ramsey et al. (2011)](http://www.genetics.org/cgi/pmidlookup?view=long&pmid=21467571).
    - [`derive_natural_simulation_parameters.py`](./simulation/derive_natural_simulation_parameters.py) derives parameters for simulating **natural** sequences
    - [`derive_dms_simulation_parameters.py`](./simulation/derive_dms_simulation_parameters.py) derives parameters for simulating **DMS** sequences. Note that experimental preferences are in the directory [`true_simulation_parameters`](./simulation/true_simulation_parameters)
    - [`true_simulation_parameters`](./simulation/true_simulation_parameters) contains all true parameters for simulation, including true *dN/dS* and entropy, amino-acid fitness, codon frequencies, and selection coefficients
    - [`simulate_alignments.py`](./simulation/simulate_alignments.py) simulates a sequence alignment, specifically on UT's (now defunct..) PhyloCluster.

- `inference/` contains all code used for mutation-selection model inference. All scripts named `*.sh` and `*.qsub` are used for submitting jobs to UT's Phylocluster, and all `*.py` scripts conduct and process inferences.

- `results/` contains all inference results.
    - [`swmutsel/`](./results/swmutsel/) contains all inference results with swMutSel for a variety of penalizations, indicated in file name. The script [./results/extract_sw_fitness.py](./results/extract_sw_fitness.py) extracts fitness values from the MLE inferences from swMutSel into separate text files for later use
    - [`phylobayes/`](./results/phylobayes/) contains all inference results with pbMutSel

- `postprocessing/` contains all code used to process, analyze, and plot data (mostly `R` scripts). All generated plots are also in this directory. **Note**: R code requires the packages cowplot, ggrepel, dplyr, tidyr, readr, grid, lme4, multcomp, and lmerTest.
    - [`calculate_inferred_quantities.py`](./postprocessing/calculate_inferred_quantities.py) calculates *dN/dS*, entropy, selection coefficient distributions, and JSD for inferences in `results/`. Resulting quantities are in the subdirectory [dataframes](./postprocessing/dataframes).
    - [`process_results.R`](./postprocessing/process_results.R) processes inference results in [dataframes](./postprocessing/dataframes) to create the final csv file [inference_results.csv](./postprocessing/dataframes/inference_results.csv)
    - [`plot_figures.R`](./postprocessing/plot_figures.R) makes all the figures in the manuscript. Figures are saved in either in the subdirectory [`maintext_figures/`](./postprocessing/maintext_figures/) or [`SI_figures/`](./postprocessing/SI_figures/)

- `universal_functions.py` is just a python module containing various functions used throughout the repository.


**\*This repository is maintained by SJS. Please file any questions/comments in [Issues](https://github.com/sjspielman/mutsel_benchmark/issues/), or contact me at stephanie.spielman@gmail.com.**
