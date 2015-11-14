##### Repository for *mutsel model benchmark paper of TBD name*, by SJS\* and COW. 

##### Contents of respository:

- **data/** contains all alignments and phylogenies in subdirectories **simulated/** and **empirical/** for simulated and empirical data, respectively. 


- **results/** contains all inference results.
    - [**raw_results/**](./results/raw_results) contains all inference results, separated into subdirectories **simulated/** and **empirical/** for simulated and empirical data, respectively. Each of these subdirectories contains several subdirectories, as follows:
        - **swmutsel/** contains all inference results with swMutSel for a variety of penalizations, indicated in file name
        - **phylobayes/** contains all inference results with PhyloBayes
        - **slac/** contains *dN/dS* inferences performed with HyPhy using SLAC
        - **derived_dnds/** contains all *dN/dS* values as predicted from either swMutSel or PhyloBayes inference and, in the case of simulated data, the true *dN/dS* values.
        - For simulated data only, the subdirectory **jsd/** contains all per-site Jensen-Shannon distances between true and inferred amino-acid frequencies.
    - [**dnds_results.csv**](./results/dnds_results.csv) contains a summary of all inferred/true dN/dS across datasets. File created with [process_results.R](./scripts/postprocessing/process_results.R).
    - [**jsd_results.csv**](./results/jsd_results.csv) contains a summary of all JSD values for simulated datasets. File created with [process_results.R](./scripts/postprocessing/process_results.R).
    
- **scripts/** contains all code used in paper.
   - [calculate_dnds.py](./scripts/calculate_dnds.py) calculates *dN/dS* values from swMutSel and PhyloBayes inferences.
   - [calculate_jsd.py](./scripts/calculate_jsd.py) calculates Jensen-Shannon distances between simulated and inferred amino-acid frequencies, for simulated data.
   - [compute_dnds_from_mutsel.py](./scripts/compute_dnds_from_mutsel.py) contains module, published originally in [Pyvolve](https://github.com/sjspielman/pyvolve), for calculating *dN/dS* from MutSel model parameters.
   - [universal_functions.py](./scripts/universal_functions.py) contains various helper functions used throughout code.
   
   - [**empirical_data_collection/**](./scripts/empirical_data_collection) contains all code used to query PANDIT, process, and save empirical alignemnts.
   - [**postprocessing/**](./scripts/postprocessing) contains R code used to process, perform statistical analyses on, and plot results.
   - [**simulation/**](./scripts/simulation) contains all code used to simulate equilibrium frequencies and alignments. The subdirectory [flib/](./scripts/simulation/flib) contains the site-specific equilibrium frequencies used for simulation, as derived from structurally-curated yeast alignments from [Ramsey et. al. (2011)](http://www.genetics.org/content/188/2/479.full.pdf).
   - [**inference/**](./scripts/inference) contains code used for HyPhy *dN/dS* inference with SLAC, and swMutSel and PhyloBayes inference. The subdirectory **hyphy/batchfiles/** contains all the base HyPhy code necessary to run HyPhy. Specifically, these batchfiles have been edited by SJS to use F1x4 frequencies during inference. Any files that differ from those distributed by HyPhy are named with "_sjs".
   - All scripts ``*.sh`` and ``*.qsub`` were used for job submission on PhyloCluster.

\*This repository is maintained by SJS. Please file any questions/comments in [Issues](https://github.com/sjspielman/mutsel_benchmark/issues/).
