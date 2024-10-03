# brains
Collection of scripts used to calculate fungi from ONT data derived from patient brain samples.

Scripts with prefixes "01_" to "03_" were designed to be used on internal HPC machines. All filepaths have been removed and will need to be updated by the user in order to work.

* "01_filter-reads" remove adapters from ONT reeads and performs a quality and length filtering
* "02_cluster" clusters the ONT reads using vsearch
* "03_assign-example" was ran in batches and enables characterisation of each cluster
* "calculateTaxa" uses bitscores to assign species or genus name to each cluster
* "visualiseBrains" generates graphics summarising the abundance profiles 


