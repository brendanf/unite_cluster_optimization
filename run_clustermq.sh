#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name unite_clust_opt
#SBATCH --account project_2003104
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 64G
#SBATCH --time 36:00:00
#SBATCH --mail-type ALL

export PATH="/projappl/project_2003156/unite_cluster_optimization/bin:$PATH"
R --vanilla -f run_clustermq.R
