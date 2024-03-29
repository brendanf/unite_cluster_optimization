#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name unite_cluster_optimization
#SBATCH --account project_2003156
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --nodes 1
#SBATCH --mem 96G
#SBATCH --time 36:00:00
#SBATCH --mail-type ALL
#SBATCH --gres=nvme:100

export OMP_STACKSIZE=8096
export OMP_THREAD_LIMIT=$SLURM_CPUS_ON_NODE
[ -z $OMP_THREAD_LIMIT ] && export OMP_THREAD_LIMIT=1
export PATH="/projappl/project_2003156/unite_cluster_optimization/bin:$PATH"
if [[ $1 == "test" ]] ; then
echo "testing outdated targets..."
R --vanilla --quiet -e 'targets::tar_outdated(callr_function=NULL)'
elif [[ $1 == "" ]] ; then
echo "building plan"
R --vanilla --quiet -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'
else
echo "building target '$1'"
R --vanilla --quiet -e "targets::tar_make($1, callr_function=NULL, reporter='timestamp')"
fi
