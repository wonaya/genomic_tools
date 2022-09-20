#!/bin/bash
#SBATCH -J tophat_BR1P      # Job Name
#SBATCH -o tophat_BR1P.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 1
#SBATCH -n 1     # Total number of mpi tasks requested
#SBATCH -p flat-quadrant        # Queue (partition) name -- normal, development, etc.
#SBATCH -t 48:00:00      # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A Thompson_Replication      # Project_name

source ~/.bashrc
export PATH=/work2/02114/wonaya/software/tophat-2.1.1.Linux_x86_64/:/work2/02114/wonaya/software/STAR/bin/Linux_x86_64/:$PATH

module load python3
module load launcher
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=tophat2.txt
export LAUNCHER_SCHED=interleaved
