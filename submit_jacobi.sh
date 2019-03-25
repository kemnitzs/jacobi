#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=jacobi_example
#SBATCH --time=00:30:00
#SBATCH --nodes=1           ### Number of Nodes
#SBATCH --ntasks-per-node=1 ### Number of tasks (MPI processes)
#SBATCH --cpus-per-task=48  ### Number of threads per task (OMP threads)

OUTPUT_IMAGES=1

module purge
module load slurm cmake/3.4.3 gcc/8.2.0 openmpi/2.1.1
module list

echo $SLURM_CPUS_PER_TASK
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi

export OMP_NUM_THREADS=$omp_threads

BASE=$(pwd)
cmake .
make

rm -r run_

### setup run ###
RUNFOLDER=run_$SLURM_JOB_ID/
mkdir -p $RUNFOLDER 
cd $RUNFOLDER

cp $BASE/bin_jacobi .
# run script
./bin_jacobi

cd $BASE

### plot results ###
if [[ $OUTPUT_IMAGES -gt 0 ]]; then
  # for python3  
  module load tensorflow/1.8.0

  STARTTIME=$(date +%s)
  python3 plot_pcolor_out.py $RUNFOLDER
  ENDTIME=$(date +%s)
  echo "It takes $(($ENDTIME - $STARTTIME)) seconds to create all images..."
fi
