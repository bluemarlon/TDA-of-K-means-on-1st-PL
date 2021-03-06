#!/bin/bash
#SBATCH --partition=general
#SBATCH --ntasks=100
#SBATCH --time=06:00:00
#SBATCH --mail-type=END                       # Event(s) that triggers email notification (BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=renjie.chen@uconn.edu

module purge
module load java/1.8.0_162
module load zlib/1.2.11
module load gcc/5.4.0-alt
module load mpi


for i in `seq 2 5`;
do
   mpirun -np 100 ./Kmeans_MPI $i land tda > mpitda${i}time.txt
done
