#!/bin/bash

#SBATCH --job-name=sp4.4x20.test	
#SBATCH --account=scw1019
#SBATCH --ntasks=25
#SBATCH --time=1-00:00

 
###
module load mpi
module load compiler/gnu/10 mpi/intel/2019/5

mpirun ./suN -i input_file -o output_file
