#!/bin/bash

#SBATCH -J test                                      
#SBATCH -N 1                                       
#SBATCH -n 4                                       
#SBATCH -o output.dat                               
#SBATCH -e test.error
#SBATCH -p backfill2                          
#SBATCH --mem-per-cpu=MaxMemPerCPU                             
#SBATCH -t  00-04:00:00                              

export myscratch=/gpfs/research/scratch/nhv17/psi4.$SLURM_JOBID
rm -rf $myscratch
mkdir -p $myscratch
export PSI_SCRATCH=$myscratch
psi4 -i input.dat -o output.dat -n 4 
rm -rf $myscratch
