#!/bin/bash
 
## Job Name
#SBATCH --job-name=run-dunex

## Allocation Defintion 
#SBATCH --account=derakhti
#SBATCH --partition=derakhti-gpu

## Resources
## Total number of nodes
#SBATCH --nodes=1

## Number of cores per node
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks=28

## Walltime
#SBATCH --time=10:00:00

## Memory per node
#SBATCH --mem=100G

## Turn on email notification 
##SBATCH --mail-type=ALL
##SBATCH --mail-user=erainvil@uw.edu

## Export all your environment variables to the batch job session 
#SBATCH --export=all

## Add appropriate compilers
module purge
module load icc_17-ompi_2.0.2
## module load icc_17-impi_2017
module load netcdf_fortran+c_4.4.1.1-icc_17
module list

## Set Environment Variables- ifort
export NETCDF_INCDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/include
export NETCDF_LIBDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/lib
export NETCDF=/sw/netcdf-fortran+c-4.4.1.1_icc-17
export NETCDF_CONFIG=/sw/netcdf-fortran+c-4.4.1.1_icc-17/bin/nc-config
export MCT_INCDIR=/gscratch/derakhti/shared/MCT/include
export MCT_LIBDIR=/gscratch/derakhti/shared/MCT/lib

## Code to be executed
## run the coawstM file
mpirun -np 28 ./coawstM ./swan_dunex.in > coawst_dunex.out
