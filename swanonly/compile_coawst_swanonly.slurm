#!/bin/bash
 
## Job Name
#SBATCH --job-name=compile-dunex

## Allocation Defintion 
#SBATCH --account=derakhti
#SBATCH --partition=derakhti-gpu

## Resources
## Total number of nodes
#SBATCH --nodes=1

## Number of cores per node
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1

## Walltime
#SBATCH --time=00:30:00

## Memory per node
#SBATCH --mem=5G

## Turn on email notification 
##SBATCH --mail-type=ALL
##SBATCH --mail-user=erainvil@uw.edu

## Export all your environment variables to the batch job session 
#SBATCH --export=all

## Add appropriate compilers
module purge
module load netcdf_fortran+c_4.4.1.1-icc_17
module load icc_17-ompi_2.0.2
module list

## Set Environment Variables
export NETCDF_INCDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/include
export NETCDF_LIBDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/lib
export NETCDF=/sw/netcdf-fortran+c-4.4.1.1_icc-17
export NETCDF_CONFIG=/sw/netcdf-fortran+c-4.4.1.1_icc-17/bin/nc-config
export MCT_INCDIR=/gscratch/derakhti/shared/MCT/include
export MCT_LIBDIR=/gscratch/derakhti/shared/MCT/lib

## Code to be executed
## Compile the coawst application
/gscratch/derakhti/users/ej/dunex_swan/RUNFILES/coawst.bash > ./build.log
echo "done"
exit
