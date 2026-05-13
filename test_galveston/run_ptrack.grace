#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE        #Do not propagate environment
#SBATCH --get-user-env     #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=ptrack     #job name
#SBATCH --time=05:00:00            #time limit
#SBATCH --ntasks=1                #number of request cpus
#SBATCH --ntasks-per-node=1        #number of cpus per node
#SBATCH --mem=32G               #request memory per node
#SBATCH --output=Example1Out.%j    #Send stdout/err to "Example1Out.[jobID]"

module purge
module load intel/2023b
module load netCDF/4.9.2
module load netCDF-Fortran/4.6.1
module load CMake/3.12.1
export NETCDF=/sw/eb/sw/netCDF/4.9.2-iimpi-2023b
export NETCDF_FORTRAN=/sw/eb/sw/netCDF-Fortran/4.6.1-iimpi-2023b

rm Example1Out.*
./ptrack3.WW > mirror.out

