#!/bin/bash --login
#PBS -N QMMM_META
# Select 1 node
#PBS -l select=1
#PBS -q short
#PBS -l walltime=00:20:00

# Replace this with your budget code
#PBS -A ta003

# Move to directory that script was submitted from
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR

module load cp2k/7.1

# Run LAMMPS
aprun -n 24 cp2k.popt -i inp.qmmm_link_atoms_metadynamics \
                      -o out.qmmm_link_atoms_metadynamics

