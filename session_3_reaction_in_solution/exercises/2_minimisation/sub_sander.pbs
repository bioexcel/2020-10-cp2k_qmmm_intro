#!/bin/bash --login
#PBS -N sander_equil
# Select 1 node
#PBS -l select=1
#PBS -q short
#PBS -l walltime=00:20:00

# Replace this with your budget code
#PBS -A ta003

# Move to directory that script was submitted from
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR

module load amber-tools/20

# Run LAMMPS
echo "Sander minimisation"
aprun -n 24 sander.MPI -O \
              -i in.classical_minimisation -o out.classical_minimisation \
              -p system.prmtop -c system.inpcrd -r system.min.r
echo "Sander thermal equilibration"
aprun -n 24 sander.MPI -O \
              -i in.classical_heating -o out.classical_heating \
              -p system.prmtop -c system.min.r -r system.md.r \
              -x system.md.nc

