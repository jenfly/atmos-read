#!/bin/bash
#PBS -N python3
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=4:ppn=1
#PBS -l walltime=150:00:00
#PBS -q default
#PBS -j oe
#PBS -m bae -M jenfly@gmail.com

cd ~/dynamics/python/atmos-read
pwd

# Write out some information on the job
echo Running on host `hostname`
echo Time is `date`

### Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cpus

# Tell me which nodes it is run on
echo " "
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`
echo " "

# The command that runs the job.
python scripts/fram/run3.py

# ------
