#!/bin/sh
#PBS -l nodes=2:WEST:ppn=8
#PBS -N phtrns_cc
                                                                                
NPROCS=`wc -l < $PBS_NODEFILE`

hostname
date


cd $PBS_O_WORKDIR

mpirun -np $NPROCS python /home/jhiskard/mylib/phonon2/exec_phtrns.py input.dat > stdout.txt
