#!/bin/bash

export OMP_NUM_THREADS=2
export MY_OCTOPUS="/work/sato/octopus_tutorial/octopus_execute/install/bin/octopus"


cp inp_gs inp

mpirun -np 1  "${MY_OCTOPUS}" >log.gs.log 2>&1
