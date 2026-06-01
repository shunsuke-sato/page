#!/bin/bash

export OMP_NUM_THREADS=4
export MY_OCTOPUS="/work/sato/octopus_tutorial/octopus_execute/install/bin/octopus"

cp inp_td_hhg inp

mpirun -np 4 "${MY_OCTOPUS}" >log.td_hhg.log 2>&1
