#!/bin/bash

export OMP_NUM_THREADS=4
export MY_OCTOPUS="/work/sato/octopus_tutorial/octopus_execute/install/bin/octopus"

cp inp_td inp

mpirun -np 4 "${MY_OCTOPUS}" >log.td.log 2>&1
