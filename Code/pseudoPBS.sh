#! /bin/bash
# PostProcessing of simulations

for PBS_ARRAYID in 400
do
  export HOOMD_WALLTIME_STOP=$((`date +%s` + 12 * 3600 - 10 * 60))
  mpirun hoomd run.py
  # run simulation
  python /nishome/students/leonardo/Dokumente/Thesis/Code/Gnan_NVE.py &
done

wait