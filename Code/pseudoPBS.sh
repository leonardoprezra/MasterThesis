#! /bin/bash
# PostProcessing of simulations

for PBS_ARRAYID in 400 500 600 700
do
  # run simulation
  python /nishome/students/leonardo/Dokumente/Thesis/Code/Gnan_Parameter_Histeresis3DSoft.py -n 32 -k $PBS_ARRAYID &
done

wait