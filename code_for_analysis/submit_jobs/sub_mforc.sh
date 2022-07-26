#!/bin/bash # shell type

#===========================================================================
# Execute this script: sh sub_mforc.sh
#===========================================================================

Dir="/nethome/yxl1496/hycom/GS_HR/filter"
DR="/nethome/yxl1496/hycom/GS_HR/filter/run"

cd $DR

for ((IDsh=40; IDsh<=106; IDsh=IDsh+1)); do
  #
  daysSh=$IDsh
  dayeSh=$(($daysSh + 1))

  # copy .f file and edit
  cat "$Dir/calc_mean_uvdp_HRA.f" >> calc_mforc\_$IDsh.f
  # replace time in .f
  sed -i "s#daysSh#$daysSh#g" calc_mforc\_$IDsh.f
  sed -i "s#dayeSh#$dayeSh#g" calc_mforc\_$IDsh.f

  # compile
  gfortran calc_mforc\_$IDsh.f -o mforc\_$IDsh
  
  # cp job-submitting script and edit
  cat "$Dir/run_meanforc_seq" >> run_mforc\_$IDsh
  sed -i "s#IDsh#$IDsh#g" run_mforc\_$IDsh

  # submit the script
  bsub < run_mforc\_$IDsh

done

