#!/bin/bash # shell type

#===========================================================================
# Execute this script: sh sub_tracflx_mod.sh
#===========================================================================

Directory="/nethome/yxl1496/HYCOM"

# which layer
for klay in {1..1}; do
 
  # which filter scale
  for ((smSh=4; smSh<=4; smSh=smSh+2)); do

    # filename to be created
    flnm=flx_S$smSh\of_Z$klay.m
    # stripped of the last 2 characters, '.m'
    jbnm="${flnm%??}" 

    # copy codes into temporary m file
    cat "$Directory/calc_tracflx_HYCOM_mod.m" >> $flnm

    # replace klaySh in mfile with klaySh
    sed -i "s#klaySh#$klay#g" $flnm

    # replace smSh in mfile with smSh
    sed -i "s#smSh#$smSh#g" $flnm 

    # submit job to cluster, -o pt.o%J
    bsub -J $jbnm -P ome -W 1:00 -q general -n 1 -R "rusage[mem=5000]" matlab -r $jbnm
  
  done
done 
