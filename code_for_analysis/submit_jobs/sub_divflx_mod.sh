#!/bin/bash # shell type

#===========================================================================
# Execute this script: sh sub_divflx_mod.sh
#===========================================================================

Directory="/nethome/yxl1496/HYCOM"

# which flx term, 0 alleddy, 1 u'<c>, 2 <u>c', 3 u'c'
tmSh=13

# if smooth the flx <F> before calc div comp
extraSh=0

# which layer
for ((klaySh=24; klaySh<=24; klaySh=klaySh+1)); do
#for klaySh in {1..30}; do

# which time interval  1-10 deeper than Z20; lp1-5 for Z02-Z19
for lpSh in {1..10}; do # 

  # no. of tracer
  for ((cSh=1; cSh<=9; cSh=cSh+2)); do

    # which filter scale
    for smSh in {4..4}; do
      # Make sure 2 or more parallel MATLAB jobs will not start at the same time
      # sleep 12s

      # filename to be created
      flnm=div_L$lpSh\_Z$klaySh\_tm$tmSh\_C$cSh\_S$smSh\_EX$extraSh.m

      # stripped of the last 2 characters, '.m'
      jbnm="${flnm%??}" 

      # copy codes into temporary m file
      cat "$Directory/calc_divflx_HYCOM_mod.m" >> $flnm

      # 1 replace lpSh in mfile with lpSh
      sed -i "s#lpSh#$lpSh#g" $flnm

      # 2 replace klaySh in mfile with klaySh
      sed -i "s#klaySh#$klaySh#g" $flnm

      # 3 replace tmSh in mfile with tmSh
      sed -i "s#tmSh#$tmSh#g" $flnm

      # 4 replace cSh in mfile with cSh
      sed -i "s#cSh#$cSh#g" $flnm

      # 5 replace smSh in mfile with smSh
      sed -i "s#smSh#$smSh#g" $flnm

      # 6 replace extraSh in mfile with extraSh
      sed -i "s#extraSh#$extraSh#g" $flnm

      # submit job  -o pt.o%J -e pt.e%J 
      bsub -J $jbnm -P ome -W 4:00 -q general -o pt.o%J -e pt.e%J -n 1 -R "rusage[mem=24000]" matlab -r $jbnm
      # bsub -J $jbnm -P ome -W 40:00 -q parallel -o pt.o%J -e pt.e%J -n 16 -R "rusage[mem=24000] span[ptile=16]" matlab -r $jbnm

    done # term
  done   # c
done     # lp
done     # z
