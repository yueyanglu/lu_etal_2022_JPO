#!/bin/bash # shell type

#===========================================================================
# Execute this script: sh sub_flx2KLambda_MOM.sh
#===========================================================================

Directory="$HOME/work_MOM"

# smooth the lhs?
extraSh=0

# fit the div or adv of flx? % 1 div; 2 adv
DivAdvSh=1

# 1 iso; 0 aniso
isoSh=0 

# numbers of tracers used to (over)-determine K
ncShs=(3 4 5 6 7 8 9 10)

# total numbers of combinations corresponding to # of tracers used (ncShs)
ncombs=(120 210 252 210 120 45 10 1)

# which time interval 
for lpSh in {1..10}; do

# which layer
for ((klaySh=24; klaySh<=24; klaySh=klaySh+1)); do

# loop over each of ncShs
for inc in {2..2}; do

ncSh=${ncShs[inc]}
ncomb=${ncombs[inc]}

# which flx term, 0 alleddy, 1 u'<c>, 2 <u>c', 3 u'c', 13 u'c
for tmSh in {13..13}; do
 
  # choose the tracer combinations to determine K
  for ((icSh=77; icSh <= 77; icSh++)); do  # 77  <= $ncombi, icSh++
  
    # which filter scale
    for ((smSh=4; smSh<=4; smSh=smSh+1)); do

      # filename to be created
      flnm=KL$DivAdvSh\_TM$tmSh\_L$lpSh\_Z$klaySh\_iso$isoSh\_P$icSh\of$ncSh\_S$smSh.m

      # stripped of the last 2 characters, '.m'
      jbnm="${flnm%??}" 

      # copy codes into temporary m file
      cat "$Directory/calc_KLambda_MOM.m" >> $flnm

      # replace lpSh in mfile with lpSh
      sed -i "s#lpSh#$lpSh#g" $flnm

      # replace klaySh in mfile with klaySh
      sed -i "s#klaySh#$klaySh#g" $flnm

      # replace icSh in mfile with icSh
      sed -i "s#icSh#$icSh#g" $flnm

      # replace smSh in mfile with smSh
      sed -i "s#smSh#$smSh#g" $flnm

      # replace extraSh in mfile with extraSh
      sed -i "s#extraSh#$extraSh#g" $flnm

      # replace ncSh in mfile with ncSh
      sed -i "s#ncSh#$ncSh#g" $flnm
      
      # replace tmSh in mfile with tmSh
      sed -i "s#tmSh#$tmSh#g" $flnm

      # replace DivAdvSh in mfile with DivAdvSh
      sed -i "s#DivAdvSh#$DivAdvSh#g" $flnm

      # replace isoSh in mfile with isoSh
      sed -i "s#isoSh#$isoSh#g" $flnm

      # submit job  -o pt.o%J -e pt.e%J 
      bsub -J $jbnm -P ome -W 0:05 -o pt.o%J -e pt.e%J -q general -n 1 -R "rusage[mem=6000]" matlab -r $jbnm
  
    done # filter scale
  done # tracer combinations
done # terms

done # inc for ncShs
done # layer
done # time interval
