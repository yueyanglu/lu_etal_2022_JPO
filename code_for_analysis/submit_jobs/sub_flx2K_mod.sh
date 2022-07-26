#!/bin/bash # shell type

#===========================================================================
# Execute this script: sh sub_flx2K_mod.sh
#===========================================================================

Directory="/nethome/yxl1496/HYCOM"

# calc K-tensor ('0') or iso K (scalar, '1') or K_red ('2')?
isoKSh=2

# if use div comp of the flx?
divSh=1

# if smooth the flx <F>?   F / <F> = -K * del_C
extraSh=0

# all numbers of tracers used to (over)-determine K
ncShs=(1 2 3 4 5 6 7 8 9 10)

# numbers of combinations corresponding to current # of tracers (ncShs)
ncombs=(10 45 120 210 252 210 120 45 10 1)

# which time interval 9
for lpSh in {1..1}; do

# which layer
for ((klaySh=24; klaySh<=24; klaySh=klaySh+1)); do

# loop over each # of comb 
for inc in {4..4}; do

ncSh=${ncShs[inc]}
ncomb=${ncombs[inc]}

# which flx term, 0 alleddy, 1 u'<c>, 2 <u>c', 3 u'c'
for tmSh in {3..3}; do

  # choose which tracer combination to determine K
  for ((icSh=77; icSh <= 77; icSh++)); do # $ncomb
  
    # which filter scale
    for ((smSh=4; smSh<=4; smSh=smSh+1)); do

      # filename to be created
      flnm=K$isoKSh\_L$lpSh\_Z$klaySh\_tm$tmSh\_P$icSh\of$ncSh\_S$smSh\_Div$divSh\_EX$extraSh.m
      # stripped of the last 2 characters, '.m'
      jbnm="${flnm%??}" 

      # copy codes into temporary m file
      cat "$Directory/calc_K_HYCOM_mod.m" >> $flnm

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

      # replace divSh in mfile with divSh
      sed -i "s#divSh#$divSh#g" $flnm

      # replace isoKSh in mfile with isoKSh
      sed -i "s#isoKSh#$isoKSh#g" $flnm

      # submit job  -o pt.o%J -e pt.e%J 
      bsub -J $jbnm -P ome -o pt.o%J -e pt.e%J -W 0:08 -q general -n 1 -R "rusage[mem=6000]" matlab -r $jbnm
  	#bsub -J $jbnm -P ome -o pt.o%J -e pt.e%J -W 8:08 -q parallel -n 16 -R "rusage[mem=6000] span[ptile=16]" matlab -r $jbnm

    done # filter scale
  done # tracer combinations
done # terms

done # inc for ncShs
done # layer
done # time interval
