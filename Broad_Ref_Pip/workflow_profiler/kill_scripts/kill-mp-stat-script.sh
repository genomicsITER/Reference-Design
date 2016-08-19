#!/bin/bash

#####################################
# Copyright (c) Intel Corporation
# PLEASE DO NOT SHARE WITHOUT NDA
#####################################

unset me
unset mpstat

export me=`whoami`
export myid=`id -u`

#  Capture and record the PID of the correct emon 'counter' process:
ps -aef | grep '^'${me}'\|^'${myid} | awk '{print $2,$8}' | grep 'mpstat' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $1}' > ./dum_PID_mpstat
sort ./dum_PID_mpstat > ./dum_PID_mpstat.sorted
mpstat_PID_mpstat=`tail -1 ./dum_PID_mpstat.sorted | awk '{printf "%8d", $1}'`
echo "mpstat_PID_mpstat: "${mpstat_PID_mpstat}
# Now kill that process
while read line; do    
  #echo $line    
  kill -s SIGUSR1 ${line}
done < dum_PID_mpstat.sorted
#kill -s SIGUSR1 ${mpstat_PID_mpstat}
rm dum_PID_mpstat*
#    perl ./sar_parser.pl -f ${SARTSV} -o ${SAROUT}.csv &&

