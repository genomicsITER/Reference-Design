#!/bin/bash

#####################################
# Copyright (c) Intel Corporation
# PLEASE DO NOT SHARE WITHOUT NDA
#####################################

unset me
unset iostat_PID

export me=`whoami`
export myid=`id -u`

#  Capture and record the PID of the correct emon 'counter' process:
#ps -aef | grep ${me} | grep 'iostat ' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $2}' > ./dum_PID_iostat
ps -aef | grep '^'${me}'\|^'${myid} | awk '{print $2,$8}' | grep 'iostat' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $1}' > ./dum_PID_iostat
sort ./dum_PID_iostat > ./dum_PID_iostat.sorted
iostat_PID=`tail -1 ./dum_PID_iostat.sorted | awk '{printf "%8d", $1}'`
echo "iostat_PID: "${iostat_PID}
while read line; do    
  #echo $line
  kill -s SIGUSR1 ${line}
done < dum_PID_iostat.sorted

# Now kill that process
#if [[ -z "${iostat_PID}" ]]; then 
#  kill -s SIGUSR1 ${iostat_PID}
#fi
rm dum_PID_iostat*
#    perl ./sar_parser.pl -f ${SARTSV} -o ${SAROUT}.csv &&

