#!/bin/bash

#####################################
# Copyright (c) Intel Corporation
# PLEASE DO NOT SHARE WITHOUT NDA
#####################################

unset me
unset free

export me=`whoami`
export myid=`id -u`

#  Capture and record the PID of the correct emon 'counter' process:
ps -aef | grep '^'${me}'\|^'${myid} | grep 'free -' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $2}' > ./dum_PID_free
sort ./dum_PID_free > ./dum_PID_free.sorted
free_PID=`tail -1 ./dum_PID_free.sorted | awk '{printf "%8d", $1}'`
echo "free_PID: "${free_PID}
while read line; do    
  kill -s SIGUSR1 ${line}
done < dum_PID_free.sorted

# Now kill that process
#if [[ -z "${free_PID}" ]]; then 
#  kill -s SIGUSR1 ${free_PID}
#fi
rm dum_PID_free*
#    perl ./sar_parser.pl -f ${SARTSV} -o ${SAROUT}.csv &&
