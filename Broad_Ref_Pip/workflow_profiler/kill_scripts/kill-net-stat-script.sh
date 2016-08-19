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
ps -aef | grep '^'${me}'\|^'${myid} | grep 'netstat ' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $2}' > ./dum_PID_netstat
sort ./dum_PID_netstat > ./dum_PID_netstat.sorted
netstat_PID=`tail -1 ./dum_PID_netstat.sorted | awk '{printf "%8d", $1}'`
echo "netstat_PID: "${netstat_PID}
while read line; do    
  kill -s SIGUSR1 ${line}
done < dum_PID_netstat.sorted
# Now kill that process
rm dum_PID_netstat*
#    perl ./sar_parser.pl -f ${SARTSV} -o ${SAROUT}.csv &&

