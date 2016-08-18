#!/bin/bash

#####################################
# Copyright (c) Intel Corporation
# PLEASE DO NOT SHARE WITHOUT NDA
#####################################

unset me
unset collectl_PID

export me=`whoami`
export myid=`id -u`

#  Capture and record the PID of the correct emon 'counter' process:
#ps -aef | grep ${me} | grep 'collectl' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $2}' > ./dum_PID_collectl
ps -aef | grep '^'${me}'\|^'${myid} | awk '{print $2,$10}' | grep 'collectl' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $1}' > ./dum_PID_collectl
sort ./dum_PID_collectl > ./dum_PID_collectl.sorted
collectl_PID=`tail -1 ./dum_PID_collectl.sorted | awk '{printf "%8d", $1}'`
echo "collectl_PID: "${collectl_PID}
while read line; do    
  #echo $line
  kill -s SIGUSR1 ${line}
  kill -s SIGTERM ${line}
done < dum_PID_collectl.sorted

# Now kill that process
#if [[ -z "${collectl_PID}" ]]; then 
#  kill -s SIGUSR1 ${collectl_PID}
#fi
rm dum_PID_collectl*
#    perl ./sar_parser.pl -f ${SARTSV} -o ${SAROUT}.csv &&

