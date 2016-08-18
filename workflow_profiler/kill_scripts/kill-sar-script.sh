#!/bin/bash

#####################################
# Copyright (c) Intel Corporation
# PLEASE DO NOT SHARE WITHOUT NDA
#####################################

unset me
unset sar_PID

export me=`whoami`
export myid=`id -u`
#echo $SARDATA
#  Capture and record the PID of the correct sar 'counter' process:
#ps -aef | grep ${me} | grep 'sar ' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $2}' > ./dum_PID_sar
ps -aef | grep '^'${me}'\|^'${myid} | awk '{print $2,$8}' | grep 'sar' | grep -v 'grep' | grep -v 'kill' | awk '{printf "%8d", $1}' > ./dum_PID_sar

sort ./dum_PID_sar > ./dum_PID_sar.sorted
sar_PID_sar=`tail -1 ./dum_PID_sar.sorted | awk '{printf "%8d", $1}'`
echo "sar_PID: "${sar_PID}
while read line; do    
  #echo $line
  kill -s SIGUSR1 ${line}
done < dum_PID_sar.sorted

# Now kill that process
#kill -s SIGUSR1 ${sar_PID}
rm dum_PID_sar*
#sadf -p ${SARDATA} -- -A > ${SARTSV}  
#perl /opt/sar_parser.pl -f ${SARTSV} -o ${SAROUT}.csv 
