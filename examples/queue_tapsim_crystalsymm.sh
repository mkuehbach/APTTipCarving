#!/bin/bash

##Script to batch process multiple TAPSIM simulations given existent TAPSIM input meshes, one thread per pillar
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de

##user required input
#start id and end id
pids=$1
pidi=$2
pide=$3

for ((pid = $pids; pid <= $pide; pid += $pidi)); do
	sh run_single_tapsim_job.sh $pid &
done
