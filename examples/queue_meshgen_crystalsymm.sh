#!/bin/bash

##runs a number of meshgen jobs uncoupled from the issuing console
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de

prefix_i="TAPSIMInputGeometry.SimID"
prefix_o="TAPSIMMesh.SimID"

simids="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24"

for simid in $simids; do
	eval "./meshgen $prefix_i.${simid}.dat $prefix_o.${simid} 1>$prefix_o.${simid}.STDOUT.txt 2>$prefix_o.${simid}.STDERR.txt &"
done
