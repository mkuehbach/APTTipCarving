#!/bin/bash

##runs a number of meshgen jobs uncoupled from the issuing console
##Markus Kuehbach, 2018/10/16, m.kuehbach at mpie.de

prefix_i="TAPSIMInputGeometry.SimID"
prefix_o="TAPSIMMesh.SimID"

simids="6604 6606 6608 6610 6612 6614"

for simid in $simids; do
	eval "./meshgen $prefix_i.${simid}.dat $prefix_o.${simid} 1>$prefix_o.${simid}.STDOUT.txt 2>$prefix_o.${simid}.STDERR.txt &"
done
