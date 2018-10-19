#!/bin/bash

##runs a number of meshgen jobs uncoupled from the issuing console
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de

prefix_i="TAPSIMInputGeometry.SimID"
prefix_o="TAPSIMMesh.SimID"

simids="1005 1010 1015 1020 1025 1030 1035 1040 1045 1050 1055 1060 1065 1070 1075 1080 1085"

for simid in $simids; do
	eval "./meshgen $prefix_i.${simid}.dat $prefix_o.${simid} 1>$prefix_o.${simid}.STDOUT.txt 2>$prefix_o.${simid}.STDERR.txt &"
done
