#!/bin/bash

##runs a number of meshgen jobs uncoupled from the issuing console
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de

prefix_i="TAPSIMInputGeometry.SimID"
prefix_o="TAPSIMMesh.SimID"

simids="2004 2006 2008 2010 2012 2014 2016 2024"

for simid in $simids; do
	eval "./meshgen $prefix_i.${simid}.dat $prefix_o.${simid} 1>$prefix_o.${simid}.STDOUT.txt 2>$prefix_o.${simid}.STDERR.txt &"
done
