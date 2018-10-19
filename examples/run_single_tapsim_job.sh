#!/bin/bash

##runs a single TAPSIM simulation given an existent TPASIM input meshgen results single-threaded
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de

simid=$1

##make folder to isolate results
output_prefix="BatchMillSimID${simid}"

##currently I am sitting in top folder
##generate therein above-mentioned simulation specific folder and go into it
workingdir="$output_prefix"
eval "mkdir -p $workingdir"
eval "cd $workingdir"

##shared utilization of cfg config and meshgen ini file as pillars sample different orientations of same crystal structure only

##copy generated input mesh over to working folder
tapsim_i="TAPSIMMesh.SimID.${simid}"
eval "cp ../$tapsim_i ../$workingdir"

##start TAPSIM simulation
tapsim_o="TAPSIMResult.SimID.${simid}"

eval "./../tapsim --threads=1 evaporation ../AlMeshCorrected.cfg $tapsim_i $tapsim_o 1>TAPSIM.SimID.${simid}.STDOUT.txt 2>TAPSIM.SimID${simid}.STDERR.txt"