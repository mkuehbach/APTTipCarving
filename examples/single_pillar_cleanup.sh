#!/bin/bash

##runs a single TAPSIM simulation given an existent TPASIM input meshgen results single-threaded
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de

simid=$1

##make folder to isolate results
output_prefix="BatchMillSimID${simid}"

##currently I am sitting in top folder go into working directory
workingdir="$output_prefix"
eval "cd $workingdir"

#delete dump files
eval "rm dump*"

#pack trajectory data then delete them
fn_o="TAPSIMMesh.SimID.${simid}.Trajectories.tar.gz"
eval "tar -czvf $fn_o trajectory_data*"
eval "rm trajectory_data*"

#pack supplementary geometry data then delete them
fn_o="TAPSIMMesh.SimID.${simid}.SupplGeometry.tar.gz"
eval "tar -czvf $fn_o grid_data* surface_data*"
eval "rm grid_data* surface_data*"

#pack geometry data then delete them
fn_o="TAPSIMMesh.SimID.${simid}.Geometry.tar.gz"
eval "tar -czvf $fn_o geometry.dat"
eval "rm geometry.dat"

#pack detector hits
fn_o="TAPSIMMesh.SimID.${simid}.DetectorHits.tar.gz"
eval "tar -czvf $fn_o TAPSIMResult.SimID.${simid}.*"

