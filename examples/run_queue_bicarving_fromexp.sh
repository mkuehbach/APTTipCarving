#!/bin/bash

##runs a number of BiCarving jobs couple to issuing console to generate meshgen pillar configurations
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de


##orientation matrix desired g11,g12,g13,g21, g22, ...
om="0.2833642 -0.70278351 0.65251599 0.28426857 0.711378811 0.64275199 -0.9159174 0.003346495 0.40136257"
exenm="bicarving_om_from_experiment"
xmlprefix="BiCarving_FromExperiment"

simid="2004"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_04nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="2006"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_06nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="2008"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_08nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="2010"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_10nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="2012"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_12nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="2014"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_14nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="2016"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_16nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="2024"
eval "mpirun -np 1 ./${exenm} ${simid} ${xmlprefix}_24nm.xml ${om} 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
