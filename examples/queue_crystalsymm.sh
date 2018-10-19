#!/bin/bash

##runs a number of BiCarving jobs couple to issuing console to generate meshgen pillar configurations
##Markus Kuehbach, 2018/10/01, m.kuehbach at mpie.de
exenm="bicarving_crystalsymm"

simid="1"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +0.0 +0.0 +0.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"

simid="2"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +0.0 +0.0 +90.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="3"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +0.0 +0.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="4"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +0.0 +0.0 +270.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="5"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +1.0 +0.0 +90.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="6"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +1.0 +0.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="7"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +1.0 +0.0 +270.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="8"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +0.0 +1.0 +90.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="9"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +0.0 +1.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="10"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +0.0 +1.0 +270.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"

simid="11"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +1.0 +0.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="12"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +0.0 +1.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="13"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +1.0 +1.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="14"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 -1.0 +0.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="15"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +0.0 -1.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="16"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +0.0 +1.0 -1.0 +180.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"

simid="17"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +1.0 +1.0 +120.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="18"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +1.0 +1.0 +240.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="19"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml -1.0 +1.0 +1.0 +120.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="20"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml -1.0 +1.0 +1.0 +240.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="21"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 -1.0 +1.0 +120.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="22"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 -1.0 +1.0 +240.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="23"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +1.0 -1.0 +120.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"
simid="24"
eval "mpirun -np 1 ./${exenm} ${simid} BiCarving.xml +1.0 +1.0 -1.0 +240.0 1>BiCarving.SimID.${simid}.STDOUT.txt 2>BiCarving.SimID.${simid}.STDERR.txt"