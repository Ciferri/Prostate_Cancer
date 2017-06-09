#!/bin/bash
rm img/*.png
octave printInputs.m
doubTime=1008
G1Dur=0.55
SDur=0.2
G2Dur=0.15
MDur=0.1
G1Distrib=0.6 
SDistrib=0.25
G2Distrib=0.075
MDistrib=0.075
alphaAlive=0
alphaTumG1=0.158
alphaTumS=0.113
alphaTumG2=0.169
alphaTumM=0.189
alphaDead=0
alphaVes=0
betaAlive=0
betaTumG1=0.051
betaTumS=0.037
betaTumG2=0.055
betaTumM=0.061
betaDead=0
betaVes=0
apopProb=0.8
apopDeadTime=234
necDeadTime=268
#for i in `seq 1 5`;
#do
./m2slv01 $doubTime $G1Dur $SDur $G2Dur $MDur $G1Distrib $SDistrib\
	  $G2Distrib $MDistrib $alphaAlive $alphaTumG1 $alphaTumS\
	  $alphaTumG2 $alphaTumM $alphaDead $alphaVes $betaAlive\
	  $betaTumG1 $betaTumS $betaTumG2 $betaTumM $betaDead\
	  $betaVes $apopProb $apopDeadTime $necDeadTime
#done
octave printOutImg.m
