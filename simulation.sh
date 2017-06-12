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
tAlphaTumG1=(0.158 0.274 0.038 0.043 0.034)
tAlphaTumS=(0.113 0.197 0.027 0.031 0.024)
tAlphaTumG2=(0.169 0.293 0.041 0.046 0.036)
tAlphaTumM=(0.189 0.327 0.045 0.052 0.040)
alphaDead=0
alphaVes=0
betaAlive=0
tBetaTumG1=(0.051 0.028 0.025 0.031 0.024)
tBetaTumS=(0.037 0.02 0.018 0.022 0.017)
tBetaTumG2=(0.055 0.03 0.027 0.033 0.026)
tBetaTumM=(0.061 0.033 0.03 0.037 0.029)
betaDead=0
betaVes=0
apopProb=0.8
apopDeadTime=234
necDeadTime=268
Vmax=15.2
Kconso=3.035
for i in `seq 0 4`;
do
    #for j in `seq 1 5`;
    #do
    ./m2slv01 $doubTime $G1Dur $SDur $G2Dur $MDur $G1Distrib\
	      $SDistrib $G2Distrib $MDistrib $alphaAlive\
	      ${tAlphaTumG1[i]} ${tAlphaTumS[i]} ${tAlphaTumG2[i]}\
	      ${tAlphaTumM[i]} $alphaDead $alphaVes $betaAlive\
	      ${tBetaTumG1[i]} ${tBetaTumS[i]} ${tBetaTumG2[i]}\
	      ${tBetaTumM[i]} $betaDead $betaVes $apopProb\
	      $apopDeadTime $necDeadTime $Vmax $Kconso
    #done
done
octave printOutImg.m
