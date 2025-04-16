#!/bin/bash

#cd /r04/atlas/amullin/ANUBIS/Simulations/MG5_aMC_v3_5_4

#/r04/atlas/amullin/ANUBIS/Simulations/MG5_aMC_v3_5_4/bin/mg5_aMC jobscript_hnl_Mathusla.txt 

#setupATLAS
#lsetup "root 6.30.02-x86_64-el9-gcc13-opt"

#rm -r MG5_aMC_v3_5_4 # remove existing directory (must make sure to have copied Events to SS dir) 
#mkdir MG5_aMC_v3_5_4_SS2
tar -xf MG5_aMC_v3_5_4.tar && mv MG5_aMC_v3_5_4 MG5_aMC_v3_5_4_SS1
cd MG5_aMC_v3_5_4_SS1

./bin/mg5_aMC jobscript_hnl_param_scan.txt
