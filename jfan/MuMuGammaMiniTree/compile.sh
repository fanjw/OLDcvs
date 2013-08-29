#!/bin/bash
#syntax="Syntax ${0} {executableName}"
#if [[ -z ${1} ]]
#then
  #echo ${syntax}
                       #  exit 1   #not used
#fi
#exeFile=${1}


g++ Selection_miniTree.C -L${ROOTSYS}lib -lTMVA -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMLP -lTreePlayer -lMinuit -pthread -lm -ldl -rdynamic -pthread -m64 -I${ROOTSYS}include -L/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/src -lToto `root-config --libs --cflags` -o Selection_miniTree.exe



#g++ Selection_miniTree.C -L`pwd` -lToto `root-config --libs --cflags` -o ${exeFile}


