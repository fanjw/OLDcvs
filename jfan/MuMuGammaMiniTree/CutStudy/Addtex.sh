#!/bin/sh

for sample in 'Run2011A16Jan2012v1' 'Run2011B16Jan2012v1'
do
  iline=0
  while [ $iline -lt 9 ]
  do
    iTxt=0
    TotalMuMuGamma=0
    TotalMuMuGammaEvent=0
    let iLine=$iline+1
    while [ $iTxt -lt 10 ] 
    do
      export line=`sed -n ${iLine}p mumugammaCutsTotal_${sample}_${iTxt}.txt`
      export MuMuGamma=`echo $line | awk -F " " '{print $3}'`
      echo $MuMuGamma
      export MuMuGammaEvent=`echo $line | awk -F " " '{print $6}'`
      echo $MuMuGammaEvent
      let TotalMuMuGamma=$TotalMuMuGamma+$MuMuGamma
      let TotalMuMuGammaEvent=$TotalMuMuGammaEvent+$MuMuGammaEvent
      echo ${TotalMuMuGamma}
      echo ${TotalMuMuGammaEvent}
      let iTxt=$iTxt+1
   done
 
      echo "TOTALnbMuMuGammaAfterID["$iline"] = "$TotalMuMuGamma"                  TOTALnbEventsAfterMuMuGammaID["$iline"] = "$TotalMuMuGammaEvent >> mumugammaCutsTotal_${sample}.txt
      let iline=$iline+1 

 done
done
