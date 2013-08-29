#!/bin/sh


for EndCaps in 1 0; do
    for r9sup in 2 1 0; do
       for Category in OneBin; do

          qsub batchJob.sh ${EndCaps} ${r9sup} ${Category}
      
       done

    done
done






























