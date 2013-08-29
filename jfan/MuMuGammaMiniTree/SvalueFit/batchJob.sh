#! /usr/local/bin/bash -l
#$ -l ct=40000
#$ -P P_cmsf
#$ -l vmem=4G
#$ -l fsize=29G
#$ -q long
#$ -l sps=1
#$ -l dcache=1
###$ -l hpss=1
#$ -N SvalueFit_FAN
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
#$ -t 1


# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
SPSDIR="/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree/FitSvalue"
cd $SPSDIR
source ${VO_CMS_SW_DIR}/cmsset_default.sh
eval `scramv1 runtime -sh`



# ADD CURRENT DIRECTORY AND LIB TO LIRARY PATH
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}:/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/src"`




# COPY EXECUTABLE TO WORKER
cp *.C ${TMPDIR}
cp *.h ${TMPDIR}
cp *.txt ${TMPDIR}




#do jobs

#EndCaps=1
#r9sup=2
#Category=\"OneBin\"


cd ${TMPDIR}
#for EndCaps in 1;do
#    for r9sup in 2;do
#      for Category in OneBin;do


         #root -b -q  'SvalueFit.C ( '${EndCaps}', '${r9sup}', '${Category}' )'
         root -b -q  'SvalueFit.C ( '${1}', '${2}', '\"${3}\"' )'

#      done
#    done
#done




