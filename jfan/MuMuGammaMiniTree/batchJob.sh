#! /usr/local/bin/bash -l
#$ -l ct=40000
#$ -P P_cmsf
#$ -l vmem=4G
#$ -l fsize=29G
#$ -q long
#$ -l sps=1
#$ -l dcache=1
###$ -l hpss=1
#$ -N Selection_FAN
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
#$ -t 1-10

syntax="${0} {parameter}"
#if [[ -z ${6} ]]
if [[ -z ${2} ]]
then
	echo ${syntax}
	exit 1
fi


ijob=`echo "${SGE_TASK_ID} - 1" | bc -ql`
#ijob=9
echo ${ijob}
echo "USER=${USER}"



# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
#echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
###export HOMEDIR=/afs/in2p3.fr/home/o/obondu
###source ${HOMEDIR}/428v2.sh
#export HOMEDIR=/afs/in2p3.fr/home/s/sgandurr
#source ${HOMEDIR}/525p3.sh
export SPSDIR=/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree
cd $SPSDIR
source ${VO_CMS_SW_DIR}/cmsset_default.sh
#export SCRAM_ARCH=slc5_amd64_gcc434
eval `scramv1 runtime -sh`



# CHECK THE ENVIRONMENT VARIABLES
echo "CHECK THE ENVIRONMENT VARIABLES"
echo "ROOTSYS :" 
echo ${ROOTSYS}



# COPY HEADER FILES TO WORKER
#echo "COPY HEADER FILES TO WORKER"
#mkdir ${TMPDIR}/interface
#cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
#if [[ ! -e ${TMPDIR}/interface ]]
#then
#  mkdir ${TMPDIR}/interface
#	cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
#fi



# COPY IpnTree LIB FILE TO WORKER
#mkdir ${TMPDIR}/lib
#cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
#echo "COPY IpnTree LIB FILE TO WORKER"
#if [[ ! -e ${TMPDIR}/lib ]]
#then
#	mkdir ${TMPDIR}/lib
#	cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
#fi



# ADD CURRENT DIRECTORY AND LIB TO LIRARY PATH
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}:/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/src"`
#LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}:${WORKDIR}/lib:${WORKDIR}"`
#echo "LD_LIBRARY_PATH"
#echo ${LD_LIBRARY_PATH}
#echo ""




# COPY EXECUTABLE TO WORKER
echo "COPY EXECUTABLE TO WORKER"
###cp ${SPSDIR}/Zmumugamma/Selection/Selection_miniTree.exe ${TMPDIR}/
cp Selection_miniTree.exe ${TMPDIR}/
cp *.C ${TMPDIR}/
cp *.h ${TMPDIR}/
cp *.txt ${TMPDIR}/
##cp -r /sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2/ ${TMPDIR}/
cp list_* ${TMPDIR}/
##cp /sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/UserCode/IpnTreeProducer/ListZmumugamma/listFiles_* ${TMPDIR}/  ###GENERER LES LISTFILES ET CHANGER LE REPERTOIRE !!!




# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}
pwd
###./Selection_miniTree.exe ${1} ${2} ${3} ${4} ${5} ${6} 2> ${2}.err | tee ${2}.out
##./Selection_miniTree.exe ${1} ${2} 20 ${ijob} ${3} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out
##./Selection_miniTree.exe ${1} ${2} 300 ${ijob} ${3} 2011 PU_S6 ${4} ${5} ${6} ${7} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out
./Selection_miniTree.exe ${1} ${2} 10 ${ijob} ${3} 2011 PU_S6 ${4} ${5} ${6} ${7} ${8} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out ###CHANGER LE SCENARIO DE PILEUP !!! 




# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPDIR}/*.root ${SPSDIR}/OutputMiniTree
mv ${TMPDIR}/${2}_part${ijob}.out ${SPSDIR}/Logfile
mv ${TMPDIR}/${2}_part${ijob}.err ${SPSDIR}/Logfile
rm ${TMPDIR}/Selection_miniTree.exe
rm ${TMPDIR}/*.C
rm ${TMPDIR}/*.h
rm ${TMPDIR}/*.txt
##rm -rf ${TMPDIR}/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2/
rm ${TMPDIR}/list_* 

#"cd ${SPSDIR}/Zmumugamma/Selection/"
#cd ${SPSDIR}/Zmumugamma/Selection/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""



exit 0


