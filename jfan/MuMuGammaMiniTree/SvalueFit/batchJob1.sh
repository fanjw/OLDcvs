

SPSDIR="/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree/FitSvalue"
cd $SPSDIR
source ${VO_CMS_SW_DIR}/cmsset_default.sh
eval `scramv1 runtime -sh`




#EndCaps=1
#r9sup=2
#Category=\"OneBin\"




for EndCaps in 1 0; do
    for r9sup in 2 1 0; do
      for Category in OneBin; do


         #root -b -q  'SvalueFit.C ( '${EndCaps}', '${r9sup}', '${Category}' )'
         root -b -q  'SvalueFit.C ( '${EndCaps}', '${r9sup}', '\"${Category}\"' )'
         #root -b -q  'SvalueFit.C ( '${1}', '${2}', '\"${3}\"' )'

      done
    done
done




