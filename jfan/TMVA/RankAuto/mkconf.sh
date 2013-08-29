#!v/bin/sh
#variables xxx
varEBdt=( myphoton_s4ratio  event_rho )

varEEdt=( myphoton_pfchargedisogood03  myphoton_pfchargedisobad03  myphoton_pfphotoniso03  myphoton_sieie  myphoton_etawidth  myphoton_r9  myphoton_s4ratio  myphoton_SCeta  event_rho myphoton_ESEffSigmaRR ) 

varEBpf=( myphoton_pfchargedisogood03 myphoton_pfphotoniso03 myphoton_pfneutraliso03 myphoton_pfchargedisobad03 event_rho )

varEEpf=( myphoton_pfchargedisogood03 myphoton_pfphotoniso03 myphoton_pfneutraliso03 myphoton_pfchargedisobad03 event_rho )

##selections
cutEBdts="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03>0 && myphoton_catind<2"
cutEBdtb="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03<=0 && myphoton_catind<2"

cutEEdts="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03>0 && myphoton_catind>=2"
cutEEdtb="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03<=0 && myphoton_catind>=2"

cutEBpfs="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03>0 && myphoton_catind<2"
cutEBpfb="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03<=0 && myphoton_catind<2"

cutEEpfs="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03>0 && myphoton_catind>=2"
cutEEpfb="myphoton_pt>25. && myphoton_MITPre>0 && myphoton_MCTruthR03<=0 && myphoton_catind>=2"

######
signalTrainFiles="/publicfs/cms/data/Hgg/taojq/h2gglobe_V13_02/GluGluToHToGG_M-123/*.root/PhotonTree"
signalWeightTrainFiles="/publicfs/cms/data/Hgg/taojq/h2gglobe_V13_02/GluGluToHToGG_M-123_GJet2BinsXS_OnlyEtaWeights.root/PhotonTree"
SigWeight="myphoton_weightFirstEta"

#backgroundTrainFiles="/publicfs/cms/data/Hgg/taojq/h2gglobe_V13_02/GJet_Pt40_Rome_2_pf/Part1/*.root/PhotonTree"
#SigWeight="myphoton_weightFirstEta"
#BkgWeight="event_NewPUweight"
backgroundTrainFiles1="/publicfs/cms/data/Hgg/taojq/h2gglobe_V13_02/GJet_Pt-20to40_pf/*.root/PhotonTree"
backgroundTrainFiles2="/publicfs/cms/data/Hgg/taojq/h2gglobe_V13_02/GJet_Pt40_Rome_2_pf/*.root/PhotonTree"
BkgWeight="event_XSweight"

signalTestFiles="/publicfs/cms/data/Hgg/taojq/h2gglobe_V13_02/GluGluToHToGG_M-124/*.root/PhotonTree"
signalWeightTestFiles="/publicfs/cms/data/Hgg/taojq/h2gglobe_V13_02/GluGluToHToGG_M-124_GJet2BinsXS_OnlyEtaWeights.root/PhotonTree"


function printconf()
{
echo "# vim: ft=conf"
   echo "variable $1"
   echo "cut_signal $2"
   echo "cut_background $3"
#   echo "signal_train_n -1"
#   echo "background_train_n -1"
   echo "signal_trainfile $4"
   echo "signalWight_trainfile $5"
   echo "weightS $6"

   echo "background_file $7"
   echo "background_file $8"
   echo "weightB $9"

   echo "signal_testfile ${10}"
   echo "signalWight_testfile ${11}"
   #echo "background_trainfile $6"
   #echo "weightB $8"
   #echo "signal_file $8"
   #echo "signal_file ${10}"
   #echo "background_file ${11}"
   #echo "signal_train_n  $8"
   #echo "background_train_n  $9"

}



function myexpand()
{
   eval echo -n "$1"
}


#echo > batch-submit.sh
### parts


declare -i Remove
declare -i NUM
declare -i m
declare -i keep
declare -i p
declare -i k
declare -a AA
declare -a AAA
  
#for d in EBdt EEdt;do
for d in EBdt;do

 AAA=( `myexpand '${var'$d'[@]}'` )
 NUM=${#AAA[*]}

 Remove=-1
 nvar=${NUM}
 for((train=0; train<NUM-1; train++));do
 #for((train=0; train<1; train++));do
    keep=0
    for((n=0; n<nvar; n++));do
      
      if [ $n -eq ${Remove} ]; then
          continue
      else
         STR="${AAA[n]}"
         AA[keep]="$STR"
         echo ${AA[keep]}
         #keep=${keep}+1
         let $[ keep +=1 ]
      fi
    done           
  

    nvar=${#AA[*]}
    echo "nvar--"${nvar} 

    unset AAA
    for((q=0; q<nvar; q++)); do
       STR="${AA[q]}"
       AAA[q]="$STR"
    done


    cutns="cut${d}s"
    cutnb="cut${d}b"
    #[[ ${d:0:1} == "U" ]] && Q=cmsq || Q=cmssq
    Q=cmssq

    for (( i=-1; i<nvar; ++i )) ;do
    #for (( i=-1; i<0; ++i )) ;do
       A="";
       for (( j=0;   j<i;    ++j )) ;do A="$A ${AA[j]}, " ;done
       for (( j=i+1; j<nvar; ++j )) ;do A="$A ${AA[j]}, " ;done

   if [ $i -ge 0 ]; then
     if [ $i -le 9 ]; then 
       if [ ${train} -le 9 ]; then
         if [ "${AA[i]}" = "myphoton_ecalisodr03-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "myphoton_hcalisodr04-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "event_rho" ]; then
           varname="${AA[i]:6:3}"
         elif [ "${AA[i]}" = "myphoton_E5x5/myphoton_mustenergy" ]; then
            varname1="${AA[i]:9:4}"
            varname2="${AA[i]:23:10}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pfmustOutEnergy/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:10}"
            varname2="${AA[i]:34:6}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pflowE/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:6}"
            varname2="${AA[i]:25:6}"
            varname="${varname1}O${varname2}"
         else
            varname="${AA[i]:9}"
         fi
         outid=`printf "${d}0%iN-%s_0%i" ${i} ${varname} ${train}`
       fi

       if [ ${train} -gt 9 ]; then  
         if [ "${AA[i]}" = "myphoton_ecalisodr03-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "myphoton_hcalisodr04-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "event_rho" ]; then
           varname="${AA[i]:6:3}"
         elif [ "${AA[i]}" = "myphoton_E5x5/myphoton_mustenergy" ]; then
            varname1="${AA[i]:9:4}"
            varname2="${AA[i]:23:10}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pfmustOutEnergy/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:10}"
            varname2="${AA[i]:34:6}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pflowE/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:6}"
            varname2="${AA[i]:25:6}"
            varname="${varname1}O${varname2}"
         else
            varname="${AA[i]:9}"
         fi
         outid=`printf "${d}0%iN-%s_%i" ${i} ${varname} ${train}`
      fi
     fi
     
     if [ $i -gt 9 ]; then  
       if [ ${train} -le 9 ]; then
         if [ "${AA[i]}" = "myphoton_ecalisodr03-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "myphoton_hcalisodr04-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "event_rho" ]; then
           varname="${AA[i]:6:3}"
         elif [ "${AA[i]}" = "myphoton_E5x5/myphoton_mustenergy" ]; then
            varname1="${AA[i]:9:4}"
            varname2="${AA[i]:23:10}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pfmustOutEnergy/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:10}"
            varname2="${AA[i]:34:6}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pflowE/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:6}"
            varname2="${AA[i]:25:6}"
            varname="${varname1}O${varname2}"
         else
            varname="${AA[i]:9}"
         fi
         outid=`printf "${d}%iN-%s_0%i" ${i} ${varname} ${train}`
       fi
       if [ ${train} -gt 9 ]; then
         if [ "${AA[i]}" = "myphoton_ecalisodr03-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "myphoton_hcalisodr04-event_rho*0.17" ]; then
           varname="${AA[i]:9:11}"
         elif [ "${AA[i]}" = "event_rho" ]; then
           varname="${AA[i]:6:3}"
         elif [ "${AA[i]}" = "myphoton_E5x5/myphoton_mustenergy" ]; then
            varname1="${AA[i]:9:4}"
            varname2="${AA[i]:23:10}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pfmustOutEnergy/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:10}"
            varname2="${AA[i]:34:6}"
            varname="${varname1}O${varname2}"
         elif [ "${AA[i]}" = "myphoton_pflowE/myphoton_Escraw" ]; then
            varname1="${AA[i]:9:6}"
            varname2="${AA[i]:25:6}"
            varname="${varname1}O${varname2}"
         else
            varname="${AA[i]:9}"
         fi
         outid=`printf "${d}%iN-%s_%i" ${i} ${varname} ${train}`
      fi
     fi
    fi

 
    if [ $i -lt 0 ]; then
      if [ ${train} -le 9 ]; then
         outid=`printf "${d}_N_0%i" ${train}`
      else 
         outid=`printf "${d}_N_%i" ${train}`
      fi
    fi

  
     #echo "Taojq--Job: ${outid}"

     #printconf "${A}" "${!cutns}" "${!cutnb}" "$signalFiles" "$signalWeightFiles" "$backgroundFiles" "$SigWeight" "$BkgWeight" > $outid.conf
     printconf "${A}" "${!cutns}" "${!cutnb}" "$signalTrainFiles" "$signalWeightTrainFiles"  "$SigWeight"  "$backgroundTrainFiles1" "$backgroundTrainFiles2"  "$BkgWeight"  "$signalTestFiles" "$signalWeightTestFiles" > $outid.conf

     #echo qsub -q $Q -N $outid -d \$PWD -V '<<<' \"root -l -b -q TMVA.C"'(\\\""$outid"\\\")'"\" >> batch-submit.sh
    
  
    #root -b -q TMVA.C\(\"$outid\"\)    
    qsub -q cmssq -N $outid -d $PWD -V <<< "root -l -b -q TMVA.C'(\"$outid\")'"

    done

    for((w=0; w<100000000; w++)); do
      
      qstat -u fanjw > check.log
      if test -s check.log; then
         sleep 20
      else
         sleep 2m
         qstat -u fanjw > check2.log
         if test -s check2.log; then
            continue
         else
            break
         fi
      fi
    done

     
    if [ ${train} -le 9 ]; then
      Rankresult=${d}_Rank_0${train}.log
      root -l -b -q eff_cmp_efficiency.C data/${d}*N*_0${train}.root > ${Rankresult} 
    else
      Rankresult=${d}_Rank_${train}.log
      root -l -b -q eff_cmp_efficiency.C data/${d}*N*_${train}.root > ${Rankresult} 
    fi
   
 
    #eval `echo 'cat Rank_EBdt.log | grep MIN'`
    declare -i removeNum
    removeNum=`cat ${Rankresult} | grep MIN | awk -F ":" '{print $2}' | awk -F " " '{print $1}'`
    removeVar=`cat ${Rankresult} | grep MIN | awk -F ":" '{print $2}' | awk -F " " '{print $2}'`
    Remove=${removeNum}
    echo "Remove---" ${Remove} "---" ${removeVar}
    echo  

    ResultAll=${d}_Result.log
    cat ${Rankresult} | grep MIN | awk -F ":" '{print $2}' | awk -F " " '{print $2}' >> ${ResultAll}
   
    let last=NUM-2 
    if [ ${train} -eq ${last} ]; then
       for((l=0; l<2; l++));do
          if [ $l -eq $Remove ]; then
            continue
          else
            echo ${AA[l]}| awk -F "_" '{print $2}' >> ${ResultAll} 
          fi
       done
    fi

    unset AA


  done
done


