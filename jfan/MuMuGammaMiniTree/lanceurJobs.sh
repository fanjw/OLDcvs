#./bin/bash

###Data###
#for sample in 'Run2012AV2' #'Run2012AV2Run5' 'Run2012CpromptV1'
for sample in 'Run2012AV2Run5' 'Run2012CpromptV1'
do
	for Zgamma2 in '0'  ##'0' '1' '2' '3'
	do
		for SetOfCorrections2 in '1.0' ##'ETHZ' 'START42_V11'
		do
			for muonsCor in '0' ##'3' 
			do
			##for HighMmumuLim in '76' '84' ##36
			##do 
			##done		
			##for HighMmumuLim in '78' '82' '84' '86' '88' ##38
			##do
			##done
			##for HighMmumuLim in '80' ##40
			##do
			##done
			##for HighMmumuLim in '82' ##42
			##do
			##done
			##for HighMmumuLim in '84' ##44
			##do
			##done
			qsub batchJob.sh ${sample} ${sample} $Zgamma2 40 80 $SetOfCorrections2 1.0 $muonsCor
                        ##hadd miniTree_${sample}_NoMuonCor_40_80_v1_partALL.root miniTree_${sample}_NoMuonCor_40_80_v1_part[0-9]*root
                        #mv miniTree_${sample}_NoMuonCor_40_80_v1_part[0-9]*root stored_miniTree/
                        #mv ${sample}_NoMuonCor_40_80_v1_part[0-9]*err stored_errput/
                        #mv ${sample}_NoMuonCor_40_80_v1_part[0-9]*out stored_output/	



##			if [ "$SetOfCorrections2" = "0" ]
##                        then	
##				qsub batchJob.sh ${sample} ${sample}_NoMuonCor_40_80_v1 $Zgamma2 40 80 $SetOfCorrections2 1.0 $muonsCor
##				##hadd miniTree_${sample}_NoMuonCor_40_80_v1_partALL.root miniTree_${sample}_40_80_v1_part[0-9]*root
##				##mv miniTree_${sample}_NoMuonCor_40_80_v1_part[0-9]*root stored_miniTree/
##                        	##mv ${sample}_NoMuonCor_40_80_v1_part[0-9]*err stored_errput/
##                        	#"mv ${sample}_NoMuonCor_40_80_v1_part[0-9]*out stored_output/
##			fi
##			if [ "$SetOfCorrections2" = "3" ]
##                        then                
##                                qsub batchJob.sh ${sample} ${sample}_RochCor_40_80_v1 $Zgamma2 40 80 $SetOfCorrections2 1.0 $muonsCor
##                                ##hadd miniTree_${sample}_RochCor_40_80_v1_partALL.root miniTree_${sample}_40_80_v1_part[0-9]*root
##                                ##mv miniTree_${sample}_RochCor_40_80_v1_part[0-9]*root stored_miniTree/
##                                ##mv ${sample}_RochCor_40_80_v1_part[0-9]*err stored_errput/
##                                #"mv ${sample}_RochCor_40_80_v1_part[0-9]*out stored_output/
##                        fi	
			done
		done
	done
done

###MC###
###for sample in '' ''
###do
###        for Zgamma2 in '1' '2'  ##'0' '1' '2' '3'
###        do
###                for SetOfCorrections2 in 'MITregression' ##'ETHZ' 'START42_V11'
###                do
###                        for muonsCor in '0' '3' ##'3' 
###                        do  
###                        ##for HighMmumuLim in '76' '84' ##36
###                        ##do 
###                        ##done          
###                        ##for HighMmumuLim in '78' '82' '84' '86' '88' ##38
###                        ##do
###                        ##done
###                        ##for HighMmumuLim in '80' ##40
###                        ##do
###                        ##done
###                        ##for HighMmumuLim in '82' ##42
###                        ##do
###                        ##done
###                        ##for HighMmumuLim in '84' ##44
###                        ##do
###                        ##done
###                        if [ "$SetOfCorrections2" = "0" ]
###                        then                
###                                qsub batchJob.sh ${sample} ${sample}_NoMuonCor_${Zgamma2}_40_80_v1 $Zgamma2 40 80 $SetOfCorrections2 1.0 $muonsCor
###                                ##hadd miniTree_${sample}_NoMuonCor_${Zgamma2}_40_80_v1_partALL.root miniTree_${sample}_40_80_v1_part[0-9]*root
###                                ##mv miniTree_${sample}_NoMuonCor_${Zgamma2}_40_80_v1_part[0-9]*root stored_miniTree/
###                                ##mv ${sample}_NoMuonCor_${Zgamma2}_40_80_v1_part[0-9]*err stored_errput/
###                                #"mv ${sample}_NoMuonCor_${Zgamma2}_40_80_v1_part[0-9]*out stored_output/
###                        fi  
###                        if [ "$SetOfCorrections2" = "3" ]
###                        then
###                                qsub batchJob.sh ${sample} ${sample}_RochCor_${Zgamma2}_40_80_v1 $Zgamma2 40 80 $SetOfCorrections2 1.0 $muonsCor
###                                ##hadd miniTree_${sample}_RochCor_${Zgamma2}_40_80_v1_partALL.root miniTree_${sample}_40_80_v1_part[0-9]*root
###                                ##mv miniTree_${sample}_RochCor_${Zgamma2}_40_80_v1_part[0-9]*root stored_miniTree/
###                                ##mv ${sample}_RochCor_${Zgamma2}_40_80_v1_part[0-9]*err stored_errput/
###                                #"mv ${sample}_RochCor_${Zgamma2}_40_80_v1_part[0-9]*out stored_output/
###                        fi      
###                        done
###                done
###        done
###done
###


###./Selection_miniTree.exe ${sample} ${sample}_${SetOfCorrections2}_Rochester_40_80_5GeV_partALL 9999 -1 0 2011 PU_S6 40 80 ${SetOfCorrections2} 1.0 0

###./Selection_miniTree.exe Run2011A-ZMu-May10ReReco-v1 Run2011A-ZMu-May10ReReco-v1 9999 0 0 2011 PU_S6 40 80 ETHZ 1.0 3
