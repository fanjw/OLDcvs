#include "miniTreeMaker.h"

void beginMacro(){



	doHLT                    = true;
	doHLTobject		 = true;
  	doMC                     = true;
  	doJetMC                  = false;
  	doMETMC                  = false;
  	doPDFInfo                = true;
  	doSignalMuMuGamma        = false;
	doLeadingPhoton		 = true;
  	doSignalTopTop           = false;
  	doPhotonConversionMC     = true;
  	doElectronConversionMC   = false;
  	doBeamSpot               = true;
  	doPrimaryVertex          = true;
  	doZeePrimaryVertex       = false;
  	doTrack                  = true;
  	doJet                    = true;
  	doMuon                   = true;
  	doElectron               = true;
  	doPhoton                 = true;
  	doCluster                = true;
  	doPhotonConversion       = true;
  	doMET                    = true;
  	doBardak                 = false;
  	doPhotonVertexCorrection = false;
  	doPhotonIsolation        = true;
	doDiphotons		 = true;
	doWorstIsolation 	 = true;


	genJets = new TClonesArray("TRootParticle", 0);
	mcParticles = new TClonesArray("TRootMCParticle", 0);
	genMETs = new TClonesArray("TRootParticle", 0);
	mcPhotons = new TClonesArray("TRootMCPhoton", 0);
	mcElectrons = new TClonesArray("TRootMCElectron", 0);
	vertices = new TClonesArray("TRootVertex", 0);
	zeeVertices = new TClonesArray("TRootVertex", 0);
	tracks = new TClonesArray("TRootTrack", 0);
	jets = new TClonesArray("TRootJet", 0);
	pflowjets = new TClonesArray("TRootJet", 0);
	sisconejets = new TClonesArray("TRootJet", 0);
	muons = new TClonesArray("TRootMuon", 0);
	electrons = new TClonesArray("TRootElectron", 0);
	photons = new TClonesArray("TRootPhoton", 0);
	clusters = new TClonesArray("TRootCluster", 0);
	superClusters = new TClonesArray("TRootSuperCluster", 0);
	conversionTracks = new TClonesArray("TRootTrack", 0);
	met = new TClonesArray("TRootMET", 0);
	HLTObjects = new TClonesArray("TRootHLTObject", 0);


  inputEventTree->SetBranchAddress("Event", &event, &event_br);
  inputEventTree->SetBranchStatus("Event", 1);
 
  inputRunTree->SetBranchAddress("runInfos", &runInfos, &run_br);
  inputRunTree->SetBranchStatus("runInfos", 1);
  if(doMC)
    {
      inputEventTree->SetBranchAddress("MCParticles", &mcParticles, &mcParticles_br);
      inputEventTree->SetBranchStatus("MCParticles", 1);
    }
  
  if(doJetMC)
    {
      inputEventTree->SetBranchAddress("genJets", &genJets, &genJets_br);
      inputEventTree->SetBranchStatus("genJets", 1);
    }
  
  if(doMETMC)
    {
      inputEventTree->SetBranchAddress("genMETs", &genMETs, &genMETs_br);
      inputEventTree->SetBranchStatus("genMETs", 1);
    }
	
  if(doSignalMuMuGamma)
    {
      inputEventTree->SetBranchAddress("MuMuGamma", &mcMuMuGammaEvent, &mcSignalMuMuGamma_br);
      inputEventTree->SetBranchStatus("MuMuGamma", 1);
    }
  
  if(doSignalTopTop)
    {
      inputEventTree->SetBranchAddress("rootMCTopTop", &mcTopTopEvent, &mcTopTopEvent_br);
      inputEventTree->SetBranchStatus("rootMCTopTop", 1);
    }
	
  if(doPhotonConversionMC)
    {
      inputEventTree->SetBranchAddress("MCPhotons", &mcPhotons, &mcPhotons_br);
      inputEventTree->SetBranchStatus("MCPhotons", 1);
    }

  if (doElectronConversionMC)
    {
      inputEventTree->SetBranchAddress("MCElectrons",&mcElectrons, &mcElectrons_br);
      inputEventTree->SetBranchStatus("MCElectrons",1);
    }
  if(doBeamSpot)
    {
      inputEventTree->SetBranchAddress("BeamSpot", &beamSpot, &beamSpot_br);
      inputEventTree->SetBranchStatus("BeamSpot", 1);
    }

  if(doPrimaryVertex)
    {
      inputEventTree->SetBranchAddress("Vertices", &vertices, &vertices_br);
      inputEventTree->SetBranchStatus("Vertices", 1);
    }

  if(doZeePrimaryVertex)
    {
      inputEventTree->SetBranchAddress("ZeeVertices", &zeeVertices, &zeeVertices_br);
      inputEventTree->SetBranchStatus("ZeeVertices", 1);
    }
  
  if(doTrack)
    {
      inputEventTree->SetBranchAddress("Tracks", &tracks, &tracks_br);
      inputEventTree->SetBranchStatus("Tracks", 1);
    }


  if(doJet)	{
    inputEventTree->SetBranchAddress("ak5CaloJets", &jets, &jets_br);
    inputEventTree->SetBranchStatus("ak5CaloJets", 1);
  }


  if(doMuon)
    {
      inputEventTree->SetBranchAddress("muons", &muons, &muons_br);
      inputEventTree->SetBranchStatus("muons", 1);
    }
  
  if(doElectron)
    {
      inputEventTree->SetBranchAddress("gsfElectrons", &electrons, &electrons_br);
      inputEventTree->SetBranchStatus("gsfElectrons", 1);
    }
  
  if(doPhoton)
    {
      inputEventTree->SetBranchAddress("photons", &photons, &photons_br);
      inputEventTree->SetBranchStatus("photons", 1);
    }
  
  if(doCluster)
    {
      inputEventTree->SetBranchAddress("BasicClusters", &clusters, &clusters_br);
      inputEventTree->SetBranchStatus("BasicClusters", 1);
      
      inputEventTree->SetBranchAddress("SuperClusters", &superClusters, &superClusters_br);
      inputEventTree->SetBranchStatus("SuperClusters", 1);
    }
  

  if(doPhotonConversion)
    {
      inputEventTree->SetBranchAddress("ConversionTracks", &conversionTracks, &conversions_br);
      inputEventTree->SetBranchStatus("ConversionTracks", 1);
    }

  if(doMET)
    {
      inputEventTree->SetBranchAddress("met", &met, &met_br);
      inputEventTree->SetBranchStatus("met", 1);
    }
  
  if(doBardak)
    {
      inputEventTree->SetBranchAddress("bardak", &bardak, &bardak_br);
      inputEventTree->SetBranchStatus("bardak", 1);
    }
  if (doHLTobject)
    {
      inputEventTree->SetBranchAddress("HLTObjects",&HLTObjects, &HLTObjects_br);
      inputEventTree->SetBranchStatus("HLTObjects", 1);
    }



                myTree_ = new TTree("diPhotons","DiPhotonsInfos");
                
                myTree_->Branch("event_number",&event_number,"event_number/I");
                myTree_->Branch("event_runNumber",&event_runNumber, "event_runNumber/I");
                myTree_->Branch("event_LumiSection",&event_LumiSection,"event_LumiSection/I");
                myTree_->Branch("event_eventPtHat",&event_eventPtHat, "event_eventPtHat/F");
                myTree_->Branch("event_processId",&event_processId, "event_processId/I");
                myTree_->Branch("event_nRecoVertex",&event_nRecoVertex,"event_nRecoVertex/I");
                myTree_->Branch("event_nGenInTimeVertex",&event_nGenInTimeVertex,"event_nGenInTimeVertex/I");
                myTree_->Branch("event_nGenOutOfTimeVertex",&event_nGenOutOfTimeVertex,"event_nGenOutOfTimeVertex/I");
                myTree_->Branch("event_nPhotons",&event_nPhotons,"event_nPhotons/I");
                myTree_->Branch("event_rho",&event_rho,"event_rho/F");              	
                myTree_->Branch("event_nPairs",&event_nPairs,"event_nPairs/I");              	


                myTree_->Branch("dipho_mgg","std::vector<float>",&Pdipho_mgg);
                myTree_->Branch("dipho_qt","std::vector<float>",&Pdipho_qt);                                 
                myTree_->Branch("dipho_ql","std::vector<float>",&Pdipho_ql);                                 
                myTree_->Branch("dipho_deltaR","std::vector<float>",&Pdipho_deltaR);                     
                myTree_->Branch("dipho_costhetastar","std::vector<float>",&Pdipho_costhetastar);   
                myTree_->Branch("dipho_eta","std::vector<float>",&Pdipho_eta);        
                myTree_->Branch("dipho_etastar","std::vector<float>",&Pdipho_etastar);         
             	myTree_->Branch("diphoMC_mgg","std::vector<float>",&PdiphoMC_mgg);
                myTree_->Branch("diphoMC_qt","std::vector<float>",&PdiphoMC_qt);                                 
                myTree_->Branch("diphoMC_ql","std::vector<float>",&PdiphoMC_ql);                                 
                myTree_->Branch("diphoMC_deltaR","std::vector<float>",&PdiphoMC_deltaR);                     
                myTree_->Branch("diphoMC_costhetastar","std::vector<float>",&PdiphoMC_costhetastar); 
                myTree_->Branch("diphoMC_eta","std::vector<float>",&PdiphoMC_eta);      
                myTree_->Branch("diphoMC_etastar","std::vector<float>",&PdiphoMC_etastar);        
	
	
                myTree_->Branch("pholead_event_processId","std::vector<int>",&pholead_event_processId); 
                myTree_->Branch("pholead_isMatchingWithMC","std::vector<int>",&pholead_isMatchingWithMC); 
                myTree_->Branch("pholead_GenId","std::vector<int>",&pholead_GenId); 
                myTree_->Branch("pholead_MotherId","std::vector<int>",&pholead_MotherId); 
                myTree_->Branch("pholead_isPromptGenPho","std::vector<int>",&pholead_isPromptGenPho); 
                myTree_->Branch("pholead_isFromQuarkGen","std::vector<int>",&pholead_isFromQuarkGen); 
                myTree_->Branch("pholead_isPi0Gen","std::vector<int>",&pholead_isPi0Gen); 
                myTree_->Branch("pholead_isEtaGen","std::vector<int>",&pholead_isEtaGen); 
                myTree_->Branch("pholead_isRhoGen","std::vector<int>",&pholead_isRhoGen); 
                myTree_->Branch("pholead_isOmegaGen","std::vector<int>",&pholead_isOmegaGen); 
                myTree_->Branch("pholead_isGenElectron","std::vector<int>",&pholead_isGenElectron); 
                myTree_->Branch("pholead_eventPassHLT_Photon10_L1R","std::vector<int>",&pholead_eventPassHLT_Photon10_L1R); 
                myTree_->Branch("pholead_eventPassHLT_Photon15_L1R","std::vector<int>",&pholead_eventPassHLT_Photon15_L1R); 
                myTree_->Branch("pholead_eventPassHLT_DoublePhoton10_L1R","std::vector<int>",&pholead_eventPassHLT_DoublePhoton10_L1R); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus1_cone02","std::vector<float>",&pholead_PromptGenIsoEnergyStatus1_cone02); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus2_cone02","std::vector<float>",&pholead_PromptGenIsoEnergyStatus2_cone02); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus1_cone03","std::vector<float>",&pholead_PromptGenIsoEnergyStatus1_cone03); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus2_cone03","std::vector<float>",&pholead_PromptGenIsoEnergyStatus2_cone03); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus1_cone035","std::vector<float>",&pholead_PromptGenIsoEnergyStatus1_cone035); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus2_cone035","std::vector<float>",&pholead_PromptGenIsoEnergyStatus2_cone035); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus1_cone04","std::vector<float>",&pholead_PromptGenIsoEnergyStatus1_cone04); 
                myTree_->Branch("pholead_PromptGenIsoEnergyStatus2_cone04","std::vector<float>",&pholead_PromptGenIsoEnergyStatus2_cone04); 
                myTree_->Branch("pholead_trueE","std::vector<float>",&pholead_trueE); 
                myTree_->Branch("pholead_truePx","std::vector<float>",&pholead_truePx); 
                myTree_->Branch("pholead_truePy","std::vector<float>",&pholead_truePy); 
                myTree_->Branch("pholead_truePz","std::vector<float>",&pholead_truePz); 
                myTree_->Branch("pholead_trueEta","std::vector<float>",&pholead_trueEta); 
                myTree_->Branch("pholead_truePhi","std::vector<float>",&pholead_truePhi); 
                myTree_->Branch("pholead_MCisConverted","std::vector<int>",&pholead_MCisConverted); 
                myTree_->Branch("pholead_MCconvEoverP","std::vector<float>",&pholead_MCconvEoverP); 
                myTree_->Branch("pholead_MCconvCotanTheta","std::vector<float>",&pholead_MCconvCotanTheta); 
                myTree_->Branch("pholead_MCconvVertexX","std::vector<float>",&pholead_MCconvVertexX); 
                myTree_->Branch("pholead_MCconvVertexY","std::vector<float>",&pholead_MCconvVertexY); 
                myTree_->Branch("pholead_MCconvVertexZ","std::vector<float>",&pholead_MCconvVertexZ); 
                myTree_->Branch("pholead_eleMCtruthBrem","std::vector<float>",&pholead_eleMCtruthBrem); 
                myTree_->Branch("pholead_eleMCtruthNBrem","std::vector<int>",&pholead_eleMCtruthNBrem); 

		myTree_->Branch("pholead_et","std::vector<float>",&pholead_et); 
		myTree_->Branch("pholead_eta","std::vector<float>",&pholead_eta);   
		myTree_->Branch("pholead_SCeta","std::vector<float>",&pholead_SCeta);
		myTree_->Branch("pholead_r9","std::vector<float>",&pholead_r9);                                             
		myTree_->Branch("pholead_cPP","std::vector<float>",&pholead_cPP);
		myTree_->Branch("pholead_cEP","std::vector<float>",&pholead_cEP);
		myTree_->Branch("pholead_cEE","std::vector<float>",&pholead_cEE);
		myTree_->Branch("pholead_r19","std::vector<float>",&pholead_r19);                                             
		myTree_->Branch("pholead_SCEraw","std::vector<float>",&pholead_SCEraw);                                             
		myTree_->Branch("pholead_eMax","std::vector<float>",&pholead_eMax);                                             
		myTree_->Branch("pholead_e2x2","std::vector<float>",&pholead_e2x2);                                             
		myTree_->Branch("pholead_e5x5","std::vector<float>",&pholead_e5x5);                                             
		myTree_->Branch("pholead_ratioSeed","std::vector<float>",&pholead_ratioSeed);
		myTree_->Branch("pholead_ratioS4","std::vector<float>",&pholead_ratioS4);
		myTree_->Branch("pholead_lambdaRatio","std::vector<float>",&pholead_lambdaRatio);
		myTree_->Branch("pholead_lamdbaDivCov","std::vector<float>",&pholead_lamdbaDivCov);
		myTree_->Branch("pholead_secondMomentMaj","std::vector<float>",&pholead_secondMomentMaj);
		myTree_->Branch("pholead_secondMomentMin","std::vector<float>",&pholead_secondMomentMin);
		myTree_->Branch("pholead_secondMomentAlpha","std::vector<float>",&pholead_secondMomentAlpha);
		myTree_->Branch("pholead_covAngle","std::vector<float>",&pholead_covAngle);
		myTree_->Branch("pholead_covAngle2","std::vector<float>",&pholead_covAngle2);
		myTree_->Branch("pholead_S9overS9minusS1S2","std::vector<float>",&pholead_S9overS9minusS1S2);
		myTree_->Branch("pholead_etawidth","std::vector<float>",&pholead_etawidth);                              
		myTree_->Branch("pholead_phiwidth","std::vector<float>",&pholead_phiwidth);                              
		myTree_->Branch("pholead_sigieta","std::vector<float>",&pholead_sigieta);                                 
		myTree_->Branch("pholead_SCbr","std::vector<float>",&pholead_SCbr);
		myTree_->Branch("pholead_seedEnergy","std::vector<float>",&pholead_seedEnergy);
		myTree_->Branch("pholead_seedTime","std::vector<float>",&pholead_seedTime);
		myTree_->Branch("pholead_SCphi","std::vector<float>",&pholead_SCphi);
		myTree_->Branch("pholead_SCEtraw","std::vector<float>",&pholead_SCEtraw);
		myTree_->Branch("pholead_SCEraw","std::vector<float>",&pholead_SCEraw);
		myTree_->Branch("pholead_SCEt","std::vector<float>",&pholead_SCEt);
		myTree_->Branch("pholead_SCr9","std::vector<float>",&pholead_SCr9);
		myTree_->Branch("pholead_SCnbBC","std::vector<int>",&pholead_SCnbBC);
		myTree_->Branch("pholead_SCnXtal","std::vector<int>",&pholead_SCnXtal);
		myTree_->Branch("pholead_isMatchingWithHLTObject","std::vector<int>",&pholead_isMatchingWithHLTObject);
		myTree_->Branch("pholead_isConverted","std::vector<int>",&pholead_isConverted);
		myTree_->Branch("pholead_NtrackConv","std::vector<int>",&pholead_NtrackConv);
		myTree_->Branch("pholead_convEoverP","std::vector<float>",&pholead_convEoverP);
		myTree_->Branch("pholead_convMass","std::vector<float>",&pholead_convMass);
		myTree_->Branch("pholead_convCotanTheta","std::vector<float>",&pholead_convCotanTheta);
		myTree_->Branch("pholead_MCconvMass","std::vector<float>",&pholead_MCconvMass);
		myTree_->Branch("pholead_convLikely","std::vector<float>",&pholead_convLikely);
		myTree_->Branch("pholead_convVertexX","std::vector<float>",&pholead_convVertexX);
		myTree_->Branch("pholead_convVertexY","std::vector<float>",&pholead_convVertexY);
		myTree_->Branch("pholead_convVertexZ","std::vector<float>",&pholead_convVertexZ);
		myTree_->Branch("pholead_fBrem","std::vector<float>",&pholead_fBrem);
		myTree_->Branch("pholead_isAspike","std::vector<int>",&pholead_isAspike);
		myTree_->Branch("pholead_xVertex","std::vector<float>",&pholead_xVertex);
		myTree_->Branch("pholead_yVertex","std::vector<float>",&pholead_yVertex);
		myTree_->Branch("pholead_zVertex","std::vector<float>",&pholead_zVertex);
 
                myTree_->Branch("pholead_HcalIso","std::vector<float>",&pholead_HcalIso);                                     
		myTree_->Branch("pholead_EcalIso","std::vector<float>",&pholead_EcalIso);
		myTree_->Branch("pholead_TrackerIso","std::vector<float>",&pholead_TrackerIso);
                myTree_->Branch("pholead_HcalIso_MIT03","std::vector<float>",&pholead_HcalIso_MIT03);                                     
		myTree_->Branch("pholead_EcalIso_MIT03","std::vector<float>",&pholead_EcalIso_MIT03);
		myTree_->Branch("pholead_TrackerIso_MIT03","std::vector<float>",&pholead_TrackerIso_MIT03);
		myTree_->Branch("pholead_HcalEcal_MIT03","std::vector<float>",&pholead_HcalEcal_MIT03);
		myTree_->Branch("pholead_AbsTrackerIso","std::vector<float>",&pholead_AbsTrackerIso);
		myTree_->Branch("pholead_HcalIsodR03","std::vector<float>",&pholead_HcalIsodR03);             
		myTree_->Branch("pholead_EcalIsodR03","std::vector<float>",&pholead_EcalIsodR03);
		myTree_->Branch("pholead_TrackerIsodR03","std::vector<float>",&pholead_TrackerIsodR03);
		myTree_->Branch("pholead_hoe","std::vector<float>",&pholead_hoe);

                myTree_->Branch("pholead_Cat","std::vector<int>",&pholead_Cat);
                myTree_->Branch("pholead_HasPixSeed","std::vector<int>",&pholead_HasPixSeed);
                myTree_->Branch("pholead_seedSeverity","std::vector<int>",&pholead_seedSeverity);
                myTree_->Branch("pholead_recoFlag","std::vector<int>",&pholead_recoFlag);
		myTree_->Branch("pholead_isEB","std::vector<int>",&pholead_isEB);
		myTree_->Branch("pholead_isEE","std::vector<int>",&pholead_isEE);








		myTree_->Branch("photrail_event_processId","std::vector<int>",&photrail_event_processId);
		myTree_->Branch("photrail_isMatchingWithMC","std::vector<int>",&photrail_isMatchingWithMC);
                myTree_->Branch("photrail_GenId","std::vector<int>",&photrail_GenId); 
                myTree_->Branch("photrail_MotherId","std::vector<int>",&photrail_MotherId); 
                myTree_->Branch("photrail_isPromptGenPho","std::vector<int>",&photrail_isPromptGenPho); 
                myTree_->Branch("photrail_isFromQuarkGen","std::vector<int>",&photrail_isFromQuarkGen); 
                myTree_->Branch("photrail_isPi0Gen","std::vector<int>",&photrail_isPi0Gen); 
                myTree_->Branch("photrail_isEtaGen","std::vector<int>",&photrail_isEtaGen); 
                myTree_->Branch("photrail_isRhoGen","std::vector<int>",&photrail_isRhoGen); 
                myTree_->Branch("photrail_isOmegaGen","std::vector<int>",&photrail_isOmegaGen); 
                myTree_->Branch("photrail_isGenElectron","std::vector<int>",&photrail_isGenElectron); 
                myTree_->Branch("photrail_eventPassHLT_Photon10_L1R","std::vector<int>",&photrail_eventPassHLT_Photon10_L1R); 
                myTree_->Branch("photrail_eventPassHLT_Photon15_L1R","std::vector<int>",&photrail_eventPassHLT_Photon15_L1R); 
                myTree_->Branch("photrail_eventPassHLT_DoublePhoton10_L1R","std::vector<int>",&photrail_eventPassHLT_DoublePhoton10_L1R); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus1_cone02","std::vector<float>",&photrail_PromptGenIsoEnergyStatus1_cone02); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus2_cone02","std::vector<float>",&photrail_PromptGenIsoEnergyStatus2_cone02); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus1_cone03","std::vector<float>",&photrail_PromptGenIsoEnergyStatus1_cone03); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus2_cone03","std::vector<float>",&photrail_PromptGenIsoEnergyStatus2_cone03); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus1_cone035","std::vector<float>",&photrail_PromptGenIsoEnergyStatus1_cone035); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus2_cone035","std::vector<float>",&photrail_PromptGenIsoEnergyStatus2_cone035); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus1_cone04","std::vector<float>",&photrail_PromptGenIsoEnergyStatus1_cone04); 
                myTree_->Branch("photrail_PromptGenIsoEnergyStatus2_cone04","std::vector<float>",&photrail_PromptGenIsoEnergyStatus2_cone04); 
                myTree_->Branch("photrail_trueE","std::vector<float>",&photrail_trueE); 
                myTree_->Branch("photrail_truePx","std::vector<float>",&photrail_truePx); 
                myTree_->Branch("photrail_truePy","std::vector<float>",&photrail_truePy); 
                myTree_->Branch("photrail_truePz","std::vector<float>",&photrail_truePz); 
                myTree_->Branch("photrail_trueEta","std::vector<float>",&photrail_trueEta); 
                myTree_->Branch("photrail_truePhi","std::vector<float>",&photrail_truePhi); 
                myTree_->Branch("photrail_MCisConverted","std::vector<int>",&photrail_MCisConverted); 
                myTree_->Branch("photrail_MCconvEoverP","std::vector<float>",&photrail_MCconvEoverP); 
                myTree_->Branch("photrail_MCconvCotanTheta","std::vector<float>",&photrail_MCconvCotanTheta); 
                myTree_->Branch("photrail_MCconvVertexX","std::vector<float>",&photrail_MCconvVertexX); 
                myTree_->Branch("photrail_MCconvVertexY","std::vector<float>",&photrail_MCconvVertexY); 
                myTree_->Branch("photrail_MCconvVertexZ","std::vector<float>",&photrail_MCconvVertexZ); 
                myTree_->Branch("photrail_eleMCtruthBrem","std::vector<float>",&photrail_eleMCtruthBrem); 
                myTree_->Branch("photrail_eleMCtruthNBrem","std::vector<int>",&photrail_eleMCtruthNBrem); 
		
                myTree_->Branch("photrail_et","std::vector<float>",&photrail_et); 
		myTree_->Branch("photrail_eta","std::vector<float>",&photrail_eta);   
		myTree_->Branch("photrail_SCeta","std::vector<float>",&photrail_SCeta);
		myTree_->Branch("photrail_r9","std::vector<float>",&photrail_r9);                                             
		myTree_->Branch("photrail_cPP","std::vector<float>",&photrail_cPP);
		myTree_->Branch("photrail_cEP","std::vector<float>",&photrail_cEP);
		myTree_->Branch("photrail_cEE","std::vector<float>",&photrail_cEE);
		myTree_->Branch("photrail_r19","std::vector<float>",&photrail_r19);                                             
		myTree_->Branch("photrail_SCEraw","std::vector<float>",&photrail_SCEraw);                                             
		myTree_->Branch("photrail_eMax","std::vector<float>",&photrail_eMax);                                             
		myTree_->Branch("photrail_e2x2","std::vector<float>",&photrail_e2x2);                                             
		myTree_->Branch("photrail_e5x5","std::vector<float>",&photrail_e5x5);                                             
		myTree_->Branch("photrail_ratioSeed","std::vector<float>",&photrail_ratioSeed);
		myTree_->Branch("photrail_ratioS4","std::vector<float>",&photrail_ratioS4);
		myTree_->Branch("photrail_lambdaRatio","std::vector<float>",&photrail_lambdaRatio);
		myTree_->Branch("photrail_lamdbaDivCov","std::vector<float>",&photrail_lamdbaDivCov);
		myTree_->Branch("photrail_secondMomentMaj","std::vector<float>",&photrail_secondMomentMaj);
		myTree_->Branch("photrail_secondMomentMin","std::vector<float>",&photrail_secondMomentMin);
		myTree_->Branch("photrail_secondMomentAlpha","std::vector<float>",&photrail_secondMomentAlpha);
		myTree_->Branch("photrail_covAngle","std::vector<float>",&photrail_covAngle);
		myTree_->Branch("photrail_covAngle2","std::vector<float>",&photrail_covAngle2);
		myTree_->Branch("photrail_S9overS9minusS1S2","std::vector<float>",&photrail_S9overS9minusS1S2);
		myTree_->Branch("photrail_etawidth","std::vector<float>",&photrail_etawidth);                              
		myTree_->Branch("photrail_phiwidth","std::vector<float>",&photrail_phiwidth);                              
		myTree_->Branch("photrail_sigieta","std::vector<float>",&photrail_sigieta);                                 
		myTree_->Branch("photrail_SCbr","std::vector<float>",&photrail_SCbr);
		myTree_->Branch("photrail_seedEnergy","std::vector<float>",&photrail_seedEnergy);
		myTree_->Branch("photrail_seedTime","std::vector<float>",&photrail_seedTime);
		myTree_->Branch("photrail_SCphi","std::vector<float>",&photrail_SCphi);
		myTree_->Branch("photrail_SCEtraw","std::vector<float>",&photrail_SCEtraw);
		myTree_->Branch("photrail_SCEraw","std::vector<float>",&photrail_SCEraw);
		myTree_->Branch("photrail_SCEt","std::vector<float>",&photrail_SCEt);
		myTree_->Branch("photrail_SCr9","std::vector<float>",&photrail_SCr9);
		myTree_->Branch("photrail_SCnbBC","std::vector<int>",&photrail_SCnbBC);
		myTree_->Branch("photrail_SCnXtal","std::vector<int>",&photrail_SCnXtal);
		myTree_->Branch("photrail_isMatchingWithHLTObject","std::vector<int>",&photrail_isMatchingWithHLTObject);
		myTree_->Branch("photrail_isConverted","std::vector<int>",&photrail_isConverted);
		myTree_->Branch("photrail_NtrackConv","std::vector<int>",&photrail_NtrackConv);
		myTree_->Branch("photrail_convEoverP","std::vector<float>",&photrail_convEoverP);
		myTree_->Branch("photrail_convMass","std::vector<float>",&photrail_convMass);
		myTree_->Branch("photrail_convCotanTheta","std::vector<float>",&photrail_convCotanTheta);
		myTree_->Branch("photrail_MCconvMass","std::vector<float>",&photrail_MCconvMass);
		myTree_->Branch("photrail_convLikely","std::vector<float>",&photrail_convLikely);
		myTree_->Branch("photrail_convVertexX","std::vector<float>",&photrail_convVertexX);
		myTree_->Branch("photrail_convVertexY","std::vector<float>",&photrail_convVertexY);
		myTree_->Branch("photrail_convVertexZ","std::vector<float>",&photrail_convVertexZ);
		myTree_->Branch("photrail_fBrem","std::vector<float>",&photrail_fBrem);
		myTree_->Branch("photrail_isAspike","std::vector<int>",&photrail_isAspike);
		myTree_->Branch("photrail_xVertex","std::vector<float>",&photrail_xVertex);
		myTree_->Branch("photrail_yVertex","std::vector<float>",&photrail_yVertex);
		myTree_->Branch("photrail_zVertex","std::vector<float>",&photrail_zVertex);

                myTree_->Branch("photrail_HcalIso","std::vector<float>",&photrail_HcalIso);                                     
		myTree_->Branch("photrail_EcalIso","std::vector<float>",&photrail_EcalIso);
		myTree_->Branch("photrail_TrackerIso","std::vector<float>",&photrail_TrackerIso);
                myTree_->Branch("photrail_HcalIso_MIT03","std::vector<float>",&photrail_HcalIso_MIT03);                                     
		myTree_->Branch("photrail_EcalIso_MIT03","std::vector<float>",&photrail_EcalIso_MIT03);
		myTree_->Branch("photrail_TrackerIso_MIT03","std::vector<float>",&photrail_TrackerIso_MIT03);
		myTree_->Branch("photrail_HcalEcal_MIT03","std::vector<float>",&photrail_HcalEcal_MIT03);
		myTree_->Branch("photrail_AbsTrackerIso","std::vector<float>",&photrail_AbsTrackerIso);
		myTree_->Branch("photrail_HcalIsodR03","std::vector<float>",&photrail_HcalIsodR03);             
		myTree_->Branch("photrail_EcalIsodR03","std::vector<float>",&photrail_EcalIsodR03);
		myTree_->Branch("photrail_TrackerIsodR03","std::vector<float>",&photrail_TrackerIsodR03);
		myTree_->Branch("photrail_hoe","std::vector<float>",&photrail_hoe);
                
                myTree_->Branch("photrail_Cat","std::vector<int>",&photrail_Cat);
                myTree_->Branch("photrail_HasPixSeed","std::vector<int>",&photrail_HasPixSeed);
                myTree_->Branch("photrail_seedSeverity","std::vector<int>",&photrail_seedSeverity);
                myTree_->Branch("photrail_recoFlag","std::vector<int>",&photrail_recoFlag);
		myTree_->Branch("photrail_isEB","std::vector<int>",&photrail_isEB);
		myTree_->Branch("photrail_isEE","std::vector<int>",&photrail_isEE);


}




void saveThisEvent(TRootEvent *theEvent, pair <TRootPhoton*, TRootPhoton*> theDiphotonPair, TClonesArray *thePhotonArray, TClonesArray *theVerticeArray){
       
           //for data
		event_number = theEvent->eventId();
		event_runNumber = theEvent->runId();
		event_LumiSection =theEvent->luminosityBlock();
		event_rho = theEvent->rho();
		event_nRecoVertex = theVerticeArray->GetEntriesFast();
		event_nPhotons = thePhotonArray->GetEntriesFast();
                event_nPairs = ngoodPairs; 
           //for MC
                event_eventPtHat = theEvent->ptHat();
		event_processId = theEvent->processID();
                //event_nGenInTimeVertex = theEvent->inTimePU_NumInteractions();
                event_nGenInTimeVertex = theEvent->nInTimePUVertices();
                //event_nGenOutOfTimeVertex = theEvent->outOfTimePU_NumInteractions();                                         
                event_nGenOutOfTimeVertex = theEvent->nOOTPUVertices();                                         



                // now calculation of the dipho kine variable 
		TLorentzVector Plead, Ptrail, Psum;
		Plead.SetPxPyPzE((theDiphotonPair.first)->Px(),(theDiphotonPair.first)->Py(),(theDiphotonPair.first)->Pz(),(theDiphotonPair.first)->Energy());
		Ptrail.SetPxPyPzE((theDiphotonPair.second)->Px(),(theDiphotonPair.second)->Py(),(theDiphotonPair.second)->Pz(),(theDiphotonPair.second)->Energy());
		Psum = Plead + Ptrail;
		Pdipho_mgg -> push_back(Psum.M());
		Pdipho_qt -> push_back(Psum.Pt());
		Pdipho_ql -> push_back(Psum.Pz());
		Pdipho_deltaR -> push_back(DeltaR((theDiphotonPair.first)->Phi(),(theDiphotonPair.second)->Phi(),(theDiphotonPair.first)->Eta(),(theDiphotonPair.second)->Eta()));
		Pdipho_costhetastar -> push_back(fabs(CosThetaStar(Plead,Ptrail)));
		Pdipho_eta -> push_back(Psum.Eta());
		Pdipho_etastar -> push_back(1.0*(Plead.Eta()-Ptrail.Eta())/2);
		int catDipho = findTheDiphoCat(*(theDiphotonPair.first), *(theDiphotonPair.second));
		cout << "passing pair : Run = " << theEvent->runId() << " LS = " << theEvent->luminosityBlock() << " Event = " << theEvent->eventId() << " SelVtx = 0  CAT4 = " << catDipho << endl;


		TRootMCParticle theLeadMC, theTrailMC;
		TLorentzVector PleadMC, PtrailMC, PsumMC;		
		if (findGenParticle(theDiphotonPair.first, mcParticles, &theLeadMC)) Ppholead_isMatchingWithMC -> push_back(1);
		else Ppholead_isMatchingWithMC -> push_back(0);

		if (findGenParticle(theDiphotonPair.second, mcParticles, &theTrailMC)) Pphotrail_isMatchingWithMC -> push_back(1);
		else Pphotrail_isMatchingWithMC -> push_back(0);


		if ((findGenParticle(theDiphotonPair.first, mcParticles, &theLeadMC))&&findGenParticle(theDiphotonPair.second, mcParticles, &theTrailMC)){
			PleadMC.SetPxPyPzE(theLeadMC.Px(),theLeadMC.Py(),theLeadMC.Pz(),theLeadMC.Energy());
			PtrailMC.SetPxPyPzE(theTrailMC.Px(),theTrailMC.Py(),theTrailMC.Pz(),theTrailMC.Energy());
			PsumMC = PleadMC + PtrailMC;
			PdiphoMC_mgg -> push_back(PsumMC.M());
	                PdiphoMC_qt -> push_back(PsumMC.Pt());
        	        PdiphoMC_ql -> push_back(PsumMC.Pz());
                	PdiphoMC_deltaR -> push_back(DeltaR(theLeadMC.Phi(),theTrailMC.Phi(),theLeadMC.Eta(),theTrailMC.Eta()));
               		PdiphoMC_costhetastar -> push_back(fabs(CosThetaStar(PleadMC,PtrailMC)));
                	PdiphoMC_eta -> push_back(PsumMC.Eta());
                	PdiphoMC_etastar -> push_back(1.0*(PleadMC.Eta()-PtrailMC.Eta())/2);
		}



               
              //now lead photon
                if (doMC) {
                    int lead_GenId=-1;int lead_MotherId=-1;int lead_isGenElectron=-1;int lead_isPromptGenPho=-1;int lead_isFromQuarkGen=-1;int lead_isPi0Gen=-1;int lead_isEtaGen=-1;int lead_isRhoGen=-1;int lead_isOmegaGen=-1;float lead_PromptGenIsoEnergyStatus1_cone02=0;float lead_PromptGenIsoEnergyStatus2_cone02=0;float lead_trueE=0;float lead_truePx=0;float lead_truePy=0;float lead_truePz=0;float lead_trueEta=0;float lead_truePhi=0;
		    

                    doGenInfo(theDiphotonPair.first, mcParticles, &(lead_GenId), &(lead_MotherId), &(lead_isGenElectron), &(lead_isPromptGenPho), &(lead_isFromQuarkGen), &(lead_isPi0Gen), &(lead_isEtaGen), &(lead_isRhoGen), &(lead_isOmegaGen), &(lead_PromptGenIsoEnergyStatus1_cone02), &(lead_PromptGenIsoEnergyStatus2_cone02),&(lead_trueE),&(lead_truePx),&(lead_truePy),&(lead_truePz),&(lead_trueEta),&(lead_truePhi), 0.2);
                    Ppholead_GenId->push_back(lead_GenId); 
                    Ppholead_MotherId->push_back(lead_MotherId); 
                    Ppholead_isGenElectron->push_back(lead_isGenElectron); 
                    Ppholead_isPromptGenPho->push_back(lead_isPromptGenPho); 
                    Ppholead_isFromQuarkGen->push_back(lead_isFromQuarkGen); 
                    Ppholead_isPi0Gen->push_back(lead_isPi0Gen); 
                    Ppholead_isEtaGen->push_back(lead_isEtaGen); 
                    Ppholead_isRhoGen->push_back(lead_isRhoGen); 
                    Ppholead_isOmegaGen->push_back(lead_isOmegaGen); 
                    Ppholead_PromptGenIsoEnergyStatus1_cone02->push_back(lead_PromptGenIsoEnergyStatus1_cone02); 
                    Ppholead_PromptGenIsoEnergyStatus2_cone02->push_back(lead_PromptGenIsoEnergyStatus2_cone02); 
                    Ppholead_trueE->push_back(lead_trueE); 
                    Ppholead_truePx->push_back(lead_truePx); 
                    Ppholead_truePy->push_back(lead_truePy); 
                    Ppholead_truePz->push_back(lead_truePz); 
                    Ppholead_trueEta->push_back(lead_trueEta); 
                    Ppholead_truePhi->push_back(lead_truePhi); 
                }

             

                Ppholead_event_processId -> push_back(theEvent->processID());
                Ppholead_et -> push_back((theDiphotonPair.first)->Et());
                Ppholead_eta -> push_back((theDiphotonPair.first)->Eta());
                Ppholead_SCeta -> push_back((theDiphotonPair.first)->superCluster()->Eta()); 
		Ppholead_r9 -> push_back((theDiphotonPair.first)->r9());
		Ppholead_cPP -> push_back((theDiphotonPair.first)->covPhiPhi());
		Ppholead_cEE -> push_back((theDiphotonPair.first)->covEtaEta());
		Ppholead_cEP -> push_back((theDiphotonPair.first)->covEtaPhi());

		Ppholead_r19 -> push_back((theDiphotonPair.first)->r19());
		Ppholead_SCEraw -> push_back((theDiphotonPair.first)->superCluster()->rawEnergy());
		Ppholead_eMax -> push_back((theDiphotonPair.first)->eMax());
		Ppholead_e2x2 -> push_back((theDiphotonPair.first)->e2x2());
		Ppholead_e5x5 -> push_back((theDiphotonPair.first)->e5x5());

		if ((theDiphotonPair.first)->superCluster()->rawEnergy()!=0) Ppholead_ratioSeed -> push_back((theDiphotonPair.first)->eMax()/(theDiphotonPair.first)->superCluster()->rawEnergy()); else Ppholead_ratioSeed -> push_back(0);

		if ((theDiphotonPair.first)->e5x5()!=0) Ppholead_ratioS4 -> push_back((theDiphotonPair.first)->e2x2() / (theDiphotonPair.first)->e5x5()); else Ppholead_ratioS4 -> push_back(0);

		if ( ((theDiphotonPair.first)->covEtaEta()+(theDiphotonPair.first)->covPhiPhi()+sqrt(((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())*((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())+4*(theDiphotonPair.first)->covEtaPhi()*(theDiphotonPair.first)->covEtaPhi()))!=0)  Ppholead_lambdaRatio -> push_back(((theDiphotonPair.first)->covEtaEta()+(theDiphotonPair.first)->covPhiPhi()-sqrt(((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())*((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())+4*(theDiphotonPair.first)->covEtaPhi()*(theDiphotonPair.first)->covEtaPhi()))/((theDiphotonPair.first)->covEtaEta()+(theDiphotonPair.first)->covPhiPhi()+sqrt(((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())*((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())+4*(theDiphotonPair.first)->covEtaPhi()*(theDiphotonPair.first)->covEtaPhi()))); else Ppholead_lambdaRatio -> push_back(0);
		//if (pholead_cEE != 0) pholead_lamdbaDivCov = (pholead_cEE+pholead_cPP-sqrt((pholead_cEE-pholead_cPP)*(pholead_cEE-pholead_cPP)+4*pholead_cEP*pholead_cEP))/pholead_cEE; else pholead_lamdbaDivCov = 0;
		if ((theDiphotonPair.first)->covEtaEta() != 0)   Ppholead_lamdbaDivCov -> push_back(((theDiphotonPair.first)->covEtaEta()+(theDiphotonPair.first)->covPhiPhi()-sqrt(((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())*((theDiphotonPair.first)->covEtaEta()-(theDiphotonPair.first)->covPhiPhi())+4*(theDiphotonPair.first)->covEtaPhi()*(theDiphotonPair.first)->covEtaPhi()))/(theDiphotonPair.first)->covEtaEta()); else Ppholead_lamdbaDivCov ->push_back(0);

	
		Ppholead_secondMomentMaj -> push_back((theDiphotonPair.first)->secondMomentMaj());
		Ppholead_secondMomentMin -> push_back((theDiphotonPair.first)->secondMomentMin());
		Ppholead_secondMomentAlpha -> push_back((theDiphotonPair.first)->secondMomentAlpha());
		Ppholead_etawidth -> push_back((theDiphotonPair.first)->superCluster()->etaWidth());
		Ppholead_phiwidth -> push_back((theDiphotonPair.first)->superCluster()->phiWidth());
		Ppholead_sigieta -> push_back((theDiphotonPair.first)->sigmaIetaIeta());
		Ppholead_SCbr -> push_back((theDiphotonPair.first)->superCluster()->phiWidth()/(theDiphotonPair.first)->superCluster()->etaWidth());
                Ppholead_seedEnergy -> push_back((theDiphotonPair.first)->superCluster()->seedEnergy());
                Ppholead_seedTime -> push_back((theDiphotonPair.first)->superCluster()->seedTime());
                Ppholead_SCphi -> push_back((theDiphotonPair.first)->superCluster()->Phi());
                Ppholead_SCEtraw -> push_back((theDiphotonPair.first)->superCluster()->rawEnergy()* sin((theDiphotonPair.first)->superCluster()->Theta()));
                Ppholead_SCEt -> push_back((theDiphotonPair.first)->superCluster()->Pt());
                Ppholead_SCr9 -> push_back((theDiphotonPair.first)->r9()); 
                Ppholead_SCnbBC -> push_back((theDiphotonPair.first)->superCluster()->nBasicClusters());
                Ppholead_SCnXtal -> push_back((theDiphotonPair.first)->superCluster()->nXtals());
                
                if (doHLTobject){
                    Ppholead_isMatchingWithHLTObject -> push_back(findMatchingWithAnHLTObjet((theDiphotonPair.first), HLTObjects, theHTLobject));
                }

                if ( (((theDiphotonPair.first)->superCluster()->seedSeverity()==5)||((theDiphotonPair.first)->superCluster()->seedSeverity()==4)||(((theDiphotonPair.first)->superCluster()->seedRecoFlag()==2)&&((theDiphotonPair.first)->superCluster()->seedEnergy()<130))||(((theDiphotonPair.first)->superCluster()->seedEnergy()>=130)&&((theDiphotonPair.first)->superCluster()->seedTime()<0)&&((theDiphotonPair.first)->superCluster()->seedRecoFlag()==2))||((theDiphotonPair.first)->sigmaIetaIeta()<0.001)||(TMath::Sqrt((theDiphotonPair.first)->covPhiPhi())<0.001))&&(theDiphotonPair.first)->isEBPho()==1){   
                       Ppholead_isAspike -> push_back(1);}
                else{
                       Ppholead_isAspike -> push_back(0);
                }
                
                Ppholead_isConverted -> push_back(0);
                if ((theDiphotonPair.first)->convNTracks() > 0 ) Ppholead_isConverted -> push_back(1);
                Ppholead_NtrackConv -> push_back((theDiphotonPair.first)->convNTracks());
                Ppholead_convEoverP -> push_back((theDiphotonPair.first)->convEoverP());
                Ppholead_convMass -> push_back((theDiphotonPair.first)->convMass());  
                Ppholead_convCotanTheta -> push_back((theDiphotonPair.first)->convCotanTheta());
                Ppholead_convLikely -> push_back((theDiphotonPair.first)->convLikely());
                Ppholead_convVertexX -> push_back((theDiphotonPair.first)->convVertex().x());
                Ppholead_convVertexY -> push_back((theDiphotonPair.first)->convVertex().y());                
                Ppholead_convVertexZ -> push_back((theDiphotonPair.first)->convVertex().z());

                int lead_isAlsoRecoAsElectron=0;float lead_fBrem=-1;float lead_momentumCorrected = -1;float lead_d0 = -1;float lead_tightEleId = -1;float lead_eleTrkIso = -1;float lead_eleEcalIso = -1;float lead_eleHcalIso = -1;float lead_eleDeltaPhiIn = -1;float lead_eleDeltaEtaIn = -1;float lead_eleHoE = -1;float lead_eleSigmaIetaIeta = -1;int lead_eleMissHits = -1;float lead_eleDistConvPartner=-1;float lead_eleDcotConvPartner=-1;
                matchWithAnElectron((theDiphotonPair.first), electrons, &lead_isAlsoRecoAsElectron, &lead_fBrem, &lead_momentumCorrected, &lead_d0,&lead_tightEleId, &lead_eleTrkIso, &lead_eleEcalIso, &lead_eleHcalIso, &lead_eleDeltaPhiIn, &lead_eleDeltaEtaIn, &lead_eleHoE, &lead_eleSigmaIetaIeta, &lead_eleMissHits, &lead_eleDistConvPartner, &lead_eleDcotConvPartner); 
                Ppholead_fBrem -> push_back(lead_fBrem);
             
                if (doPhotonConversionMC){
                      int lead_MCisConverted=-1;float lead_MCconvEoverP=-1;float lead_MCconvMass=-1;float lead_MCconvCotanTheta=-1;float lead_MCconvVertexX=-1;float lead_MCconvVertexY=-1;float lead_MCconvVertexZ=-1;

                      findConversionMCtruth((theDiphotonPair.first), mcPhotons, lead_MCisConverted, lead_MCconvEoverP, lead_MCconvMass, lead_MCconvCotanTheta, lead_MCconvVertexX, lead_MCconvVertexY, lead_MCconvVertexZ);
                      Ppholead_MCisConverted -> push_back(lead_MCisConverted);
                      Ppholead_MCconvEoverP -> push_back(lead_MCconvEoverP);
                      Ppholead_MCconvMass -> push_back(lead_MCconvMass);      
                      Ppholead_convCotanTheta -> push_back(lead_MCconvCotanTheta);
                      Ppholead_MCconvVertexX -> push_back(lead_MCconvVertexX);
                      Ppholead_MCconvVertexY -> push_back(lead_MCconvVertexY); 
                      Ppholead_MCconvVertexZ -> push_back(lead_MCconvVertexZ);
                }


                Ppholead_xVertex -> push_back((theDiphotonPair.first)->vx());
                Ppholead_yVertex -> push_back((theDiphotonPair.first)->vy());
                Ppholead_zVertex -> push_back((theDiphotonPair.first)->vz());
		Ppholead_HcalIso -> push_back((theDiphotonPair.first)->dR04IsolationHcalRecHit());
		Ppholead_EcalIso -> push_back((theDiphotonPair.first)->dR04IsolationEcalRecHit());
		Ppholead_TrackerIso -> push_back((theDiphotonPair.first)->dR04IsolationHollowTrkCone());
		Ppholead_HcalIso_MIT03 -> push_back((theDiphotonPair.first)->dR03IsolationHcalRecHit()-0.005*Plead.Et());
		Ppholead_EcalIso_MIT03 -> push_back((theDiphotonPair.first)->dR03IsolationEcalRecHit()-0.012*Plead.Et());
		Ppholead_TrackerIso_MIT03 -> push_back((theDiphotonPair.first)->dR03IsolationHollowTrkCone()-0.002*Plead.Et());
		Ppholead_HcalEcal_MIT03 -> push_back((theDiphotonPair.first)->dR03IsolationEcalRecHit()+(theDiphotonPair.first)->dR03IsolationHcalRecHit()-theEvent->rho()*0.17);
		TRootVertex* theBestVertexLead= (TRootVertex*) vertices->At(0);
                Ppholead_AbsTrackerIso -> push_back(localtrackIsolation(*theBestVertexLead, tracks, *(theDiphotonPair.first), *beamSpot,0.3));
		Ppholead_HcalIsodR03 -> push_back((theDiphotonPair.first)->dR03IsolationHcalRecHit());
		Ppholead_EcalIsodR03 -> push_back((theDiphotonPair.first)->dR03IsolationEcalRecHit());
		Ppholead_TrackerIsodR03 -> push_back((theDiphotonPair.first)->dR03IsolationHollowTrkCone());
		Ppholead_hoe -> push_back((theDiphotonPair.first)->hoe());

                Ppholead_Cat-> push_back(findThePhoCat((theDiphotonPair.first)));
		Ppholead_HasPixSeed -> push_back((theDiphotonPair.first)->hasPixelSeed());
		Ppholead_seedSeverity -> push_back((theDiphotonPair.first)->superCluster()->seedSeverity()); 
		Ppholead_recoFlag -> push_back((theDiphotonPair.first)->superCluster()->seedRecoFlag());
		Ppholead_isEB -> push_back((theDiphotonPair.first)->isEBPho());
		Ppholead_isEE -> push_back((theDiphotonPair.first)->isEEPho());





              //now trail photon
                if (doMC) {
                    int trail_GenId=-1;int trail_MotherId=-1;int trail_isGenElectron=-1;int trail_isPromptGenPho=-1;int trail_isFromQuarkGen=-1;int trail_isPi0Gen=-1;int trail_isEtaGen=-1;int trail_isRhoGen=-1;int trail_isOmegaGen=-1;float trail_PromptGenIsoEnergyStatus1_cone02=0;float trail_PromptGenIsoEnergyStatus2_cone02=0;float trail_trueE=0;float trail_truePx=0;float trail_truePy=0;float trail_truePz=0;float trail_trueEta=0;float trail_truePhi=0;
                    
                    doGenInfo(theDiphotonPair.second, mcParticles, &(trail_GenId), &(trail_MotherId), &(trail_isGenElectron), &(trail_isPromptGenPho), &(trail_isFromQuarkGen), &(trail_isPi0Gen), &(trail_isEtaGen), &(trail_isRhoGen), &(trail_isOmegaGen), &(trail_PromptGenIsoEnergyStatus1_cone02), &(trail_PromptGenIsoEnergyStatus2_cone02),&(trail_trueE),&(trail_truePx),&(trail_truePy),&(trail_truePz),&(trail_trueEta),&(trail_truePhi), 0.2);
                    Pphotrail_GenId->push_back(trail_GenId); 
                    Pphotrail_MotherId->push_back(trail_MotherId); 
                    Pphotrail_isGenElectron->push_back(trail_isGenElectron); 
                    Pphotrail_isPromptGenPho->push_back(trail_isPromptGenPho); 
                    Pphotrail_isFromQuarkGen->push_back(trail_isFromQuarkGen); 
                    Pphotrail_isPi0Gen->push_back(trail_isPi0Gen); 
                    Pphotrail_isEtaGen->push_back(trail_isEtaGen); 
                    Pphotrail_isRhoGen->push_back(trail_isRhoGen); 
                    Pphotrail_isOmegaGen->push_back(trail_isOmegaGen); 
                    Pphotrail_PromptGenIsoEnergyStatus1_cone02->push_back(trail_PromptGenIsoEnergyStatus1_cone02); 
                    Pphotrail_PromptGenIsoEnergyStatus2_cone02->push_back(trail_PromptGenIsoEnergyStatus2_cone02); 
                    Pphotrail_trueE->push_back(trail_trueE); 
                    Pphotrail_truePx->push_back(trail_truePx); 
                    Pphotrail_truePy->push_back(trail_truePy); 
                    Pphotrail_truePz->push_back(trail_truePz); 
                    Pphotrail_trueEta->push_back(trail_trueEta); 
                    Pphotrail_truePhi->push_back(trail_truePhi); 
                }


                Pphotrail_event_processId -> push_back(theEvent->processID());
                Pphotrail_et -> push_back((theDiphotonPair.second)->Et());
                Pphotrail_eta -> push_back((theDiphotonPair.second)->Eta());
                Pphotrail_SCeta -> push_back((theDiphotonPair.second)->superCluster()->Eta()); 
		Pphotrail_r9 -> push_back((theDiphotonPair.second)->r9());
		Pphotrail_cPP -> push_back((theDiphotonPair.second)->covPhiPhi());
		Pphotrail_cEE -> push_back((theDiphotonPair.second)->covEtaEta());
		Pphotrail_cEP -> push_back((theDiphotonPair.second)->covEtaPhi());

		Pphotrail_r19 -> push_back((theDiphotonPair.second)->r19());
		Pphotrail_SCEraw -> push_back((theDiphotonPair.second)->superCluster()->rawEnergy());
		Pphotrail_eMax -> push_back((theDiphotonPair.second)->eMax());
		Pphotrail_e2x2 -> push_back((theDiphotonPair.second)->e2x2());
		Pphotrail_e5x5 -> push_back((theDiphotonPair.second)->e5x5());

		if ((theDiphotonPair.second)->superCluster()->rawEnergy()!=0)  Pphotrail_ratioSeed ->push_back((theDiphotonPair.second)->eMax()/(theDiphotonPair.second)->superCluster()->rawEnergy()); else Pphotrail_ratioSeed ->push_back(0);

		if ((theDiphotonPair.second)->e5x5()!=0) Pphotrail_ratioS4 -> push_back( (theDiphotonPair.second)->e2x2() / (theDiphotonPair.second)->e5x5()); else Pphotrail_ratioS4 -> push_back(0);

                if ( ((theDiphotonPair.second)->covEtaEta()+(theDiphotonPair.second)->covPhiPhi()+sqrt(((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())*((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())+4*(theDiphotonPair.second)->covEtaPhi()*(theDiphotonPair.second)->covEtaPhi()))!=0)  Pphotrail_lambdaRatio -> push_back(((theDiphotonPair.second)->covEtaEta()+(theDiphotonPair.second)->covPhiPhi()-sqrt(((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())*((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())+4*(theDiphotonPair.second)->covEtaPhi()*(theDiphotonPair.second)->covEtaPhi()))/((theDiphotonPair.second)->covEtaEta()+(theDiphotonPair.second)->covPhiPhi()+sqrt(((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())*((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())+4*(theDiphotonPair.second)->covEtaPhi()*(theDiphotonPair.second)->covEtaPhi()))); else Pphotrail_lambdaRatio -> push_back(0);

		if ((theDiphotonPair.second)->covEtaEta() != 0) Pphotrail_lamdbaDivCov -> push_back(((theDiphotonPair.second)->covEtaEta()+(theDiphotonPair.second)->covPhiPhi()-sqrt(((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())*((theDiphotonPair.second)->covEtaEta()-(theDiphotonPair.second)->covPhiPhi())+4*(theDiphotonPair.second)->covEtaPhi()*(theDiphotonPair.second)->covEtaPhi()))/(theDiphotonPair.second)->covEtaEta()); else Pphotrail_lamdbaDivCov -> push_back(0);
	
		Pphotrail_secondMomentMaj -> push_back((theDiphotonPair.second)->secondMomentMaj());
		Pphotrail_secondMomentMin -> push_back((theDiphotonPair.second)->secondMomentMin());
		Pphotrail_secondMomentAlpha -> push_back((theDiphotonPair.second)->secondMomentAlpha());
		Pphotrail_etawidth -> push_back((theDiphotonPair.second)->superCluster()->etaWidth());
		Pphotrail_phiwidth -> push_back((theDiphotonPair.second)->superCluster()->phiWidth());
		Pphotrail_sigieta -> push_back((theDiphotonPair.second)->sigmaIetaIeta());
		Pphotrail_SCbr -> push_back((theDiphotonPair.second)->superCluster()->phiWidth()/(theDiphotonPair.second)->superCluster()->etaWidth());
                Pphotrail_seedEnergy -> push_back((theDiphotonPair.second)->superCluster()->seedEnergy());
                Pphotrail_seedTime -> push_back((theDiphotonPair.second)->superCluster()->seedTime());
                Pphotrail_SCphi -> push_back((theDiphotonPair.second)->superCluster()->Phi());
                Pphotrail_SCEtraw -> push_back((theDiphotonPair.second)->superCluster()->rawEnergy()* sin((theDiphotonPair.second)->superCluster()->Theta()));
                Pphotrail_SCEt -> push_back((theDiphotonPair.second)->superCluster()->Pt());
                Pphotrail_SCr9 -> push_back((theDiphotonPair.second)->r9()); 
                Pphotrail_SCnbBC -> push_back((theDiphotonPair.second)->superCluster()->nBasicClusters());
                Pphotrail_SCnXtal -> push_back((theDiphotonPair.second)->superCluster()->nXtals());

                if (doHLTobject){
                    Pphotrail_isMatchingWithHLTObject -> push_back(findMatchingWithAnHLTObjet((theDiphotonPair.second), HLTObjects, theHTLobject));
                }

                if ( (((theDiphotonPair.second)->superCluster()->seedSeverity()==5)||((theDiphotonPair.second)->superCluster()->seedSeverity()==4)||(((theDiphotonPair.second)->superCluster()->seedRecoFlag()==2)&&((theDiphotonPair.second)->superCluster()->seedEnergy()<130))||(((theDiphotonPair.second)->superCluster()->seedEnergy()>=130)&&((theDiphotonPair.second)->superCluster()->seedTime()<0)&&((theDiphotonPair.second)->superCluster()->seedRecoFlag()==2))||((theDiphotonPair.second)->sigmaIetaIeta()<0.001)||(TMath::Sqrt((theDiphotonPair.second)->covPhiPhi())<0.001))&&(theDiphotonPair.second)->isEBPho()==1){   
                       Pphotrail_isAspike -> push_back(1);}
                else{
                       Pphotrail_isAspike -> push_back(0);
                }

                Pphotrail_isConverted -> push_back(0);
                if ((theDiphotonPair.second)->convNTracks() > 0 ) Pphotrail_isConverted -> push_back(1);
                Pphotrail_NtrackConv -> push_back((theDiphotonPair.second)->convNTracks());
                Pphotrail_convEoverP -> push_back((theDiphotonPair.second)->convEoverP());
                Pphotrail_convMass -> push_back((theDiphotonPair.second)->convMass());  
                Pphotrail_convCotanTheta -> push_back((theDiphotonPair.second)->convCotanTheta());
                Pphotrail_convLikely -> push_back((theDiphotonPair.second)->convLikely());
                Pphotrail_convVertexX -> push_back((theDiphotonPair.second)->convVertex().x());
                Pphotrail_convVertexY -> push_back((theDiphotonPair.second)->convVertex().y());                
                Pphotrail_convVertexZ -> push_back((theDiphotonPair.second)->convVertex().z());

                int trail_isAlsoRecoAsElectron=0;float trail_fBrem=-1;float trail_momentumCorrected = -1;float trail_d0 = -1;float trail_tightEleId = -1;float trail_eleTrkIso = -1;float trail_eleEcalIso = -1;float trail_eleHcalIso = -1;float trail_eleDeltaPhiIn = -1;float trail_eleDeltaEtaIn = -1;float trail_eleHoE = -1;float trail_eleSigmaIetaIeta = -1;int trail_eleMissHits = -1;float trail_eleDistConvPartner=-1;float trail_eleDcotConvPartner=-1;
                matchWithAnElectron((theDiphotonPair.second), electrons, &trail_isAlsoRecoAsElectron, &trail_fBrem, &trail_momentumCorrected, &trail_d0,&trail_tightEleId, &trail_eleTrkIso, &trail_eleEcalIso, &trail_eleHcalIso, &trail_eleDeltaPhiIn, &trail_eleDeltaEtaIn, &trail_eleHoE, &trail_eleSigmaIetaIeta, &trail_eleMissHits, &trail_eleDistConvPartner, &trail_eleDcotConvPartner); 
                Pphotrail_fBrem -> push_back(trail_fBrem);

                if (doPhotonConversionMC){
                      int trail_MCisConverted=-1;float trail_MCconvEoverP=-1;float trail_MCconvMass=-1;float trail_MCconvCotanTheta=-1;float trail_MCconvVertexX=-1;float trail_MCconvVertexY=-1;float trail_MCconvVertexZ=-1;

                      findConversionMCtruth((theDiphotonPair.second), mcPhotons, trail_MCisConverted, trail_MCconvEoverP, trail_MCconvMass, trail_MCconvCotanTheta, trail_MCconvVertexX, trail_MCconvVertexY, trail_MCconvVertexZ);
                      Pphotrail_MCisConverted -> push_back(trail_MCisConverted);
                      Pphotrail_MCconvEoverP -> push_back(trail_MCconvEoverP);
                      Pphotrail_MCconvMass -> push_back(trail_MCconvMass);      
                      Pphotrail_convCotanTheta -> push_back(trail_MCconvCotanTheta);
                      Pphotrail_MCconvVertexX -> push_back(trail_MCconvVertexX);
                      Pphotrail_MCconvVertexY -> push_back(trail_MCconvVertexY); 
                      Pphotrail_MCconvVertexZ -> push_back(trail_MCconvVertexZ);
                }


                Pphotrail_xVertex -> push_back((theDiphotonPair.second)->vx());
                Pphotrail_yVertex -> push_back((theDiphotonPair.second)->vy());
                Pphotrail_zVertex -> push_back((theDiphotonPair.second)->vz());
		Pphotrail_HcalIso -> push_back((theDiphotonPair.second)->dR04IsolationHcalRecHit());
		Pphotrail_EcalIso -> push_back((theDiphotonPair.second)->dR04IsolationEcalRecHit());
		Pphotrail_TrackerIso -> push_back((theDiphotonPair.second)->dR04IsolationHollowTrkCone());
		Pphotrail_HcalIso_MIT03 -> push_back((theDiphotonPair.second)->dR03IsolationHcalRecHit()-0.005*Ptrail.Et());
		Pphotrail_EcalIso_MIT03 -> push_back((theDiphotonPair.second)->dR03IsolationEcalRecHit()-0.012*Ptrail.Et());
		Pphotrail_TrackerIso_MIT03 -> push_back((theDiphotonPair.second)->dR03IsolationHollowTrkCone()-0.002*Ptrail.Et());
		Pphotrail_HcalEcal_MIT03 -> push_back((theDiphotonPair.second)->dR03IsolationEcalRecHit()+(theDiphotonPair.second)->dR03IsolationHcalRecHit()-theEvent->rho()*0.17);
		TRootVertex* theBestVertexTrail= (TRootVertex*) vertices->At(0);
                Pphotrail_AbsTrackerIso -> push_back(localtrackIsolation(*theBestVertexTrail, tracks, *(theDiphotonPair.second), *beamSpot,0.3));
		Pphotrail_HcalIsodR03 -> push_back((theDiphotonPair.second)->dR03IsolationHcalRecHit());
		Pphotrail_EcalIsodR03 -> push_back((theDiphotonPair.second)->dR03IsolationEcalRecHit());
		Pphotrail_TrackerIsodR03 -> push_back((theDiphotonPair.second)->dR03IsolationHollowTrkCone());
		Pphotrail_hoe -> push_back((theDiphotonPair.second)->hoe());

                Pphotrail_Cat-> push_back(findThePhoCat((theDiphotonPair.second)));
		Pphotrail_HasPixSeed -> push_back((theDiphotonPair.second)->hasPixelSeed());
		Pphotrail_seedSeverity -> push_back((theDiphotonPair.second)->superCluster()->seedSeverity()); 
		Pphotrail_recoFlag -> push_back((theDiphotonPair.second)->superCluster()->seedRecoFlag());
		Pphotrail_isEB -> push_back((theDiphotonPair.second)->isEBPho());
		Pphotrail_isEE -> push_back((theDiphotonPair.second)->isEEPho());


		//myTree_->Fill();
    }


void endMacro(){
//	myTree_->Write();
	myFile->Write();
	myFile->Close();
}


//miniTreeMaker(){
int main(){
	cout << "coucou" << endl;
	gSystem->Load("/sps/cms/jfan/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/src/libToto.so");
	inputEventTree->Add("dcap://ccdcapcms.in2p3.fr:22125//pnfs/in2p3.fr/data/cms/t2data//store/user/jfan/MCproduction423Fall11/GJetPt20/MC_EG_goodVtx_noscrapping_62_2_Ufn.root");


	myFile=new TFile("theMiniTree.root","RECREATE");


	beginMacro();
   
       
	int NbEvents = inputEventTree->GetEntries();	cout << "NbEvents = " << NbEvents << endl;
	//NbEvents = 10000;
        int NgoodEvents = 0;
        int NgoodPairs = 0;
        //int ngoodPairs = 0;
	int NbHLT20 = 0;
	int NbPhotons = 0;
	int NbPhotonsCleaned = 0;
	int NbPhotonsAccept = 0;
	int NbPhotonsCutEt = 0;
	int NbPhotonsPb = 0;

	int phoPair = 0;
	int phoPairLeadCiC0 = 0;	
	int phoPairLeadCiC1 = 0;	
	int phoPairLeadCiC2 = 0;	
	int phoPairLeadCiC3 = 0;	
	int phoPairLeadCiC4 = 0;	
	int phoPairLeadCiC5 = 0;	
	int phoPairLeadCiC6 = 0;	
	int phoPairLeadCiC7 = 0;	
	int phoPairTrailCiC0 = 0;	
	int phoPairTrailCiC1 = 0;	
	int phoPairTrailCiC2 = 0;	
	int phoPairTrailCiC3 = 0;	
	int phoPairTrailCiC4 = 0;	
	int phoPairTrailCiC5 = 0;	
	int phoPairTrailCiC6 = 0;	
	int phoPairTrailCiC7 = 0;
	int phoAfterCutEt = 0;
	int phoGoodPair = 0;	



	cout << "nomFichier=" << inputEventTree->GetFile()->GetName() << " nbEntries=" << NbEvents << endl;
	for (int ievt  = 0 ; ievt < NbEvents ; ievt++){

            inputEventTree->GetEvent(ievt);


            int NbPhoton = photons->GetEntriesFast();
            if (NbPhoton < 2) continue;
            ngoodPairs = accumulation(NbPhoton);
            cout << "nb of photons " << NbPhoton << "       nb of photon pairs " << ngoodPairs <<endl;
            NgoodPairs = NgoodPairs + ngoodPairs;




 
            Pdipho_mgg->clear();
            Pdipho_mgg->reserve(ngoodPairs);
            Pdipho_qt->clear();
            Pdipho_qt-> reserve(ngoodPairs);
            Pdipho_ql->clear();
            Pdipho_ql-> reserve(ngoodPairs);
            Pdipho_deltaR->clear();
            Pdipho_deltaR-> reserve(ngoodPairs);
            Pdipho_costhetastar->clear();
            Pdipho_costhetastar-> reserve(ngoodPairs);
            Pdipho_eta->clear();
            Pdipho_eta-> reserve(ngoodPairs);
            Pdipho_etastar->clear();
            Pdipho_etastar-> reserve(ngoodPairs);
            PdiphoMC_mgg->clear();
            PdiphoMC_mgg-> reserve(ngoodPairs);
            PdiphoMC_qt->clear();
            PdiphoMC_qt-> reserve(ngoodPairs);
            PdiphoMC_ql->clear();
            PdiphoMC_ql-> reserve(ngoodPairs);
            PdiphoMC_deltaR->clear();
            PdiphoMC_deltaR-> reserve(ngoodPairs);
            PdiphoMC_costhetastar->clear();
            PdiphoMC_costhetastar-> reserve(ngoodPairs);
            PdiphoMC_eta->clear();
            PdiphoMC_eta-> reserve(ngoodPairs);
            PdiphoMC_etastar->clear();
            PdiphoMC_etastar-> reserve(ngoodPairs);
      
            Ppholead_event_processId->clear(); 
            Ppholead_event_processId-> reserve(ngoodPairs); 
            Ppholead_isMatchingWithMC->clear(); 
            Ppholead_isMatchingWithMC-> reserve(ngoodPairs); 
            Ppholead_GenId->clear();
            Ppholead_GenId-> reserve(ngoodPairs);            
            Ppholead_MotherId->clear();
            Ppholead_MotherId-> reserve(ngoodPairs);            
            Ppholead_isPromptGenPho->clear();
            Ppholead_isPromptGenPho-> reserve(ngoodPairs);            
            Ppholead_isFromQuarkGen->clear();
            Ppholead_isFromQuarkGen-> reserve(ngoodPairs);            
            Ppholead_isPi0Gen->clear();
            Ppholead_isPi0Gen-> reserve(ngoodPairs);            
            Ppholead_isEtaGen->clear();
            Ppholead_isEtaGen-> reserve(ngoodPairs);            
            Ppholead_isRhoGen->clear();
            Ppholead_isRhoGen-> reserve(ngoodPairs);            
            Ppholead_isOmegaGen->clear();
            Ppholead_isOmegaGen-> reserve(ngoodPairs);            
            Ppholead_isGenElectron->clear();
            Ppholead_isGenElectron-> reserve(ngoodPairs);            
            Ppholead_eventPassHLT_Photon10_L1R->clear();
            Ppholead_eventPassHLT_Photon10_L1R-> reserve(ngoodPairs);            
            Ppholead_eventPassHLT_Photon15_L1R->clear();
            Ppholead_eventPassHLT_Photon15_L1R-> reserve(ngoodPairs);            
            Ppholead_eventPassHLT_DoublePhoton10_L1R->clear();
            Ppholead_eventPassHLT_DoublePhoton10_L1R-> reserve(ngoodPairs);            
            Ppholead_PromptGenIsoEnergyStatus1_cone02->clear();
            Ppholead_PromptGenIsoEnergyStatus1_cone02-> reserve(ngoodPairs);            
            Ppholead_PromptGenIsoEnergyStatus2_cone02->clear();
            Ppholead_PromptGenIsoEnergyStatus2_cone02-> reserve(ngoodPairs);            
            Ppholead_PromptGenIsoEnergyStatus1_cone03->clear();
            Ppholead_PromptGenIsoEnergyStatus1_cone03-> reserve(ngoodPairs);            
            Ppholead_PromptGenIsoEnergyStatus2_cone03->clear();
            Ppholead_PromptGenIsoEnergyStatus2_cone03-> reserve(ngoodPairs);            
            Ppholead_PromptGenIsoEnergyStatus1_cone035->clear();
            Ppholead_PromptGenIsoEnergyStatus2_cone035-> reserve(ngoodPairs);            
            Ppholead_PromptGenIsoEnergyStatus1_cone04->clear();
            Ppholead_PromptGenIsoEnergyStatus1_cone04-> reserve(ngoodPairs);            
            Ppholead_PromptGenIsoEnergyStatus2_cone04->clear();
            Ppholead_PromptGenIsoEnergyStatus2_cone04-> reserve(ngoodPairs);            
            Ppholead_trueE->clear();
            Ppholead_trueE-> reserve(ngoodPairs);            
            Ppholead_truePx->clear();
            Ppholead_truePx-> reserve(ngoodPairs);            
            Ppholead_truePy->clear();
            Ppholead_truePy-> reserve(ngoodPairs);            
            Ppholead_truePz->clear();
            Ppholead_truePz-> reserve(ngoodPairs);            
            Ppholead_trueEta->clear();
            Ppholead_trueEta-> reserve(ngoodPairs);            
            Ppholead_truePhi->clear();
            Ppholead_truePhi-> reserve(ngoodPairs);            
            Ppholead_MCisConverted->clear();
            Ppholead_MCisConverted-> reserve(ngoodPairs);            
            Ppholead_MCconvEoverP->clear();
            Ppholead_MCconvEoverP-> reserve(ngoodPairs);            
            Ppholead_MCconvCotanTheta->clear();
            Ppholead_MCconvCotanTheta-> reserve(ngoodPairs);            
            Ppholead_MCconvVertexX->clear();
            Ppholead_MCconvVertexX-> reserve(ngoodPairs);            
            Ppholead_MCconvVertexY->clear();
            Ppholead_MCconvVertexY-> reserve(ngoodPairs);            
            Ppholead_MCconvVertexZ->clear();
            Ppholead_MCconvVertexZ-> reserve(ngoodPairs);            
            Ppholead_eleMCtruthBrem->clear();
            Ppholead_eleMCtruthBrem-> reserve(ngoodPairs);            
            Ppholead_eleMCtruthNBrem->clear();
            Ppholead_eleMCtruthNBrem-> reserve(ngoodPairs);            

            Ppholead_et->clear(); 
            Ppholead_et-> reserve(ngoodPairs); 
            Ppholead_eta->clear();
            Ppholead_eta-> reserve(ngoodPairs);
            Ppholead_SCeta->clear(); 
            Ppholead_SCeta-> reserve(ngoodPairs); 
            Ppholead_r9->clear(); 
            Ppholead_r9-> reserve(ngoodPairs); 
            Ppholead_cPP->clear(); 
            Ppholead_cPP-> reserve(ngoodPairs); 
            Ppholead_cEP->clear();
            Ppholead_cEP-> reserve(ngoodPairs);
            Ppholead_cEE->clear(); 
            Ppholead_cEE-> reserve(ngoodPairs); 
            Ppholead_r19->clear(); 
            Ppholead_r19-> reserve(ngoodPairs); 
            Ppholead_SCEraw->clear(); 
            Ppholead_SCEraw-> reserve(ngoodPairs); 
            Ppholead_eMax->clear(); 
            Ppholead_eMax-> reserve(ngoodPairs); 
            Ppholead_e2x2->clear(); 
            Ppholead_e2x2-> reserve(ngoodPairs); 
            Ppholead_e5x5->clear(); 
            Ppholead_e5x5-> reserve(ngoodPairs); 
            Ppholead_ratioSeed->clear(); 
            Ppholead_ratioSeed-> reserve(ngoodPairs); 
            Ppholead_ratioS4->clear(); 
            Ppholead_ratioS4-> reserve(ngoodPairs); 
            Ppholead_lambdaRatio->clear(); 
            Ppholead_lambdaRatio-> reserve(ngoodPairs); 
            Ppholead_lamdbaDivCov->clear(); 
            Ppholead_lamdbaDivCov-> reserve(ngoodPairs); 
            Ppholead_secondMomentMaj->clear(); 
            Ppholead_secondMomentMaj-> reserve(ngoodPairs); 
            Ppholead_secondMomentMin->clear(); 
            Ppholead_secondMomentMin-> reserve(ngoodPairs); 
            Ppholead_secondMomentAlpha->clear(); 
            Ppholead_secondMomentAlpha-> reserve(ngoodPairs); 
            Ppholead_covAngle->clear(); 
            Ppholead_covAngle-> reserve(ngoodPairs); 
            Ppholead_covAngle2->clear(); 
            Ppholead_covAngle2-> reserve(ngoodPairs); 
            Ppholead_S9overS9minusS1S2->clear(); 
            Ppholead_S9overS9minusS1S2-> reserve(ngoodPairs); 
            Ppholead_etawidth->clear(); 
            Ppholead_etawidth-> reserve(ngoodPairs); 
            Ppholead_phiwidth->clear(); 
            Ppholead_phiwidth-> reserve(ngoodPairs); 
            Ppholead_sigieta->clear(); 
            Ppholead_sigieta-> reserve(ngoodPairs); 
            Ppholead_SCbr->clear(); 
            Ppholead_SCbr-> reserve(ngoodPairs); 
            Ppholead_seedEnergy->clear();
            Ppholead_seedEnergy-> reserve(ngoodPairs);
            Ppholead_seedTime->clear();
            Ppholead_seedTime-> reserve(ngoodPairs);            
            Ppholead_SCphi->clear();
            Ppholead_SCphi-> reserve(ngoodPairs);            
            Ppholead_SCEtraw->clear();
            Ppholead_SCEtraw-> reserve(ngoodPairs);            
            Ppholead_SCEt->clear();
            Ppholead_SCEt-> reserve(ngoodPairs);            
            Ppholead_SCr9->clear();
            Ppholead_SCr9-> reserve(ngoodPairs);            
            Ppholead_SCnbBC->clear();
            Ppholead_SCnbBC-> reserve(ngoodPairs);            
            Ppholead_SCnXtal->clear();
            Ppholead_SCnXtal-> reserve(ngoodPairs);            
            Ppholead_isMatchingWithHLTObject->clear();
            Ppholead_isMatchingWithHLTObject-> reserve(ngoodPairs);            
            Ppholead_isConverted->clear();
            Ppholead_isConverted-> reserve(ngoodPairs);            
            Ppholead_NtrackConv->clear();
            Ppholead_NtrackConv-> reserve(ngoodPairs);            
            Ppholead_convEoverP->clear();
            Ppholead_convEoverP-> reserve(ngoodPairs);            
            Ppholead_convMass->clear();
            Ppholead_convMass-> reserve(ngoodPairs);            
            Ppholead_convCotanTheta->clear();
            Ppholead_convCotanTheta-> reserve(ngoodPairs);            
            Ppholead_MCconvMass->clear();
            Ppholead_MCconvMass-> reserve(ngoodPairs);            
            Ppholead_convLikely->clear();
            Ppholead_convLikely-> reserve(ngoodPairs);            
            Ppholead_convVertexX->clear();
            Ppholead_convVertexX-> reserve(ngoodPairs);            
            Ppholead_convVertexY->clear();
            Ppholead_convVertexY-> reserve(ngoodPairs);            
            Ppholead_convVertexZ->clear();
            Ppholead_convVertexZ-> reserve(ngoodPairs);
            Ppholead_fBrem->clear();
            Ppholead_fBrem-> reserve(ngoodPairs);            
            Ppholead_isAspike->clear();
            Ppholead_isAspike-> reserve(ngoodPairs);            
            Ppholead_xVertex->clear();
            Ppholead_xVertex-> reserve(ngoodPairs);            
            Ppholead_yVertex->clear();
            Ppholead_yVertex-> reserve(ngoodPairs);            
            Ppholead_zVertex->clear();
            Ppholead_zVertex-> reserve(ngoodPairs);            

            Ppholead_HcalIso->clear(); 
            Ppholead_HcalIso-> reserve(ngoodPairs); 
            Ppholead_EcalIso->clear(); 
            Ppholead_EcalIso-> reserve(ngoodPairs); 
            Ppholead_TrackerIso->clear(); 
            Ppholead_TrackerIso-> reserve(ngoodPairs); 
            Ppholead_HcalIso_MIT03->clear(); 
            Ppholead_HcalIso_MIT03-> reserve(ngoodPairs); 
            Ppholead_EcalIso_MIT03->clear(); 
            Ppholead_EcalIso_MIT03-> reserve(ngoodPairs); 
            Ppholead_TrackerIso_MIT03->clear(); 
            Ppholead_TrackerIso_MIT03-> reserve(ngoodPairs); 
            Ppholead_HcalEcal_MIT03->clear(); 
            Ppholead_HcalEcal_MIT03-> reserve(ngoodPairs); 
            Ppholead_AbsTrackerIso->clear(); 
            Ppholead_AbsTrackerIso-> reserve(ngoodPairs); 
            Ppholead_HcalIsodR03->clear(); 
            Ppholead_HcalIsodR03-> reserve(ngoodPairs); 
            Ppholead_EcalIsodR03->clear(); 
            Ppholead_EcalIsodR03-> reserve(ngoodPairs); 
            Ppholead_TrackerIsodR03->clear(); 
            Ppholead_TrackerIsodR03-> reserve(ngoodPairs); 
            Ppholead_HcalIsoPerso->clear(); 
            Ppholead_HcalIsoPerso-> reserve(ngoodPairs); 
            Ppholead_EcalIsoPerso->clear(); 
            Ppholead_EcalIsoPerso-> reserve(ngoodPairs); 
            Ppholead_hoe->clear(); 
            Ppholead_hoe-> reserve(ngoodPairs); 

            Ppholead_Cat->clear(); 
            Ppholead_Cat -> reserve(ngoodPairs);
            Ppholead_HasPixSeed->clear(); 
            Ppholead_HasPixSeed -> reserve(ngoodPairs);
            Ppholead_seedSeverity->clear(); 
            Ppholead_seedSeverity-> reserve(ngoodPairs); 
            Ppholead_recoFlag->clear(); 
            Ppholead_recoFlag-> reserve(ngoodPairs); 
            Ppholead_isEB->clear(); 
            Ppholead_isEB-> reserve(ngoodPairs); 
            Ppholead_isEE->clear(); 
            Ppholead_isEE-> reserve(ngoodPairs); 
           



            Pphotrail_event_processId->clear(); 
            Pphotrail_event_processId-> reserve(ngoodPairs); 
            Pphotrail_isMatchingWithMC->clear(); 
            Pphotrail_isMatchingWithMC-> reserve(ngoodPairs); 
            Pphotrail_GenId->clear();
            Pphotrail_GenId-> reserve(ngoodPairs);            
            Pphotrail_MotherId->clear();
            Pphotrail_MotherId-> reserve(ngoodPairs);            
            Pphotrail_isPromptGenPho->clear();
            Pphotrail_isPromptGenPho-> reserve(ngoodPairs);            
            Pphotrail_isFromQuarkGen->clear();
            Pphotrail_isFromQuarkGen-> reserve(ngoodPairs);            
            Pphotrail_isPi0Gen->clear();
            Pphotrail_isPi0Gen-> reserve(ngoodPairs);            
            Pphotrail_isEtaGen->clear();
            Pphotrail_isEtaGen-> reserve(ngoodPairs);            
            Pphotrail_isRhoGen->clear();
            Pphotrail_isRhoGen-> reserve(ngoodPairs);            
            Pphotrail_isOmegaGen->clear();
            Pphotrail_isOmegaGen-> reserve(ngoodPairs);            
            Pphotrail_isGenElectron->clear();
            Pphotrail_isGenElectron-> reserve(ngoodPairs);            
            Pphotrail_eventPassHLT_Photon10_L1R->clear();
            Pphotrail_eventPassHLT_Photon10_L1R-> reserve(ngoodPairs);            
            Pphotrail_eventPassHLT_Photon15_L1R->clear();
            Pphotrail_eventPassHLT_Photon15_L1R-> reserve(ngoodPairs);            
            Pphotrail_eventPassHLT_DoublePhoton10_L1R->clear();
            Pphotrail_eventPassHLT_DoublePhoton10_L1R-> reserve(ngoodPairs);            
            Pphotrail_PromptGenIsoEnergyStatus1_cone02->clear();
            Pphotrail_PromptGenIsoEnergyStatus1_cone02-> reserve(ngoodPairs);            
            Pphotrail_PromptGenIsoEnergyStatus2_cone02->clear();
            Pphotrail_PromptGenIsoEnergyStatus2_cone02-> reserve(ngoodPairs);            
            Pphotrail_PromptGenIsoEnergyStatus1_cone03->clear();
            Pphotrail_PromptGenIsoEnergyStatus1_cone03-> reserve(ngoodPairs);            
            Pphotrail_PromptGenIsoEnergyStatus2_cone03->clear();
            Pphotrail_PromptGenIsoEnergyStatus2_cone03-> reserve(ngoodPairs);            
            Pphotrail_PromptGenIsoEnergyStatus1_cone035->clear();
            Pphotrail_PromptGenIsoEnergyStatus2_cone035-> reserve(ngoodPairs);            
            Pphotrail_PromptGenIsoEnergyStatus1_cone04->clear();
            Pphotrail_PromptGenIsoEnergyStatus1_cone04-> reserve(ngoodPairs);            
            Pphotrail_PromptGenIsoEnergyStatus2_cone04->clear();
            Pphotrail_PromptGenIsoEnergyStatus2_cone04-> reserve(ngoodPairs);            
            Pphotrail_trueE->clear();
            Pphotrail_trueE-> reserve(ngoodPairs);            
            Pphotrail_truePx->clear();
            Pphotrail_truePx-> reserve(ngoodPairs);            
            Pphotrail_truePy->clear();
            Pphotrail_truePy-> reserve(ngoodPairs);            
            Pphotrail_truePz->clear();
            Pphotrail_truePz-> reserve(ngoodPairs);            
            Pphotrail_trueEta->clear();
            Pphotrail_trueEta-> reserve(ngoodPairs);            
            Pphotrail_truePhi->clear();
            Pphotrail_truePhi-> reserve(ngoodPairs);            
            Pphotrail_MCisConverted->clear();
            Pphotrail_MCisConverted-> reserve(ngoodPairs);            
            Pphotrail_MCconvEoverP->clear();
            Pphotrail_MCconvEoverP-> reserve(ngoodPairs);            
            Pphotrail_MCconvCotanTheta->clear();
            Pphotrail_MCconvCotanTheta-> reserve(ngoodPairs);            
            Pphotrail_MCconvMass->clear();
            Pphotrail_MCconvMass-> reserve(ngoodPairs);            
            Pphotrail_MCconvVertexX->clear();
            Pphotrail_MCconvVertexX-> reserve(ngoodPairs);            
            Pphotrail_MCconvVertexY->clear();
            Pphotrail_MCconvVertexY-> reserve(ngoodPairs);            
            Pphotrail_MCconvVertexZ->clear();
            Pphotrail_MCconvVertexZ-> reserve(ngoodPairs);            
            Pphotrail_eleMCtruthBrem->clear();
            Pphotrail_eleMCtruthBrem-> reserve(ngoodPairs);            
            Pphotrail_eleMCtruthNBrem->clear();

            Pphotrail_et->clear(); 
            Pphotrail_et-> reserve(ngoodPairs); 
            Pphotrail_eta->clear();
            Pphotrail_eta-> reserve(ngoodPairs);
            Pphotrail_SCeta->clear(); 
            Pphotrail_SCeta-> reserve(ngoodPairs); 
            Pphotrail_r9->clear(); 
            Pphotrail_r9-> reserve(ngoodPairs); 
            Pphotrail_cPP->clear(); 
            Pphotrail_cPP-> reserve(ngoodPairs); 
            Pphotrail_cEP->clear();
            Pphotrail_cEP-> reserve(ngoodPairs);
            Pphotrail_cEE->clear(); 
            Pphotrail_cEE-> reserve(ngoodPairs); 
            Pphotrail_r19->clear(); 
            Pphotrail_r19-> reserve(ngoodPairs); 
            Pphotrail_SCEraw->clear(); 
            Pphotrail_SCEraw-> reserve(ngoodPairs); 
            Pphotrail_eMax->clear(); 
            Pphotrail_eMax-> reserve(ngoodPairs); 
            Pphotrail_e2x2->clear(); 
            Pphotrail_e2x2-> reserve(ngoodPairs); 
            Pphotrail_e5x5->clear(); 
            Pphotrail_e5x5-> reserve(ngoodPairs); 
            Pphotrail_ratioSeed->clear(); 
            Pphotrail_ratioSeed-> reserve(ngoodPairs); 
            Pphotrail_ratioS4->clear(); 
            Pphotrail_ratioS4-> reserve(ngoodPairs); 
            Pphotrail_lambdaRatio->clear(); 
            Pphotrail_lambdaRatio-> reserve(ngoodPairs); 
            Pphotrail_lamdbaDivCov->clear(); 
            Pphotrail_lamdbaDivCov-> reserve(ngoodPairs); 
            Pphotrail_secondMomentMaj->clear(); 
            Pphotrail_secondMomentMaj-> reserve(ngoodPairs); 
            Pphotrail_secondMomentMin->clear(); 
            Pphotrail_secondMomentMin-> reserve(ngoodPairs); 
            Pphotrail_secondMomentAlpha->clear(); 
            Pphotrail_secondMomentAlpha-> reserve(ngoodPairs); 
            Pphotrail_covAngle->clear(); 
            Pphotrail_covAngle-> reserve(ngoodPairs); 
            Pphotrail_covAngle2->clear(); 
            Pphotrail_covAngle2-> reserve(ngoodPairs); 
            Pphotrail_S9overS9minusS1S2->clear(); 
            Pphotrail_S9overS9minusS1S2-> reserve(ngoodPairs); 
            Pphotrail_etawidth->clear(); 
            Pphotrail_etawidth-> reserve(ngoodPairs); 
            Pphotrail_phiwidth->clear(); 
            Pphotrail_phiwidth-> reserve(ngoodPairs); 
            Pphotrail_sigieta->clear(); 
            Pphotrail_sigieta-> reserve(ngoodPairs); 
            Pphotrail_SCbr->clear(); 
            Pphotrail_SCbr-> reserve(ngoodPairs); 
            Pphotrail_seedEnergy->clear();
            Pphotrail_seedEnergy-> reserve(ngoodPairs);
            Pphotrail_seedTime->clear();
            Pphotrail_seedTime-> reserve(ngoodPairs);            
            Pphotrail_SCphi->clear();
            Pphotrail_SCphi-> reserve(ngoodPairs);            
            Pphotrail_SCEtraw->clear();
            Pphotrail_SCEtraw-> reserve(ngoodPairs);            
            Pphotrail_SCEt->clear();
            Pphotrail_SCEt-> reserve(ngoodPairs);            
            Pphotrail_SCr9->clear();
            Pphotrail_SCr9-> reserve(ngoodPairs);            
            Pphotrail_SCnbBC->clear();
            Pphotrail_SCnbBC-> reserve(ngoodPairs);            
            Pphotrail_SCnXtal->clear();
            Pphotrail_SCnXtal-> reserve(ngoodPairs);            
            Pphotrail_isMatchingWithHLTObject->clear();
            Pphotrail_isMatchingWithHLTObject-> reserve(ngoodPairs);            
            Pphotrail_isConverted->clear();
            Pphotrail_isConverted-> reserve(ngoodPairs);            
            Pphotrail_NtrackConv->clear();
            Pphotrail_NtrackConv-> reserve(ngoodPairs);            
            Pphotrail_convEoverP->clear();
            Pphotrail_convEoverP-> reserve(ngoodPairs);            
            Pphotrail_convMass->clear();
            Pphotrail_convMass-> reserve(ngoodPairs);            
            Pphotrail_convCotanTheta->clear();
            Pphotrail_convCotanTheta-> reserve(ngoodPairs);            
            Pphotrail_convLikely->clear();
            Pphotrail_convLikely-> reserve(ngoodPairs);            
            Pphotrail_convVertexX->clear();
            Pphotrail_convVertexX-> reserve(ngoodPairs);            
            Pphotrail_convVertexY->clear();
            Pphotrail_convVertexY-> reserve(ngoodPairs);            
            Pphotrail_convVertexZ->clear();
            Pphotrail_convVertexZ-> reserve(ngoodPairs);
            Pphotrail_fBrem->clear();
            Pphotrail_fBrem-> reserve(ngoodPairs);            
            Pphotrail_isAspike->clear();
            Pphotrail_isAspike-> reserve(ngoodPairs);            
            Pphotrail_xVertex->clear();
            Pphotrail_xVertex-> reserve(ngoodPairs);            
            Pphotrail_yVertex->clear();
            Pphotrail_yVertex-> reserve(ngoodPairs);            
            Pphotrail_zVertex->clear();
            Pphotrail_zVertex-> reserve(ngoodPairs);            

            Pphotrail_HcalIso->clear(); 
            Pphotrail_HcalIso-> reserve(ngoodPairs); 
            Pphotrail_EcalIso->clear(); 
            Pphotrail_EcalIso-> reserve(ngoodPairs); 
            Pphotrail_TrackerIso->clear(); 
            Pphotrail_TrackerIso-> reserve(ngoodPairs); 
            Pphotrail_HcalIso_MIT03->clear(); 
            Pphotrail_HcalIso_MIT03-> reserve(ngoodPairs); 
            Pphotrail_EcalIso_MIT03->clear(); 
            Pphotrail_EcalIso_MIT03-> reserve(ngoodPairs); 
            Pphotrail_TrackerIso_MIT03->clear(); 
            Pphotrail_TrackerIso_MIT03-> reserve(ngoodPairs); 
            Pphotrail_HcalEcal_MIT03->clear(); 
            Pphotrail_HcalEcal_MIT03-> reserve(ngoodPairs); 
            Pphotrail_AbsTrackerIso->clear(); 
            Pphotrail_AbsTrackerIso-> reserve(ngoodPairs); 
            Pphotrail_HcalIsodR03->clear(); 
            Pphotrail_HcalIsodR03-> reserve(ngoodPairs); 
            Pphotrail_EcalIsodR03->clear(); 
            Pphotrail_EcalIsodR03-> reserve(ngoodPairs); 
            Pphotrail_TrackerIsodR03->clear(); 
            Pphotrail_TrackerIsodR03-> reserve(ngoodPairs); 
            Pphotrail_HcalIsoPerso->clear(); 
            Pphotrail_HcalIsoPerso-> reserve(ngoodPairs); 
            Pphotrail_EcalIsoPerso->clear(); 
            Pphotrail_EcalIsoPerso-> reserve(ngoodPairs); 
            Pphotrail_hoe->clear(); 
            Pphotrail_hoe-> reserve(ngoodPairs);
 
            Pphotrail_Cat->clear(); 
            Pphotrail_Cat -> reserve(ngoodPairs);
            Pphotrail_HasPixSeed->clear(); 
            Pphotrail_HasPixSeed -> reserve(ngoodPairs);
            Pphotrail_seedSeverity->clear(); 
            Pphotrail_seedSeverity-> reserve(ngoodPairs); 
            Pphotrail_recoFlag->clear(); 
            Pphotrail_recoFlag-> reserve(ngoodPairs); 
            Pphotrail_isEB->clear(); 
            Pphotrail_isEB-> reserve(ngoodPairs); 
            Pphotrail_isEE->clear(); 
            Pphotrail_isEE-> reserve(ngoodPairs); 





	    //inputEventTree->GetEvent(ievt);
	    if( (ievt%10==0 && ievt<=100)  || (ievt%100==0 && ievt<=1000)   || (ievt%1000==0 && ievt>1000)  )
	    {
		  cout <<"Analyzing "<< ievt << "th event: " << endl;
	    }
	      
           
            if (doHLT){
			
			if (nbHlt > 0) {if (event->hltAccept(ListWantedHLTnames[0])) dipho_HLT_bit0 = 1; else dipho_HLT_bit0 = 0;}
			if (nbHlt > 1) {if (event->hltAccept(ListWantedHLTnames[1])) dipho_HLT_bit1 = 1; else dipho_HLT_bit1 = 0;}
			if (nbHlt > 2) {if (event->hltAccept(ListWantedHLTnames[2])) dipho_HLT_bit2 = 1; else dipho_HLT_bit2 = 0;}
			if (nbHlt > 3) {if (event->hltAccept(ListWantedHLTnames[3])) dipho_HLT_bit3 = 1; else dipho_HLT_bit3 = 0;}
			if (nbHlt > 4) {if (event->hltAccept(ListWantedHLTnames[4])) dipho_HLT_bit4 = 1; else dipho_HLT_bit4 = 0;}
			if (nbHlt > 5) {if (event->hltAccept(ListWantedHLTnames[5])) dipho_HLT_bit5 = 1; else dipho_HLT_bit5 = 0;}
			if (nbHlt > 6) {if (event->hltAccept(ListWantedHLTnames[6])) dipho_HLT_bit6 = 1; else dipho_HLT_bit6 = 0;}
			if (nbHlt > 7) {if (event->hltAccept(ListWantedHLTnames[7])) dipho_HLT_bit7 = 1; else dipho_HLT_bit7 = 0;}
			if (nbHlt > 8) {if (event->hltAccept(ListWantedHLTnames[8])) dipho_HLT_bit8 = 1; else dipho_HLT_bit8 = 0;}
			if (nbHlt > 9) {if (event->hltAccept(ListWantedHLTnames[9])) dipho_HLT_bit9 = 1; else dipho_HLT_bit9 = 0;}
			if (nbHlt > 10) {if (event->hltAccept(ListWantedHLTnames[10])) dipho_HLT_bit10 = 1; else dipho_HLT_bit10 = 0;}
			if (nbHlt > 11) {if (event->hltAccept(ListWantedHLTnames[11])) dipho_HLT_bit11 = 1; else dipho_HLT_bit11 = 0;}
			if (nbHlt > 12) {if (event->hltAccept(ListWantedHLTnames[12])) dipho_HLT_bit12 = 1; else dipho_HLT_bit12 = 0;}
			if (nbHlt > 13) {if (event->hltAccept(ListWantedHLTnames[13])) dipho_HLT_bit13 = 1; else dipho_HLT_bit13 = 0;}
			if (nbHlt > 14) {if (event->hltAccept(ListWantedHLTnames[14])) dipho_HLT_bit14 = 1; else dipho_HLT_bit14 = 0;}
			if (nbHlt > 15) {if (event->hltAccept(ListWantedHLTnames[15])) dipho_HLT_bit15 = 1; else dipho_HLT_bit15 = 0;}
			if (nbHlt > 16) {if (event->hltAccept(ListWantedHLTnames[16])) dipho_HLT_bit16 = 1; else dipho_HLT_bit16 = 0;}
			if (nbHlt > 17) {if (event->hltAccept(ListWantedHLTnames[17])) dipho_HLT_bit17 = 1; else dipho_HLT_bit17 = 0;}
			if (nbHlt > 18) {if (event->hltAccept(ListWantedHLTnames[18])) dipho_HLT_bit18 = 1; else dipho_HLT_bit18 = 0;}
			if (nbHlt > 19) {if (event->hltAccept(ListWantedHLTnames[19])) dipho_HLT_bit19 = 1; else dipho_HLT_bit19 = 0;}
			if (nbHlt > 20) {if (event->hltAccept(ListWantedHLTnames[20])) dipho_HLT_bit20 = 1; else dipho_HLT_bit20 = 0;}
			if (nbHlt > 21) {if (event->hltAccept(ListWantedHLTnames[21])) dipho_HLT_bit21 = 1; else dipho_HLT_bit21 = 0;}
			if (nbHlt > 22) {if (event->hltAccept(ListWantedHLTnames[22])) dipho_HLT_bit22 = 1; else dipho_HLT_bit22 = 0;}
			if (nbHlt > 23) {if (event->hltAccept(ListWantedHLTnames[23])) dipho_HLT_bit23 = 1; else dipho_HLT_bit23 = 0;}
			if (nbHlt > 24) {if (event->hltAccept(ListWantedHLTnames[24])) dipho_HLT_bit24 = 1; else dipho_HLT_bit24 = 0;}
			if (nbHlt > 25) {if (event->hltAccept(ListWantedHLTnames[25])) dipho_HLT_bit25 = 1; else dipho_HLT_bit25 = 0;}
			if (nbHlt > 26) {if (event->hltAccept(ListWantedHLTnames[26])) dipho_HLT_bit26 = 1; else dipho_HLT_bit26 = 0;}
			if (nbHlt > 27) {if (event->hltAccept(ListWantedHLTnames[27])) dipho_HLT_bit27 = 1; else dipho_HLT_bit27 = 0;}
			if (nbHlt > 28) {if (event->hltAccept(ListWantedHLTnames[28])) dipho_HLT_bit28 = 1; else dipho_HLT_bit28 = 0;}
			if (nbHlt > 29) {if (event->hltAccept(ListWantedHLTnames[29])) dipho_HLT_bit29 = 1; else dipho_HLT_bit29 = 0;}
			if (nbHlt > 30) {if (event->hltAccept(ListWantedHLTnames[30])) dipho_HLT_bit30 = 1; else dipho_HLT_bit30 = 0;}
			if (nbHlt > 31) {if (event->hltAccept(ListWantedHLTnames[31])) dipho_HLT_bit31 = 1; else dipho_HLT_bit31 = 0;}
			if (nbHlt > 32) {if (event->hltAccept(ListWantedHLTnames[32])) dipho_HLT_bit32 = 1; else dipho_HLT_bit32 = 0;}
			if (nbHlt > 33) {if (event->hltAccept(ListWantedHLTnames[33])) dipho_HLT_bit33 = 1; else dipho_HLT_bit33 = 0;}
			if (nbHlt > 34) {if (event->hltAccept(ListWantedHLTnames[34])) dipho_HLT_bit34 = 1; else dipho_HLT_bit34 = 0;}
			if (nbHlt > 35) {if (event->hltAccept(ListWantedHLTnames[35])) dipho_HLT_bit35 = 1; else dipho_HLT_bit35 = 0;}
			if (nbHlt > 36) {if (event->hltAccept(ListWantedHLTnames[36])) dipho_HLT_bit36 = 1; else dipho_HLT_bit36 = 0;}
			if (nbHlt > 37) {if (event->hltAccept(ListWantedHLTnames[37])) dipho_HLT_bit37 = 1; else dipho_HLT_bit37 = 0;}
			if (nbHlt > 38) {if (event->hltAccept(ListWantedHLTnames[38])) dipho_HLT_bit38 = 1; else dipho_HLT_bit38 = 0;}
			if (nbHlt > 39) {if (event->hltAccept(ListWantedHLTnames[39])) dipho_HLT_bit39 = 1; else dipho_HLT_bit39 = 0;}
			if (nbHlt > 40) {if (event->hltAccept(ListWantedHLTnames[40])) dipho_HLT_bit40 = 1; else dipho_HLT_bit40 = 0;}
			if (nbHlt > 41) {if (event->hltAccept(ListWantedHLTnames[41])) dipho_HLT_bit41 = 1; else dipho_HLT_bit41 = 0;}
			if (nbHlt > 42) {if (event->hltAccept(ListWantedHLTnames[42])) dipho_HLT_bit42 = 1; else dipho_HLT_bit42 = 0;}
			if (nbHlt > 43) {if (event->hltAccept(ListWantedHLTnames[43])) dipho_HLT_bit43 = 1; else dipho_HLT_bit43 = 0;}
			if (nbHlt > 44) {if (event->hltAccept(ListWantedHLTnames[44])) dipho_HLT_bit44 = 1; else dipho_HLT_bit44 = 0;}
			if (nbHlt > 45) {if (event->hltAccept(ListWantedHLTnames[45])) dipho_HLT_bit45 = 1; else dipho_HLT_bit45 = 0;}
		}





		vector<pair <TRootPhoton*, TRootPhoton*> > theDiphotonPairs; //lead, trail
                pair <TRootPhoton*, TRootPhoton*> theDiphotonPair;

		for (int iphoton=0; iphoton< NbPhoton ; iphoton++){
			NbPhotons++;
			TRootPhoton *myphoton = (TRootPhoton*) photons->At(iphoton);
			//if (!(photonPassingPreselection(myphoton))) continue;
			for (int iphoton2=(iphoton+1) ; iphoton2 < NbPhoton ; iphoton2++){
				TRootPhoton *myphoton2 = (TRootPhoton*) photons->At(iphoton2);
				//if (!(photonPassingPreselection(myphoton2))) continue;
				//pair <TRootPhoton*, TRootPhoton*> theDiphotonPair;
				if (myphoton->Et()>myphoton2->Et()) theDiphotonPair = make_pair(myphoton, myphoton2);
				else theDiphotonPair = make_pair(myphoton2, myphoton);
                                saveThisEvent(event, theDiphotonPair, photons, vertices);          
				theDiphotonPairs.push_back(theDiphotonPair);
                          
			}
                        
		}
              
                myTree_->Fill();
                NgoodEvents++;
                cout << endl;

}

                cout << "nb of good events = " << NgoodEvents << endl;                      
                cout << "nb of good pairs = "  << NgoodPairs << endl;         
                cout << "nb of total photons = " << NbPhotons <<endl;
























/*
		// now run over diphoton to see if pass CiC	
		float theMaxSumEt = 0; int  theMaxSumEtInd = -10; int nbOfGood = 0;
		for (int i = 0 ; i < theDiphotonPairs.size() ; i++){
			int leadCiCStop, trailCiCStop;
			phoPair++;
		TLorentzVector Plead, Ptrail, Psum;
		Plead.SetPxPyPzE((theDiphotonPairs[i].first)->Px(),(theDiphotonPairs[i].first)->Py(),(theDiphotonPairs[i].first)->Pz(),(theDiphotonPairs[i].first)->Energy());
		Ptrail.SetPxPyPzE((theDiphotonPairs[i].second)->Px(),(theDiphotonPairs[i].second)->Py(),(theDiphotonPairs[i].second)->Pz(),(theDiphotonPairs[i].second)->Energy());
		Psum = Plead + Ptrail;
			int theDiphoCat = findTheDiphoCat(*(theDiphotonPairs[i].first), *(theDiphotonPairs[i].second));
			int phoLeadCat = findThePhoCat(*(theDiphotonPairs[i].first)); 
			cout << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
			cout << "Pair event Number=" << event->eventId() << "   runNumber=" << event->runId() << "   LS=" << event->luminosityBlock() <<"    cat of the pair=" << theDiphoCat << "   selectedVertex=0 "<< endl;
			bool leadPassingCiC = photonIsPassingCIC(*(theDiphotonPairs[i].first), vertices, tracks, *beamSpot, electrons, &leadCiCStop, phoLeadCat);
		//	cout << leadCiCStop << endl;
			switch (leadCiCStop){
				case 7 : phoPairLeadCiC7++; 
				case 6 : phoPairLeadCiC6++; 
				case 5 : phoPairLeadCiC5++;
				case 4 : phoPairLeadCiC4++;
				case 3 : phoPairLeadCiC3++;
				case 2 : phoPairLeadCiC2++;
				case 1 : phoPairLeadCiC1++;
				case 0 : phoPairLeadCiC0++;
			}
			int phoTrailCat = findThePhoCat(*(theDiphotonPairs[i].second)); 
			bool trailPassingCiC = photonIsPassingCIC(*(theDiphotonPairs[i].second), vertices, tracks, *beamSpot, electrons, &trailCiCStop, phoTrailCat);	
			switch (trailCiCStop){
				case 7 : phoPairTrailCiC7++;
				case 6 : phoPairTrailCiC6++;
				case 5 : phoPairTrailCiC5++;
				case 4 : phoPairTrailCiC4++;
				case 3 : phoPairTrailCiC3++;
				case 2 : phoPairTrailCiC2++;
				case 1 : phoPairTrailCiC1++;
				case 0 : phoPairTrailCiC0++; 
			}
			if (!(leadPassingCiC&&((theDiphotonPairs[i].first)->Et()>40))) continue;//test if lead pass CiC
			if (!(trailPassingCiC&&((theDiphotonPairs[i].second)->Et()>30))) continue; // tets if trail pass CiC
			phoAfterCutEt++;
	//		cout << "le second passe " << endl;
			nbOfGood++;
			float sumEt = (theDiphotonPairs[i].first)->Et() + (theDiphotonPairs[i].second)->Et();
			if (sumEt > theMaxSumEt) {
			//	cout << "Et lead " << (theDiphotonPairs[i].first)->Et() << " photon trail " << (theDiphotonPairs[i].second)->Et() << endl;
				theMaxSumEt = sumEt;
				theMaxSumEtInd = i;
			}
			//cout << "good pair in the event  " << nbOfGood << endl;
		}
		if ( theMaxSumEtInd >= 0 ) {
			phoGoodPair++;
			saveThisEvent(event, theDiphotonPairs[theMaxSumEtInd], photons, vertices);
		}
	}
	cout << "phopair " << phoPair << endl;
	cout << "phoPairLeadCiC0 = " << phoPairLeadCiC0 << endl;
	cout << "phoPairLeadCiC1 = " << phoPairLeadCiC1 << endl;
	cout << "phoPairLeadCiC2 = " << phoPairLeadCiC2 << endl;
	cout << "phoPairLeadCiC3 = " << phoPairLeadCiC3 << endl;
	cout << "phoPairLeadCiC4 = " << phoPairLeadCiC4 << endl;
	cout << "phoPairLeadCiC5 = " << phoPairLeadCiC5 << endl;
	cout << "phoPairLeadCiC6 = " << phoPairLeadCiC6 << endl;
	cout << "phoPairLeadCiC7 = " << phoPairLeadCiC7 << endl;
	cout << "phoPairTrailCiC0 = " << phoPairTrailCiC0 << endl;
	cout << "phoPairTrailCiC1 = " << phoPairTrailCiC1 << endl;
	cout << "phoPairTrailCiC2 = " << phoPairTrailCiC2 << endl;
	cout << "phoPairTrailCiC3 = " << phoPairTrailCiC3 << endl;
	cout << "phoPairTrailCiC4 = " << phoPairTrailCiC4 << endl;
	cout << "phoPairTrailCiC5 = " << phoPairTrailCiC5 << endl;
	cout << "phoPairTrailCiC6 = " << phoPairTrailCiC6 << endl;
	cout << "phoPairTrailCiC7 = " << phoPairTrailCiC7 << endl;
	cout << "phoAfterCutEt = " << phoAfterCutEt << endl;
	cout << "phoGoodPair = " << phoGoodPair << endl;
*/



	endMacro();
}
