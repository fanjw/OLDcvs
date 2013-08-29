#include "functions.h"
#include "functions.cc"


#include "TF1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"
#include <utility>

#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootBardak.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootBeamSpot.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootBeamStatus.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootCluster.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootDummyEvent.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootEcalRecHit.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootElectron.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootEvent.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootJet.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootMCParticle.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootMCPhoton.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootMET.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootMuon.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootParticle.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootPhoton.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootRun.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootSignalEvent.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootSuperCluster.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootTopTop.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootTrack.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootVertex.h"
#include "/sps/cms/jfan/CMSSW_5_2_4_patch4/src/Morgan/IpnTreeProducer/interface/TRootHLTObject.h"


TFile *myFile;// = new TFile("theMiniTree.root","RECREATE");
TTree *myTree_;
TTree *myTree;
TChain *inputEventTree = new TChain("eventTree");
TChain *inputRunTree = new TChain("runTree");

//string ListWantedHLTnames[13] = {"HLT_DoublePhoton33_v1","HLT_Photon125_NoSpikeFilter_v1","HLT_Photon20_R9Id_Photon18_R9Id_v1","HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1","HLT_Photon26_CaloIdL_IsoVL_Photon18_v1","HLT_Photon26_IsoVL_Photon18_IsoVL_v1","HLT_Photon26_IsoVL_Photon18_v1","HLT_Photon26_Photon18_v1","HLT_Photon30_CaloIdVL_IsoL_v1","HLT_Photon30_CaloIdVL_v1","HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1","HLT_Photon75_CaloIdVL_IsoL_v1","HLT_Photon75_CaloIdVL_v1"};

//string ListWantedHLTnames[18] = {"HLT_DoublePhoton33_v2","HLT_Photon125_NoSpikeFilter_v2","HLT_Photon20_CaloIdVL_IsoL_v1","HLT_Photon20_R9Id_Photon18_R9Id_v2","HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2","HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1","HLT_Photon26_CaloIdL_IsoVL_Photon18_v2","HLT_Photon26_IsoVL_Photon18_IsoVL_v2","HLT_Photon26_IsoVL_Photon18_v2","HLT_Photon26_Photon18_v2","HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1","HLT_Photon30_CaloIdVL_IsoL_v2","HLT_Photon30_CaloIdVL_v2","HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2","HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1","HLT_Photon50_CaloIdVL_IsoL_v1","HLT_Photon75_CaloIdVL_IsoL_v2","HLT_Photon75_CaloIdVL_v2"};

//string ListWantedHLTnames[18] = {"HLT_DoublePhoton33_v3","HLT_Photon125_NoSpikeFilter_v3","HLT_Photon20_CaloIdVL_IsoL_v2","HLT_Photon20_R9Id_Photon18_R9Id_v3","HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3","HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2","HLT_Photon26_CaloIdL_IsoVL_Photon18_v3","HLT_Photon26_IsoVL_Photon18_IsoVL_v3","HLT_Photon26_IsoVL_Photon18_v3","HLT_Photon26_Photon18_v3","HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2","HLT_Photon30_CaloIdVL_IsoL_v3","HLT_Photon30_CaloIdVL_v3","HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3","HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2","HLT_Photon50_CaloIdVL_IsoL_v2","HLT_Photon75_CaloIdVL_IsoL_v3","HLT_Photon75_CaloIdVL_v3"};

/*string ListWantedHLTnames[34] = { // Run2011 1e33 v1.3 
"HLT_DoubleEle33_CaloIdL_v1",
"HLT_DoubleEle33_v1",
"HLT_DoublePhoton33_HEVT_v1",
"HLT_DoublePhoton33_v4",
"HLT_DoublePhoton40_MR150_v1",
"HLT_DoublePhoton40_R014_MR150_v1",
"HLT_DoublePhoton50_v1",
"HLT_DoublePhoton5_IsoVL_CEP_v3",
"HLT_DoublePhoton60_v1",
"HLT_Photon125_v1",
"HLT_Photon200_NoHE_v1",
"HLT_Photon20_CaloIdVL_IsoL_v3",
"HLT_Photon20_R9Id_Photon18_R9Id_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v4",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v4",
"HLT_Photon26_IsoVL_Photon18_v4",
"HLT_Photon26_Photon18_v4",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3",
"HLT_Photon26_R9Id_Photon18_R9Id_v1",
"HLT_Photon30_CaloIdVL_IsoL_v4",
"HLT_Photon30_CaloIdVL_v4",
"HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v1",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3",
"HLT_Photon36_IsoVL_Photon22_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1",
"HLT_Photon50_CaloIdVL_IsoL_v3",
"HLT_Photon50_CaloIdVL_v1",
"HLT_Photon75_CaloIdVL_IsoL_v4",
"HLT_Photon75_CaloIdVL_v4",
"HLT_Photon90_CaloIdVL_IsoL_v1",
"HLT_Photon90_CaloIdVL_IsoL_v1"
};*/


/*string ListWantedHLTnames[38] = { // Run2011 1e33 v2.3
"HLT_DoubleEle33_CaloIdL_v2",
"HLT_DoubleEle33_v2",
"HLT_DoublePhoton33_HEVT_v2",
"HLT_DoublePhoton33_v5",
"HLT_DoublePhoton40_MR150_v3",
"HLT_DoublePhoton40_R014_MR150_v3",
"HLT_DoublePhoton50_v2",
"HLT_DoublePhoton5_IsoVL_CEP_v4",
"HLT_DoublePhoton60_v2",
"HLT_Photon125_v2",
"HLT_Photon200_NoHE_v2",
"HLT_Photon20_CaloIdVL_IsoL_v4",
"HLT_Photon20_R9Id_Photon18_R9Id_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v5",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v5",
"HLT_Photon26_IsoVL_Photon18_v5",
"HLT_Photon26_Photon18_v5",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4",
"HLT_Photon26_R9Id_Photon18_R9Id_v2",
"HLT_Photon30_CaloIdVL_IsoL_v5",
"HLT_Photon30_CaloIdVL_v5",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v2",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4",
"HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1",
"HLT_Photon36_IsoVL_Photon22_v2",
"HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_R9Id_Photon22_R9Id_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2",
"HLT_Photon50_CaloIdVL_IsoL_v4",
"HLT_Photon50_CaloIdVL_v2",
"HLT_Photon75_CaloIdVL_IsoL_v5",
"HLT_Photon75_CaloIdVL_v5",
"HLT_Photon90_CaloIdVL_IsoL_v2",
"HLT_Photon90_CaloIdVL_v2"
};*/


/*string ListWantedHLTnames[38] = { // Run2011 1e33 v2.4
"HLT_DoubleEle33_CaloIdL_v2",
"HLT_DoubleEle33_v2",
"HLT_DoublePhoton33_HEVT_v2",
"HLT_DoublePhoton33_v5",
"HLT_DoublePhoton40_MR150_v3",
"HLT_DoublePhoton40_R014_MR150_v3",
"HLT_DoublePhoton50_v2",
"HLT_DoublePhoton5_IsoVL_CEP_v4",
"HLT_DoublePhoton60_v2",
"HLT_Photon125_v2",
"HLT_Photon200_NoHE_v2",
"HLT_Photon20_CaloIdVL_IsoL_v4",
"HLT_Photon20_R9Id_Photon18_R9Id_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v5",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v5",
"HLT_Photon26_IsoVL_Photon18_v5",
"HLT_Photon26_Photon18_v5",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4",
"HLT_Photon26_R9Id_Photon18_R9Id_v2",
"HLT_Photon30_CaloIdVL_IsoL_v5",
"HLT_Photon30_CaloIdVL_v5",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v2",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4",
"HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1",
"HLT_Photon36_IsoVL_Photon22_v2",
"HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_R9Id_Photon22_R9Id_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2",
"HLT_Photon50_CaloIdVL_IsoL_v4",
"HLT_Photon50_CaloIdVL_v2",
"HLT_Photon75_CaloIdVL_IsoL_v5",
"HLT_Photon75_CaloIdVL_v5",
"HLT_Photon90_CaloIdVL_IsoL_v2",
"HLT_Photon90_CaloIdVL_v2"
};*/

//string ListWantedHLTnames[12] = {"HLT_DoublePhoton22_L1R_v1","HLT_DoublePhoton17_SingleIsol_L1R_v1","HLT_Photon20_Cleaned_L1R","HLT_Photon30_Cleaned_L1R","HLT_Photon40_Isol_Cleaned_L1R_v1","HLT_DoublePhoton5_CEP_L1R_v3","HLT_Photon110_NoHE_Cleaned_L1R_v1","HLT_Photon17_Isol_SC17HE_L1R_v1","HLT_Photon22_SC22HE_L1R_v1","HLT_Photon40_CaloId_Cleaned_L1R_v1","HLT_Photon50_Cleaned_L1R_v1","HLT_Photon70_Cleaned_L1R_v1"};

string ListWantedHLTnames[44] = { // Monte Carlo 1e33
"HLT_DoubleEle33_CaloIdL_v3",
"HLT_DoubleEle33_v3",
"HLT_DoubleEle45_CaloIdL_v2",
"HLT_DoublePhoton33_HEVT_v3",
"HLT_DoublePhoton38_HEVT_v2",
"HLT_DoublePhoton40_MR150_v4",
"HLT_DoublePhoton40_R014_MR150_v4",	
"HLT_DoublePhoton5_IsoVL_CEP_v5",
"HLT_DoublePhoton60_v3",
"HLT_DoublePhoton80_v1",
"HLT_Photon135_v1",
"HLT_Photon200_NoHE_v3",
"HLT_Photon20_CaloIdVL_IsoL_v5",
"HLT_Photon20_R9Id_Photon18_R9Id_v6",
"HLT_Photon225_NoHE_v1",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v6",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v6",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v6",
"HLT_Photon26_IsoVL_Photon18_v6",
"HLT_Photon26_Photon18_v6",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v5",	
"HLT_Photon26_R9Id_Photon18_R9Id_v5",
"HLT_Photon30_CaloIdVL_IsoL_v6",
"HLT_Photon30_CaloIdVL_v6",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v2",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v2",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v3",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v5",
"HLT_Photon36_CaloIdVL_Photon22_CaloIdVL_v1",
"HLT_Photon36_IsoVL_Photon22_v3",
"HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v2",
"HLT_Photon36_R9Id_Photon22_R9Id_v2",
"HLT_Photon400_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v3",
"HLT_Photon44_CaloIdL_Photon34_CaloIdL_v1",
"HLT_Photon48_CaloIdL_Photon38_CaloIdL_v1",
"HLT_Photon50_CaloIdVL_IsoL_v5",
"HLT_Photon50_CaloIdVL_v3",
"HLT_Photon75_CaloIdVL_IsoL_v6",	
"HLT_Photon75_CaloIdVL_v6",
"HLT_Photon90_CaloIdVL_IsoL_v3",
"HLT_Photon90_CaloIdVL_v3"
};
int nbHlt = 44; 







float secondPhotonCut = 20.0;
float IsoTrk = 6.0;
float IsoEcal = 4.0;
float IsoHcal = 4.0;


TString theHTLobject = "hltPhoton26IsoVLTrackIsoFilter";

  bool doHLT;
  bool doHLTobject;
  bool doMC;
  bool doJetMC;
  bool doMETMC;
  bool doPDFInfo;
  bool doDiphotons;
  bool doLeadingPhoton;
  bool doSignalMuMuGamma;
  bool doSignalTopTop;
  bool doPhotonConversionMC;
  bool doBeamSpot;
  bool doPrimaryVertex;
  bool doZeePrimaryVertex;
  bool doTrack; 
  bool doJet;
  bool doMuon;
  bool doElectron;
  bool doPhoton;
  bool doCluster;
  bool doPhotonConversion;
  bool doMET;
  bool doBardak;
  bool doPhotonVertexCorrection;
  bool doPhotonIsolation;   
  bool doElectronConversionMC;
  bool doWorstIsolation; 

int dipho_HLT_bit0;
int dipho_HLT_bit1;
int dipho_HLT_bit2;
int dipho_HLT_bit3;
int dipho_HLT_bit4;
int dipho_HLT_bit5;
int dipho_HLT_bit6;
int dipho_HLT_bit7;
int dipho_HLT_bit8;
int dipho_HLT_bit9;
int dipho_HLT_bit10;
int dipho_HLT_bit11;
int dipho_HLT_bit12;
int dipho_HLT_bit13;
int dipho_HLT_bit14;
int dipho_HLT_bit15;
int dipho_HLT_bit16;
int dipho_HLT_bit17;
int dipho_HLT_bit18;
int dipho_HLT_bit19;
int dipho_HLT_bit20;
int dipho_HLT_bit21;
int dipho_HLT_bit22;
int dipho_HLT_bit23;
int dipho_HLT_bit24;
int dipho_HLT_bit25;
int dipho_HLT_bit26;
int dipho_HLT_bit27;
int dipho_HLT_bit28;
int dipho_HLT_bit29;
int dipho_HLT_bit30;
int dipho_HLT_bit31;
int dipho_HLT_bit32;
int dipho_HLT_bit33;
int dipho_HLT_bit34;
int dipho_HLT_bit35;
int dipho_HLT_bit36;
int dipho_HLT_bit37;
int dipho_HLT_bit38;
int dipho_HLT_bit39;
int dipho_HLT_bit40;
int dipho_HLT_bit41;
int dipho_HLT_bit42;
int dipho_HLT_bit43;
int dipho_HLT_bit44;
int dipho_HLT_bit45;
// event infos

// miniTree de sortie

//here the event infos
int event_number;    
int event_runNumber;                                    
int event_LumiSection; 
float event_eventPtHat; 
int event_processId;
int event_nRecoVertex;
int event_nGenInTimeVertex;
int event_nGenOutOfTimeVertex;
int event_nPhotons;
float  event_rho;
int event_nPairs;
int ngoodPairs;



// diphoton kine
std::vector<float> dipho_mgg;
std::vector<float> *Pdipho_mgg = &dipho_mgg;
std::vector<float> dipho_qt;
std::vector<float> *Pdipho_qt = &dipho_qt;
std::vector<float> dipho_ql;
std::vector<float> *Pdipho_ql = &dipho_ql;
std::vector<float> dipho_deltaR;
std::vector<float> *Pdipho_deltaR = &dipho_deltaR;
std::vector<float> dipho_costhetastar;
std::vector<float> *Pdipho_costhetastar = &dipho_costhetastar;
std::vector<float> dipho_eta;
std::vector<float> *Pdipho_eta = &dipho_eta;
std::vector<float> dipho_etastar;
std::vector<float> *Pdipho_etastar = &dipho_etastar;

// MC diphoton kine
std::vector<float> diphoMC_mgg;
std::vector<float> *PdiphoMC_mgg = &diphoMC_mgg;
std::vector<float> diphoMC_qt;
std::vector<float> *PdiphoMC_qt = &diphoMC_qt;
std::vector<float> diphoMC_ql;
std::vector<float> *PdiphoMC_ql = &diphoMC_ql;
std::vector<float> diphoMC_deltaR;
std::vector<float> *PdiphoMC_deltaR = &diphoMC_deltaR;
std::vector<float> diphoMC_costhetastar;
std::vector<float> *PdiphoMC_costhetastar = &diphoMC_costhetastar;
std::vector<float> diphoMC_eta;
std::vector<float> *PdiphoMC_eta = &diphoMC_eta;
std::vector<float> diphoMC_etastar;
std::vector<float> *PdiphoMC_etastar = &diphoMC_etastar;




//phoLead
//relatedMC
std::vector<int> pholead_isMatchingWithMC;
std::vector<int> *Ppholead_isMatchingWithMC = &pholead_isMatchingWithMC;
std::vector<int> pholead_GenId;
std::vector<int> *Ppholead_GenId = &pholead_GenId;
std::vector<int> pholead_MotherId;
std::vector<int> *Ppholead_MotherId = &pholead_MotherId;
std::vector<int> pholead_isPromptGenPho;
std::vector<int> *Ppholead_isPromptGenPho = &pholead_isPromptGenPho;
std::vector<int> pholead_isFromQuarkGen;
std::vector<int> *Ppholead_isFromQuarkGen = &pholead_isFromQuarkGen;
std::vector<int> pholead_isPi0Gen;
std::vector<int> *Ppholead_isPi0Gen = &pholead_isPi0Gen;
std::vector<int> pholead_isEtaGen;
std::vector<int> *Ppholead_isEtaGen = &pholead_isEtaGen;
std::vector<int> pholead_isRhoGen;
std::vector<int> *Ppholead_isRhoGen = &pholead_isRhoGen;
std::vector<int> pholead_isOmegaGen;
std::vector<int> *Ppholead_isOmegaGen = &pholead_isOmegaGen;
std::vector<int> pholead_isGenElectron;
std::vector<int> *Ppholead_isGenElectron = &pholead_isGenElectron;
std::vector<int> pholead_eventPassHLT_Photon10_L1R;
std::vector<int> *Ppholead_eventPassHLT_Photon10_L1R = &pholead_eventPassHLT_Photon10_L1R;
std::vector<int> pholead_eventPassHLT_Photon15_L1R;
std::vector<int> *Ppholead_eventPassHLT_Photon15_L1R = &pholead_eventPassHLT_Photon15_L1R;
std::vector<int> pholead_eventPassHLT_DoublePhoton10_L1R;
std::vector<int> *Ppholead_eventPassHLT_DoublePhoton10_L1R = &pholead_eventPassHLT_DoublePhoton10_L1R;
std::vector<float> pholead_PromptGenIsoEnergyStatus1_cone02;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus1_cone02 = &pholead_PromptGenIsoEnergyStatus1_cone02;
std::vector<float> pholead_PromptGenIsoEnergyStatus2_cone02;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus2_cone02 = &pholead_PromptGenIsoEnergyStatus2_cone02;
std::vector<float> pholead_PromptGenIsoEnergyStatus1_cone03;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus1_cone03 = &pholead_PromptGenIsoEnergyStatus1_cone03;
std::vector<float> pholead_PromptGenIsoEnergyStatus2_cone03;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus2_cone03 = &pholead_PromptGenIsoEnergyStatus2_cone03;
std::vector<float> pholead_PromptGenIsoEnergyStatus1_cone035;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus1_cone035 = &pholead_PromptGenIsoEnergyStatus1_cone035;
std::vector<float> pholead_PromptGenIsoEnergyStatus2_cone035;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus2_cone035 = &pholead_PromptGenIsoEnergyStatus2_cone035;
std::vector<float> pholead_PromptGenIsoEnergyStatus1_cone04;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus1_cone04 = &pholead_PromptGenIsoEnergyStatus1_cone04;
std::vector<float> pholead_PromptGenIsoEnergyStatus2_cone04;
std::vector<float> *Ppholead_PromptGenIsoEnergyStatus2_cone04 = &pholead_PromptGenIsoEnergyStatus2_cone04;
std::vector<float> pholead_trueE;
std::vector<float> *Ppholead_trueE = &pholead_trueE;
std::vector<float> pholead_truePx;
std::vector<float> *Ppholead_truePx = &pholead_truePx;
std::vector<float> pholead_truePy;
std::vector<float> *Ppholead_truePy = &pholead_truePy;
std::vector<float> pholead_truePz;
std::vector<float> *Ppholead_truePz = &pholead_truePz;
std::vector<float> pholead_trueEta;
std::vector<float> *Ppholead_trueEta = &pholead_trueEta;
std::vector<float> pholead_truePhi;
std::vector<float> *Ppholead_truePhi = &pholead_truePhi;
std::vector<int> pholead_MCisConverted;
std::vector<int> *Ppholead_MCisConverted = &pholead_MCisConverted;
std::vector<float> pholead_MCconvEoverP;
std::vector<float> *Ppholead_MCconvEoverP = &pholead_MCconvEoverP;
std::vector<float> pholead_MCconvCotanTheta;
std::vector<float> *Ppholead_MCconvCotanTheta = &pholead_MCconvCotanTheta;
std::vector<float> pholead_MCconvVertexX;
std::vector<float> *Ppholead_MCconvVertexX = &pholead_MCconvVertexX;
std::vector<float> pholead_MCconvVertexY;
std::vector<float> *Ppholead_MCconvVertexY = &pholead_MCconvVertexY;
std::vector<float> pholead_MCconvVertexZ;
std::vector<float> *Ppholead_MCconvVertexZ = &pholead_MCconvVertexZ;
std::vector<float> pholead_eleMCtruthBrem;
std::vector<float> *Ppholead_eleMCtruthBrem = &pholead_eleMCtruthBrem;
std::vector<int> pholead_eleMCtruthNBrem;
std::vector<int> *Ppholead_eleMCtruthNBrem = &pholead_eleMCtruthNBrem;




// kinevars
std::vector<float> pholead_et;
std::vector<float> *Ppholead_et = &pholead_et;
std::vector<float> pholead_eta;
std::vector<float> *Ppholead_eta = &pholead_eta;
std::vector<float> pholead_SCeta;
std::vector<float> *Ppholead_SCeta = &pholead_SCeta;
// cluster shape vars
std::vector<float> pholead_r9;
std::vector<float> *Ppholead_r9 = &pholead_r9;
std::vector<float> pholead_cPP;
std::vector<float> *Ppholead_cPP = &pholead_cPP;
std::vector<float> pholead_cEP;
std::vector<float> *Ppholead_cEP = &pholead_cEP;
std::vector<float> pholead_cEE;
std::vector<float> *Ppholead_cEE = &pholead_cEE;
std::vector<float> pholead_r19;
std::vector<float> *Ppholead_r19 = &pholead_r19;
std::vector<float> pholead_SCEraw;
std::vector<float> *Ppholead_SCEraw = &pholead_SCEraw;
std::vector<float> pholead_eMax;
std::vector<float> *Ppholead_eMax = &pholead_eMax;
std::vector<float> pholead_e2x2;
std::vector<float> *Ppholead_e2x2 = &pholead_e2x2;
std::vector<float> pholead_e5x5;
std::vector<float> *Ppholead_e5x5 = &pholead_e5x5;
std::vector<float> pholead_ratioSeed;
std::vector<float> *Ppholead_ratioSeed = &pholead_ratioSeed;
std::vector<float> pholead_ratioS4;
std::vector<float> *Ppholead_ratioS4 = &pholead_ratioS4;
std::vector<float> pholead_lambdaRatio;
std::vector<float> *Ppholead_lambdaRatio = &pholead_lambdaRatio;
std::vector<float> pholead_lamdbaDivCov;
std::vector<float> *Ppholead_lamdbaDivCov = &pholead_lamdbaDivCov;
std::vector<float> pholead_secondMomentMaj;
std::vector<float> *Ppholead_secondMomentMaj = &pholead_secondMomentMaj;
std::vector<float> pholead_secondMomentMin;
std::vector<float> *Ppholead_secondMomentMin = &pholead_secondMomentMin;
std::vector<float> pholead_secondMomentAlpha;
std::vector<float> *Ppholead_secondMomentAlpha = &pholead_secondMomentAlpha;
std::vector<float> pholead_covAngle;
std::vector<float> *Ppholead_covAngle = &pholead_covAngle;
std::vector<float> pholead_covAngle2;
std::vector<float> *Ppholead_covAngle2 = &pholead_covAngle2;
std::vector<float> pholead_S9overS9minusS1S2;
std::vector<float> *Ppholead_S9overS9minusS1S2 = &pholead_S9overS9minusS1S2;
std::vector<float> pholead_etawidth;
std::vector<float> *Ppholead_etawidth = &pholead_etawidth;
std::vector<float> pholead_phiwidth;
std::vector<float> *Ppholead_phiwidth = &pholead_phiwidth;
std::vector<float> pholead_sigieta;
std::vector<float> *Ppholead_sigieta = &pholead_sigieta;
std::vector<float> pholead_SCbr;
std::vector<float> *Ppholead_SCbr = &pholead_SCbr;
std::vector<float> pholead_seedEnergy;
std::vector<float> *Ppholead_seedEnergy = &pholead_seedEnergy;
std::vector<float> pholead_seedTime;
std::vector<float> *Ppholead_seedTime = &pholead_seedTime;
std::vector<float> pholead_SCphi;
std::vector<float> *Ppholead_SCphi = &pholead_SCphi;
std::vector<float> pholead_SCEtraw;
std::vector<float> *Ppholead_SCEtraw = &pholead_SCEtraw;
std::vector<float> pholead_SCEt;
std::vector<float> *Ppholead_SCEt = &pholead_SCEt;
std::vector<float> pholead_SCr9;
std::vector<float> *Ppholead_SCr9 = &pholead_SCr9;
std::vector<int> pholead_SCnbBC;
std::vector<int> *Ppholead_SCnbBC = &pholead_SCnbBC;
std::vector<int> pholead_SCnXtal;
std::vector<int> *Ppholead_SCnXtal = &pholead_SCnXtal;
std::vector<int> pholead_isMatchingWithHLTObject;
std::vector<int> *Ppholead_isMatchingWithHLTObject = &pholead_isMatchingWithHLTObject;
std::vector<int> pholead_isConverted;
std::vector<int> *Ppholead_isConverted = &pholead_isConverted;
std::vector<int> pholead_NtrackConv;
std::vector<int> *Ppholead_NtrackConv = &pholead_NtrackConv;
std::vector<float> pholead_convEoverP;
std::vector<float> *Ppholead_convEoverP = &pholead_convEoverP;
std::vector<float> pholead_convMass;
std::vector<float> *Ppholead_convMass = &pholead_convMass;
std::vector<float> pholead_convCotanTheta;
std::vector<float> *Ppholead_convCotanTheta = &pholead_convCotanTheta;
std::vector<float> pholead_MCconvMass;
std::vector<float> *Ppholead_MCconvMass = &pholead_MCconvMass;
std::vector<float> pholead_convLikely;
std::vector<float> *Ppholead_convLikely = &pholead_convLikely;
std::vector<float> pholead_convVertexX;
std::vector<float> *Ppholead_convVertexX = &pholead_convVertexX;
std::vector<float> pholead_convVertexY;
std::vector<float> *Ppholead_convVertexY = &pholead_convVertexY;
std::vector<float> pholead_convVertexZ;
std::vector<float> *Ppholead_convVertexZ = &pholead_convVertexZ;
std::vector<float> pholead_fBrem;
std::vector<float> *Ppholead_fBrem = &pholead_fBrem;
std::vector<int> pholead_isAspike;
std::vector<int> *Ppholead_isAspike = &pholead_isAspike;
std::vector<float> pholead_xVertex ;
std::vector<float> *Ppholead_xVertex  = &pholead_xVertex ;
std::vector<float> pholead_yVertex ;
std::vector<float> *Ppholead_yVertex  = &pholead_yVertex ;
std::vector<float> pholead_zVertex ;
std::vector<float> *Ppholead_zVertex  = &pholead_zVertex ;



//isolation variables
std::vector<float> pholead_HcalIso;
std::vector<float> *Ppholead_HcalIso = &pholead_HcalIso;
std::vector<float> pholead_EcalIso;
std::vector<float> *Ppholead_EcalIso = &pholead_EcalIso;
std::vector<float> pholead_TrackerIso;
std::vector<float> *Ppholead_TrackerIso = &pholead_TrackerIso;
std::vector<float> pholead_HcalIso_MIT03;
std::vector<float> *Ppholead_HcalIso_MIT03 = &pholead_HcalIso_MIT03;
std::vector<float> pholead_EcalIso_MIT03;
std::vector<float> *Ppholead_EcalIso_MIT03 = &pholead_EcalIso_MIT03;
std::vector<float> pholead_TrackerIso_MIT03;
std::vector<float> *Ppholead_TrackerIso_MIT03 = &pholead_TrackerIso_MIT03;
std::vector<float> pholead_HcalEcal_MIT03;
std::vector<float> *Ppholead_HcalEcal_MIT03 = &pholead_HcalEcal_MIT03;
std::vector<float> pholead_AbsTrackerIso;
std::vector<float> *Ppholead_AbsTrackerIso = &pholead_AbsTrackerIso;
std::vector<float> pholead_HcalIsodR03;
std::vector<float> *Ppholead_HcalIsodR03 = &pholead_HcalIsodR03;
std::vector<float> pholead_EcalIsodR03;
std::vector<float> *Ppholead_EcalIsodR03 = &pholead_EcalIsodR03;
std::vector<float> pholead_TrackerIsodR03;
std::vector<float> *Ppholead_TrackerIsodR03 = &pholead_TrackerIsodR03;
std::vector<float> pholead_HcalIsoPerso;
std::vector<float> *Ppholead_HcalIsoPerso = &pholead_HcalIsoPerso;
std::vector<float> pholead_EcalIsoPerso;
std::vector<float> *Ppholead_EcalIsoPerso = &pholead_EcalIsoPerso;
std::vector<float> pholead_hoe;
std::vector<float> *Ppholead_hoe = &pholead_hoe;

//other variables
std::vector<int> pholead_Cat;
std::vector<int> *Ppholead_Cat = &pholead_Cat;
std::vector<int> pholead_HasPixSeed;
std::vector<int> *Ppholead_HasPixSeed = &pholead_HasPixSeed;
std::vector<int> pholead_seedSeverity;
std::vector<int> *Ppholead_seedSeverity = &pholead_seedSeverity;
std::vector<int> pholead_recoFlag;
std::vector<int> *Ppholead_recoFlag = &pholead_recoFlag;
std::vector<int> pholead_isEB;
std::vector<int> *Ppholead_isEB = &pholead_isEB;
std::vector<int> pholead_isEE;
std::vector<int> *Ppholead_isEE = &pholead_isEE;






//phoTrail
//relatedMC
std::vector<int> photrail_isMatchingWithMC;
std::vector<int> *Pphotrail_isMatchingWithMC = &photrail_isMatchingWithMC;
std::vector<int> photrail_GenId;
std::vector<int> *Pphotrail_GenId = &photrail_GenId;
std::vector<int> photrail_MotherId;
std::vector<int> *Pphotrail_MotherId = &photrail_MotherId;
std::vector<int> photrail_isPromptGenPho;
std::vector<int> *Pphotrail_isPromptGenPho = &photrail_isPromptGenPho;
std::vector<int> photrail_isFromQuarkGen;
std::vector<int> *Pphotrail_isFromQuarkGen = &photrail_isFromQuarkGen;
std::vector<int> photrail_isPi0Gen;
std::vector<int> *Pphotrail_isPi0Gen = &photrail_isPi0Gen;
std::vector<int> photrail_isEtaGen;
std::vector<int> *Pphotrail_isEtaGen = &photrail_isEtaGen;
std::vector<int> photrail_isRhoGen;
std::vector<int> *Pphotrail_isRhoGen = &photrail_isRhoGen;
std::vector<int> photrail_isOmegaGen;
std::vector<int> *Pphotrail_isOmegaGen = &photrail_isOmegaGen;
std::vector<int> photrail_isGenElectron;
std::vector<int> *Pphotrail_isGenElectron = &photrail_isGenElectron;
std::vector<int> photrail_eventPassHLT_Photon10_L1R;
std::vector<int> *Pphotrail_eventPassHLT_Photon10_L1R = &photrail_eventPassHLT_Photon10_L1R;
std::vector<int> photrail_eventPassHLT_Photon15_L1R;
std::vector<int> *Pphotrail_eventPassHLT_Photon15_L1R = &photrail_eventPassHLT_Photon15_L1R;
std::vector<int> photrail_eventPassHLT_DoublePhoton10_L1R;
std::vector<int> *Pphotrail_eventPassHLT_DoublePhoton10_L1R = &photrail_eventPassHLT_DoublePhoton10_L1R;
std::vector<float> photrail_PromptGenIsoEnergyStatus1_cone02;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus1_cone02 = &photrail_PromptGenIsoEnergyStatus1_cone02;
std::vector<float> photrail_PromptGenIsoEnergyStatus2_cone02;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus2_cone02 = &photrail_PromptGenIsoEnergyStatus2_cone02;
std::vector<float> photrail_PromptGenIsoEnergyStatus1_cone03;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus1_cone03 = &photrail_PromptGenIsoEnergyStatus1_cone03;
std::vector<float> photrail_PromptGenIsoEnergyStatus2_cone03;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus2_cone03 = &photrail_PromptGenIsoEnergyStatus2_cone03;
std::vector<float> photrail_PromptGenIsoEnergyStatus1_cone035;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus1_cone035 = &photrail_PromptGenIsoEnergyStatus1_cone035;
std::vector<float> photrail_PromptGenIsoEnergyStatus2_cone035;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus2_cone035 = &photrail_PromptGenIsoEnergyStatus2_cone035;
std::vector<float> photrail_PromptGenIsoEnergyStatus1_cone04;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus1_cone04 = &photrail_PromptGenIsoEnergyStatus1_cone04;
std::vector<float> photrail_PromptGenIsoEnergyStatus2_cone04;
std::vector<float> *Pphotrail_PromptGenIsoEnergyStatus2_cone04 = &photrail_PromptGenIsoEnergyStatus2_cone04;
std::vector<float> photrail_trueE;
std::vector<float> *Pphotrail_trueE = &photrail_trueE;
std::vector<float> photrail_truePx;
std::vector<float> *Pphotrail_truePx = &photrail_truePx;
std::vector<float> photrail_truePy;
std::vector<float> *Pphotrail_truePy = &photrail_truePy;
std::vector<float> photrail_truePz;
std::vector<float> *Pphotrail_truePz = &photrail_truePz;
std::vector<float> photrail_trueEta;
std::vector<float> *Pphotrail_trueEta = &photrail_trueEta;
std::vector<float> photrail_truePhi;
std::vector<float> *Pphotrail_truePhi = &photrail_truePhi;
std::vector<int> photrail_MCisConverted;
std::vector<int> *Pphotrail_MCisConverted = &photrail_MCisConverted;
std::vector<float> photrail_MCconvEoverP;
std::vector<float> *Pphotrail_MCconvEoverP = &photrail_MCconvEoverP;
std::vector<float> photrail_MCconvCotanTheta;
std::vector<float> *Pphotrail_MCconvCotanTheta = &photrail_MCconvCotanTheta;
std::vector<float> photrail_MCconvMass;
std::vector<float> *Pphotrail_MCconvMass = &photrail_MCconvMass;
std::vector<float> photrail_MCconvVertexX;
std::vector<float> *Pphotrail_MCconvVertexX = &photrail_MCconvVertexX;
std::vector<float> photrail_MCconvVertexY;
std::vector<float> *Pphotrail_MCconvVertexY = &photrail_MCconvVertexY;
std::vector<float> photrail_MCconvVertexZ;
std::vector<float> *Pphotrail_MCconvVertexZ = &photrail_MCconvVertexZ;
std::vector<float> photrail_eleMCtruthBrem;
std::vector<float> *Pphotrail_eleMCtruthBrem = &photrail_eleMCtruthBrem;
std::vector<int> photrail_eleMCtruthNBrem;
std::vector<int> *Pphotrail_eleMCtruthNBrem = &photrail_eleMCtruthNBrem;



// kinevars
std::vector<float> photrail_et;
std::vector<float> *Pphotrail_et = &photrail_et;
std::vector<float> photrail_eta;
std::vector<float> *Pphotrail_eta = &photrail_eta;
std::vector<float> photrail_SCeta;
std::vector<float> *Pphotrail_SCeta = &photrail_SCeta;



// cluster shape vars
std::vector<float> photrail_r9;
std::vector<float> *Pphotrail_r9 = &photrail_r9;
std::vector<float> photrail_cPP;
std::vector<float> *Pphotrail_cPP = &photrail_cPP;
std::vector<float> photrail_cEP;
std::vector<float> *Pphotrail_cEP = &photrail_cEP;
std::vector<float> photrail_cEE;
std::vector<float> *Pphotrail_cEE = &photrail_cEE;
std::vector<float> photrail_r19;
std::vector<float> *Pphotrail_r19 = &photrail_r19;
std::vector<float> photrail_SCEraw;
std::vector<float> *Pphotrail_SCEraw = &photrail_SCEraw;
std::vector<float> photrail_eMax;
std::vector<float> *Pphotrail_eMax = &photrail_eMax;
std::vector<float> photrail_e2x2;
std::vector<float> *Pphotrail_e2x2 = &photrail_e2x2;
std::vector<float> photrail_e5x5;
std::vector<float> *Pphotrail_e5x5 = &photrail_e5x5;
std::vector<float> photrail_ratioSeed;
std::vector<float> *Pphotrail_ratioSeed = &photrail_ratioSeed;
std::vector<float> photrail_ratioS4;
std::vector<float> *Pphotrail_ratioS4 = &photrail_ratioS4;
std::vector<float> photrail_lambdaRatio;
std::vector<float> *Pphotrail_lambdaRatio = &photrail_lambdaRatio;
std::vector<float> photrail_lamdbaDivCov;
std::vector<float> *Pphotrail_lamdbaDivCov = &photrail_lamdbaDivCov;
std::vector<float> photrail_secondMomentMaj;
std::vector<float> *Pphotrail_secondMomentMaj = &photrail_secondMomentMaj;
std::vector<float> photrail_secondMomentMin;
std::vector<float> *Pphotrail_secondMomentMin = &photrail_secondMomentMin;
std::vector<float> photrail_secondMomentAlpha;
std::vector<float> *Pphotrail_secondMomentAlpha = &photrail_secondMomentAlpha;
std::vector<float> photrail_covAngle;
std::vector<float> *Pphotrail_covAngle = &photrail_covAngle;
std::vector<float> photrail_covAngle2;
std::vector<float> *Pphotrail_covAngle2 = &photrail_covAngle2;
std::vector<float> photrail_S9overS9minusS1S2;
std::vector<float> *Pphotrail_S9overS9minusS1S2 = &photrail_S9overS9minusS1S2;
std::vector<float> photrail_etawidth;
std::vector<float> *Pphotrail_etawidth = &photrail_etawidth;
std::vector<float> photrail_phiwidth;
std::vector<float> *Pphotrail_phiwidth = &photrail_phiwidth;
std::vector<float> photrail_sigieta;
std::vector<float> *Pphotrail_sigieta = &photrail_sigieta;
std::vector<float> photrail_SCbr;
std::vector<float> *Pphotrail_SCbr = &photrail_SCbr;
std::vector<float> photrail_seedEnergy;
std::vector<float> *Pphotrail_seedEnergy = &photrail_seedEnergy;
std::vector<float> photrail_seedTime;
std::vector<float> *Pphotrail_seedTime = &photrail_seedTime;
std::vector<float> photrail_SCphi;
std::vector<float> *Pphotrail_SCphi = &photrail_SCphi;
std::vector<float> photrail_SCEtraw;
std::vector<float> *Pphotrail_SCEtraw = &photrail_SCEtraw;
std::vector<float> photrail_SCEt;
std::vector<float> *Pphotrail_SCEt = &photrail_SCEt;
std::vector<float> photrail_SCr9;
std::vector<float> *Pphotrail_SCr9 = &photrail_SCr9;
std::vector<int> photrail_SCnbBC;
std::vector<int> *Pphotrail_SCnbBC = &photrail_SCnbBC;
std::vector<int> photrail_SCnXtal;
std::vector<int> *Pphotrail_SCnXtal = &photrail_SCnXtal;
std::vector<int> photrail_isMatchingWithHLTObject;
std::vector<int> *Pphotrail_isMatchingWithHLTObject = &photrail_isMatchingWithHLTObject;
std::vector<int> photrail_isConverted;
std::vector<int> *Pphotrail_isConverted = &photrail_isConverted;
std::vector<int> photrail_NtrackConv;
std::vector<int> *Pphotrail_NtrackConv = &photrail_NtrackConv;
std::vector<float> photrail_convEoverP;
std::vector<float> *Pphotrail_convEoverP = &photrail_convEoverP;
std::vector<float> photrail_convMass;
std::vector<float> *Pphotrail_convMass = &photrail_convMass;
std::vector<float> photrail_convCotanTheta;
std::vector<float> *Pphotrail_convCotanTheta = &photrail_convCotanTheta;
std::vector<float> photrail_convLikely;
std::vector<float> *Pphotrail_convLikely = &photrail_convLikely;
std::vector<float> photrail_convVertexX;
std::vector<float> *Pphotrail_convVertexX = &photrail_convVertexX;
std::vector<float> photrail_convVertexY;
std::vector<float> *Pphotrail_convVertexY = &photrail_convVertexY;
std::vector<float> photrail_convVertexZ;
std::vector<float> *Pphotrail_convVertexZ = &photrail_convVertexZ;
std::vector<float> photrail_fBrem;
std::vector<float> *Pphotrail_fBrem = &photrail_fBrem;
std::vector<int> photrail_isAspike;
std::vector<int> *Pphotrail_isAspike = &photrail_isAspike;
std::vector<float> photrail_xVertex ;
std::vector<float> *Pphotrail_xVertex  = &photrail_xVertex ;
std::vector<float> photrail_yVertex ;
std::vector<float> *Pphotrail_yVertex  = &photrail_yVertex ;
std::vector<float> photrail_zVertex ;
std::vector<float> *Pphotrail_zVertex  = &photrail_zVertex ;




//isolation varia>bles
std::vector<float> photrail_HcalIso;
std::vector<float> *Pphotrail_HcalIso = &photrail_HcalIso;
std::vector<float> photrail_EcalIso;
std::vector<float> *Pphotrail_EcalIso = &photrail_EcalIso;
std::vector<float> photrail_TrackerIso;
std::vector<float> *Pphotrail_TrackerIso = &photrail_TrackerIso;
std::vector<float> photrail_HcalIso_MIT03;
std::vector<float> *Pphotrail_HcalIso_MIT03 = &photrail_HcalIso_MIT03;
std::vector<float> photrail_EcalIso_MIT03;
std::vector<float> *Pphotrail_EcalIso_MIT03 = &photrail_EcalIso_MIT03;
std::vector<float> photrail_TrackerIso_MIT03;
std::vector<float> *Pphotrail_TrackerIso_MIT03 = &photrail_TrackerIso_MIT03;
std::vector<float> photrail_HcalEcal_MIT03;
std::vector<float> *Pphotrail_HcalEcal_MIT03 = &photrail_HcalEcal_MIT03;
std::vector<float> photrail_AbsTrackerIso;
std::vector<float> *Pphotrail_AbsTrackerIso = &photrail_AbsTrackerIso;
std::vector<float> photrail_HcalIsodR03;
std::vector<float> *Pphotrail_HcalIsodR03 = &photrail_HcalIsodR03;
std::vector<float> photrail_EcalIsodR03;
std::vector<float> *Pphotrail_EcalIsodR03 = &photrail_EcalIsodR03;
std::vector<float> photrail_TrackerIsodR03;
std::vector<float> *Pphotrail_TrackerIsodR03 = &photrail_TrackerIsodR03;
std::vector<float> photrail_HcalIsoPerso;
std::vector<float> *Pphotrail_HcalIsoPerso = &photrail_HcalIsoPerso;
std::vector<float> photrail_EcalIsoPerso;
std::vector<float> *Pphotrail_EcalIsoPerso = &photrail_EcalIsoPerso;
std::vector<float> photrail_hoe;
std::vector<float> *Pphotrail_hoe = &photrail_hoe;


//other variables
std::vector<int> photrail_Cat;
std::vector<int> *Pphotrail_Cat = &photrail_Cat;
std::vector<int> photrail_HasPixSeed;
std::vector<int> *Pphotrail_HasPixSeed = &photrail_HasPixSeed;
std::vector<int> photrail_seedSeverity;
std::vector<int> *Pphotrail_seedSeverity = &photrail_seedSeverity;
std::vector<int> photrail_recoFlag;
std::vector<int> *Pphotrail_recoFlag = &photrail_recoFlag;
std::vector<int> photrail_isEB;
std::vector<int> *Pphotrail_isEB = &photrail_isEB;
std::vector<int> photrail_isEE;
std::vector<int> *Pphotrail_isEE = &photrail_isEE;


TBranch* event_br = 0;
TBranch* run_br = 0;
TBranch* mcParticles_br = 0;
TBranch* genJets_br = 0;
TBranch* genMETs_br = 0;
TBranch* mcSignalMuMuGamma_br = 0;
TBranch* mcTopTopEvent_br = 0;
TBranch* mcPhotons_br = 0;
TBranch* mcElectrons_br = 0;
TBranch* beamSpot_br = 0;
TBranch* vertices_br = 0;
TBranch* zeeVertices_br = 0;
TBranch* tracks_br = 0;
TBranch* jets_br = 0;
TBranch* pflowjets_br = 0;
TBranch* sisconejets_br = 0;
TBranch* muons_br = 0;
TBranch* electrons_br = 0;
TBranch* photons_br = 0;
TBranch* clusters_br = 0;
TBranch* superClusters_br = 0;
TBranch* conversions_br = 0;
TBranch* met_br = 0;
TBranch* bardak_br = 0;
TBranch* HLTObjects_br = 0;

TRootEvent* event = 0;
TRootRun* runInfos = 0;
TClonesArray* mcParticles; // = new TClonesArray("TRootMCParticle", 0);
TClonesArray* genJets; // = new TClonesArray("TRootParticle", 0);
TClonesArray* genMETs; // = new TClonesArray("TRootParticle", 0);
TRootSignalEvent* mcMuMuGammaEvent = 0;
TRootSignalEvent* mcTopTopEvent = 0;
TClonesArray* mcPhotons; // = new TClonesArray("TRootMCPhoton", 0);
TClonesArray* mcElectrons;
TRootBeamSpot* beamSpot = 0;
TClonesArray* vertices; // = new TClonesArray("TRootVertex", 0);
TClonesArray* zeeVertices; // = new TClonesArray("TRootVertex", 0);
TClonesArray* tracks; // = new TClonesArray("TRootTrack", 0);
TClonesArray* jets; // = new TClonesArray("TRootJet", 0);
TClonesArray* pflowjets; // = new TClonesArray("TRootJet", 0);
TClonesArray* sisconejets; // = new TClonesArray("TRootJet", 0);
TClonesArray* muons; /// = new TClonesArray("TRootMuon", 0);
TClonesArray* electrons; // = new TClonesArray("TRootElectron", 0);
TClonesArray* photons; // = new TClonesArray("TRootPhoton", 0);
TClonesArray* clusters; // = new TClonesArray("TRootCluster", 0);
TClonesArray* superClusters; // = new TClonesArray("TRootSuperCluster", 0);
TClonesArray* conversionTracks; // = new TClonesArray("TRootTrack", 0);
TClonesArray* met; // = new TClonesArray("TRootMET", 0);
TRootBardak* bardak = 0;
TRootSuperCluster* mysc;
TClonesArray* HLTObjects;
