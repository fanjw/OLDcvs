#include <TH1D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <TPostScript.h>
#include <TCut.h>

//////////////////////////
using namespace std; 
using namespace RooFit;


//const static int BinTotal= 160; //80; //1; //32;//80;
//const static double BinXLow=100; //80; //400.0; //100.;
//const static double BinXHig=180; //180.0;  //180; //700.0;//160.;
/*
const static double Kpp=1.3;
const static double KppMad=1.15;
const static double Kpf=1.3;
const static double Kff=1.0;
const static double Kdy=1.15;
*/
const static double Kpp=1.0;
const static double KppMad=1.0;
const static double Kpf=1.0; 
const static double Kff=1.5;  //MITpre Kff=1.5
const static double Kdy=1.0;

/*
const static TString TreeName="DiPhotonTree";
const static int NPlots=21;
const static TString Object[NPlots] = { 
  "Diphoton_Mass", "Diphoton_PT", "Diphoton_Eta",
  "Diphoton_ETmaxOverMass","Diphoton_ETminOverMass","Diphoton_NormalizedPL","Diphoton_PhotonDeltaEta","Diphoton_SumPToverJetPT",
  "Diphoton_LeadMVAOutput", "Diphoton_SubLeadMVAOutput",
  "Diphoton_LeadTaoMVAOutput", "Diphoton_SubLeadTaoMVAOutput",
  "Diphoton_CosThetaStar","Diphoton_AlphaL/Diphoton_AlphaT","Diphoton_CosHiggsJet",
  "Diphoton_dmom", "Diphoton_dmom_wrongvtx", "Diphoton_CosDeltaPhi", "Diphoton_VtxProb", "Diphoton_MVAOutput","Diphoton_MVAOutputTao"
};
const static TString HisTitle[NPlots] = { 
    "M_{#gamma#gamma}", "pT_{#gamma#gamma}", "#eta_{#gamma#gamma}",
    "E_{Tmax}^{#gamma}/M_{#gamma#gamma}","E_{Tmin}^{#gamma}/M_{#gamma#gamma}", "P_{L}^{#gamma#gamma}/1000","#Delta#eta_{#gamma#gamma}","(E_{T}^{#gamma1}+E_{T}^{#gamma2})/(P_{T}^{jet1}+P_{T}^{jet2})",
    "MVA_{#gamma1}", "MVA_{#gamma2}",
    "Tao MVA_{#gamma1}", "Tao MVA_{#gamma2}",
    "cos#theta^{*}_{#gamma#gamma}", "#alpha_{L}/#alpha_{T}","cos(Higgs,Jet)",
    "#sigmam/m_{#gamma#gamma} (good vertex)", "#sigmam/m_{#gamma#gamma} (worst vertex)", "cos(#Delta#phi(#gamma#gamma))", "Vertex probability", "event MVA output", "event MVA output (My PhotonID MVA)"                                       
}; 
const static int BinTotal[NPlots] =   { 
  40,  50, 60,
  100, 100, 100, 100, 100,
  100, 100,
  100, 100,
  100, 100, 100,
  60, 60, 100, 80, 50, 50
};
const static double BinXLow[NPlots] = { 
  100, 50., -6.,
  0.2, 0.2, 0., 0., 0.5,
  -0.6, -0.6,
  -0.6, -0.6,
  0.,  0., 0.0,
  0., 0.,  -1.0, 0.2, -1.0, -1.0
};
const static double BinXHig[NPlots] = { 
  180, 250., 6.0,
  1.2, 0.7, 1., 3.0, 1.1,
  0.4, 0.4,
  0.4, 0.4,
  1.0, 10.0, 1.0,
  0.06, 0.06, 1.0, 1.0, 1.0, 1.0
};
//========================
//TCut * TaoCut=new TCut( Form("event_DiPhotonCiC4Level>=1  && Diphoton_Mass>100.  && Diphoton_Mass<180.") );
TCut * TaoCut=new TCut( Form("event_DiPhotonMIT>=0  && Diphoton_Mass>100. && Diphoton_Mass<180.") );
//TCut * TaoCutMC=new TCut(  Form("((event_weight>1e-9 && event_weight<1e3)?event_weight:0.0)*(event_DiPhotonCiC4Level>=1 && Diphoton_Mass>100.  && Diphoton_Mass<180.)") );
TCut * TaoCutMC=new TCut(  Form("((event_weight>1e-9 && event_weight<1e3)?event_weight:0.0)*(event_DiPhotonMIT>=0 && Diphoton_Mass>100. && Diphoton_Mass<180.)") );
*/
///////////////////////////////////////////////////////////////////

const static TString TreeName="PhotonTree";
const static int NPlots=30; //16;
const static TString Object[NPlots] = { 
  "myphoton_pt", "myphoton_r9", "myphoton_sieie","myphoton_r19", "myphoton_cep",
  "myphoton_SCeta","myphoton_phiwidth", "event_Nvtx",
  "myphoton_s4ratio"," myphoton_lambdaratio", "myphoton_etawidth", "myphoton_brem","myphoton_maxoraw",
  "myphoton_lambdadivcov", "myphoton_smaj", "myphoton_MVAOutput", //"myphoton_NNclustershape"
   "myphoton_hoe", 
  "myphoton_isoOverEt", "myphoton_badisoOverEt", "myphoton_trkisooet", "myphoton_drtotk",
  "myphoton_trkisohollowdr03","myphoton_ecalisodr03","myphoton_hcalisodr03", 
  "myphoton_pfsumisogood03", "myphoton_pfsumisobad03", "myphoton_pfneutraliso03", "myphoton_pfphotoniso03", "myphoton_pfchargedisogood03", "myphoton_pfchargedisobad03"
}; 
const static TString HisTitle[NPlots] = { 
  "pT_{#gamma}", "R_{9}", "#sigma_{i#etai#eta}", "S_{1}/S_{9}", "Covariance(#eta,#phi)",
  "#eta_{#gamma SC}","#phi_{Width}", "Number of vertices",
  "S_{2x2Max}/S_{5x5}", "(#sigma_{#eta#eta}+#sigma_{#phi#phi}-#lambda)/(#sigma_{#eta#eta}+#sigma_{#phi#phi}+#lambda)","#eta_{Width}", "#phi_{Width}/#eta_{Width}", "E_{seed}/E_{SCraw}",
  "#lambda^{-}/C_{#eta#eta}","M_{MAJ}^{2}", "O_{MVA}", "H/E",
  "#sumET/Et^{#gamma} (selected vertex)", "#sumET/Et^{#gamma} (worst vertex)", "#sumPT^{trk}/Et^{#gamma} (selected vertex)", "#DeltaR_{eSc,eTrk}", 
  "#sumpT (track, #DeltaR=0.3)", "#sumET (ECAL, #DeltaR=0.3)", "#sumET (HCAL, #DeltaR=0.3)",
  "pf #sumET (#DeltaR=0.3, selected vertex)", "pf #sumET (#DeltaR=0.3, worst vertex)", "pf #sumET (Neutral, #DeltaR=0.3)","pf #sumET (#gamma, #DeltaR=0.3)", "pf #sumET (Charged,#DeltaR=0.3, selected vertex)", "pf #sumET (Charged, #DeltaR=0.3, worst vertex)"
};
const static int BinTotal[NPlots] =   { 
  100, 100, 100, 50, 80,
  50,  150, 50,
  60, 50, 60, 50, 50,
  //100, 50, 120, 50, 50,
  115, 50, 60, 100,
  100, 100, 100, 100, 
  100, 100, 100, 
  100, 100, 100, 100, 100, 100
};
const static double BinXLow[NPlots] = { 
  15.0, 0.0, 0.001, 0.0, -0.0004,
  -2.5, 0.0, -0.5, 
  0.4, 0.0, 0.0, 0.0, 0.0,
  0.25, 0.0, -0.9, 0.0,
  -1.0,-1.0,-1.0,-0.1,
  -0.1,-0.1,-0.1,
  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1
  //0.25, -1.05, 0.0
};


const static double BinXHig[NPlots] = { 
  215.0, 1.0, 0.041, 1.0, 0.0004,
  2.5,  0.15, 49.5,  
  1.0,  1.0, 0.03, 10.0, 1.0,
  //1.0,  1.0, 0.06, 10.0, 1.0, 
  2.05,  3.0, 0.3,   0.2,
  19.9,  19.9, 19.9,  0.9,
  19.9,  19.9, 19.9,
  19.9,  19.9, 19.9, 19.9,  19.9, 19.9
  //1.05,  1.95, 1.0   
};

//========================

//TCut * TaoCut=new TCut("event_DiPhotonMIT>=0 && myphoton_IndexInDiphoton>=0 && event_DiPhotonMass>100. && event_DiPhotonMass<180."); 
//TCut * TaoCutMC=new TCut(" ((event_weight>1e-9 && event_weight<1e2)?event_weight:0.0)  * (event_DiPhotonMIT>=0 && myphoton_IndexInDiphoton>=0 && event_DiPhotonMass>100. && event_DiPhotonMass<180.) ");


TCut * TaoCutEB=new TCut("event_DiPhotonMIT>=0  && myphoton_IndexInDiphoton>=0 && event_DiPhotonMass>100. && event_DiPhotonMass<180. && myphoton_catind<=1");
TCut * TaoCutEE=new TCut("event_DiPhotonMIT>=0  && myphoton_IndexInDiphoton>=0 && event_DiPhotonMass>100. && event_DiPhotonMass<180. && myphoton_catind>=2");
TCut * TaoCutMCEB=new TCut(" ((event_weight>1e-9 && event_weight<1e2)?event_weight:0.0)  * (event_DiPhotonMIT>=0  && myphoton_IndexInDiphoton>=0 && event_DiPhotonMass>100. && event_DiPhotonMass<180. && myphoton_catind<=1) ");
TCut * TaoCutMCEE=new TCut(" ((event_weight>1e-9 && event_weight<1e2)?event_weight:0.0)  * (event_DiPhotonMIT>=0  && myphoton_IndexInDiphoton>=0 && event_DiPhotonMass>100. && event_DiPhotonMass<180. && myphoton_catind>=2) ");

//---------------------------------------------
//const static double RecorDataInteLumi=1.146; //fb-1
const static double RecorDataInteLumi=4.763;
///======Summer11 MC samples===================
//////////////////////////
void Fun_do_NLLFit(int part, TFile *myFile) {


  TChain *Data_Tree=new TChain(TreeName);
  Data_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/PhotonRun2011A_Clean/*.root");
  Data_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/PhotonRun2011B_Clean/*.root");


  //===Gamma Born Jet=====
  TChain *dpj_Tree=new TChain(TreeName);
  dpj_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/DiPhotonJets/*.root");
   
  //===Gamma Box===== 
  TChain *box10_Tree=new TChain(TreeName);
  TChain *box25_Tree=new TChain(TreeName);
  TChain *box250_Tree=new TChain(TreeName);
  box10_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/BoxPt10to25/*.root");
  box25_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/BoxPt25to250/*.root");
  box250_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/BoxPt250/*.root");

  //===Gamma Jet====
  TChain *gjpp_Tree=new TChain(TreeName);
  gjpp_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/GJet_Pt-20_pp/*.root");

  TChain *gjpf_Tree=new TChain(TreeName);
  gjpf_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/GJet_Pt-20_pf/*.root");

   //===QCD EM============
  TChain *qcd30pp_Tree=new TChain(TreeName);
  qcd30pp_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/QCDPt30to40_pp/*.root");

  TChain *qcd30pf_Tree=new TChain(TreeName);
  qcd30pf_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/QCDPt30to40_pf/*.root");

  TChain *qcd30ff_Tree=new TChain(TreeName);
  qcd30ff_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/QCDPt30to40_ff/*.root");

  TChain *qcd40pp_Tree=new TChain(TreeName);
  qcd40pp_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/QCDPt40_pp/*.root");

  TChain *qcd40pf_Tree=new TChain(TreeName);
  qcd40pf_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/QCDPt40_pf/*.root");

  TChain *qcd40ff_Tree=new TChain(TreeName);
  qcd40ff_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/QCDPt40_ff/*.root");

  //=====DY ee==========
  TChain *dy_Tree=new TChain(TreeName);
  dy_Tree->Add("/farm/t1/taojq/HggAna/h2gglobe_V11_04_05_reduction_mva/DYJetsToLL_M50/*.root");
  //=================================
 
 


 
  //for(int i=0; i<NPlots; i++){
  for(int i=8; i<9; i++){
  //for(int i=NPlots-2; i<NPlots; i++){
  //for(int i=0; i<4; i++){
  //for(int i=NPlots-1; i<NPlots; i++){

        TString nomDATA;
        if(part == 0)  nomDATA = Object[i] + "_EB_DATA";
        else if (part == 1) nomDATA = Object[i] + "_EE_DATA";
        
        TString nomMC;
        if(part == 0)  nomMC = Object[i] + "_EB_MC";
        else if (part == 1) nomMC = Object[i] + "_EE_MC";
        
        
        cout << "on va plot     "<< nomDATA <<endl;
        TH1F *hdata = new TH1F(nomDATA, "", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_Data = Object[i] + ">>" + nomDATA;
        if(part == 0)  Data_Tree->Draw(variable_Data, *TaoCutEB);
        else if(part == 1)  Data_Tree->Draw(variable_Data, *TaoCutEE);

        double Ntot=double(hdata->Integral()); double rms=hdata->GetRMS();
        cout<<"Taojq: Ndata  = "<<Ntot<<" with RMS="<<rms<<" with Inte Lumi "<<RecorDataInteLumi<<endl;
      
       


        //RooRealVar x("x","x", BinXLow[i], BinXHig[i]);
        RooRealVar x("x","x", 0.8, 0.96);

        RooDataHist dataHist("dataHist","the data histo ",x, hdata);
        RooHistPdf theHistoPdf("theHistoPdf","the Histo pdf",x,dataHist);
        RooAbsReal*  nll;
        
        //TH1D *hLL = new TH1D("h","", 30, 0.895, 1.195);
        TH1D *hLL = new TH1D("h","", 20, 0.9895, 1.0095);
        float minNLL = 100000;
        float theGoodI;


        cout << "on va plot     "<< nomMC <<endl;
    for (int j = -10 ; j < 10 ; j++){
                  
        TString VariableToDraw = Form("%f*", 1.0+0.001*j) +Object[i];
        cout << VariableToDraw <<endl;



        TH1D *hmc=new TH1D(nomMC,"",BinTotal[i],BinXLow[i],BinXHig[i]);
         //===PP Madgraph: Gamma Born Jets=====
        TString hpp_dpj_tempNom;
        if (part ==0 )  hpp_dpj_tempNom = "hpp_dpj_temp_EB";
        else if (part ==1 )   hpp_dpj_tempNom = "hpp_dpj_temp_EE"; 
        TH1F *hpp_dpj_temp = new TH1F(hpp_dpj_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_dpj = VariableToDraw + ">>" + hpp_dpj_tempNom;
        if (part == 0 ) dpj_Tree->Draw(variable_dpj, *TaoCutMCEB);
        else if (part == 1) dpj_Tree->Draw(variable_dpj, *TaoCutMCEE);
        TH1F *hpp = (TH1F*)gDirectory->Get(hpp_dpj_tempNom);
        c1->Clear();
        hpp->Scale(KppMad);

        //===PP  Pythia:  Gamma Box=====
        TString hpp_box10_tempNom;
        if (part ==0 ) hpp_box10_tempNom = "hpp_box10_temp_EB";
        else if (part ==1 ) hpp_box10_tempNom = "hpp_box10_temp_EE";
        TH1F *hpp_box10_temp = new TH1F(hpp_box10_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_box10 = VariableToDraw + ">>" + hpp_box10_tempNom;
        if (part == 0)  box10_Tree->Draw(variable_box10, *TaoCutMCEB);
        else if (part ==1 )  box10_Tree->Draw(variable_box10, *TaoCutMCEE);
        TH1F *hpp_box10 = (TH1F*)gDirectory->Get(hpp_box10_tempNom);
        c1->Clear();

        TString hpp_box25_tempNom;
        if (part ==0 ) hpp_box25_tempNom = "hpp_box25_temp_EB"
        else if (part ==1 ) hpp_box25_tempNom = "hpp_box25_temp_EE";
        TH1F *hpp_box25_temp = new TH1F(hpp_box25_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_box25 = VariableToDraw + ">>" + hpp_box25_tempNom;
        if (part == 0)  box25_Tree->Draw(variable_box25, *TaoCutMCEB);
        else if (part == 1)box25_Tree->Draw(variable_box25, *TaoCutMCEE);
        TH1F *hpp_box25 = (TH1F*)gDirectory->Get(hpp_box25_tempNom);
        c1->Clear();

        TString hpp_box250_tempNom;
        if (part ==0 ) hpp_box250_tempNom = "hpp_box250_temp_EB";
        else if (part ==1 ) hpp_box250_tempNom = "hpp_box250_temp_EE";
        TH1F *hpp_box250_temp = new TH1F(hpp_box250_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_box250 = VariableToDraw + ">>" + hpp_box250_tempNom;
        if (part == 0) box250_Tree->Draw(variable_box250, *TaoCutMCEB);
        else if (part ==1 ) box250_Tree->Draw(variable_box250, *TaoCutMCEE);
        TH1F *hpp_box250 = (TH1F*)gDirectory->Get(hpp_box250_tempNom);
        c1->Clear();

        hpp_box10->Scale(Kpp);
        hpp_box25->Scale(Kpp);
        hpp_box250->Scale(Kpp);
        hpp->Add(hpp_box10, 1.0);
        hpp->Add(hpp_box25, 1.0);
        hpp->Add(hpp_box250, 1.0);

         //===PP Pythia: Gamma Jet=====
        TString hpp_gjpp_tempNom;
        if (part ==0 ) hpp_gjpp_tempNom = "hpp_gjpp_temp_EB";
        else if (part == 1) hpp_gjpp_tempNom = "hpp_gjpp_temp_EE";
        TH1F *hpp_gjpp_temp = new TH1F(hpp_gjpp_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_gjpp = VariableToDraw + ">>" + hpp_gjpp_tempNom;
        if (part == 0)  gjpp_Tree->Draw(variable_gjpp, *TaoCutMCEB);
        else if (part == 1)  gjpp_Tree->Draw(variable_gjpp, *TaoCutMCEE);
        TH1F *hpp_gjpp = (TH1F*)gDirectory->Get(hpp_gjpp_tempNom);
        c1->Clear();
        hpp_gjpp->Scale(Kpp);
        hpp->Add(hpp_gjpp, 1.0);

         //===PP Pythia: QCD EM=====
        TString hpp_qcd30pp_tempNom;
        if (part ==0 ) hpp_qcd30pp_tempNom = "hpp_qcd30pp_temp_EB";
        else if (part == 1)  hpp_qcd30pp_tempNom = "hpp_qcd30pp_temp_EE";
        TH1F *hpp_qcd30pp_temp = new TH1F(hpp_qcd30pp_tempNom, "", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_qcd30pp = VariableToDraw + ">>" + hpp_qcd30pp_tempNom;
        if (part == 0)  qcd30pp_Tree->Draw(variable_qcd30pp, *TaoCutMCEB);
        else if (part ==1 ) qcd30pp_Tree->Draw(variable_qcd30pp, *TaoCutMCEE);
        TH1F *hpp_qcd30pp = (TH1F*)gDirectory->Get(hpp_qcd30pp_tempNom);
        c1->Clear();

        TString hpp_qcd40pp_tempNom;
        if (part ==0 )  hpp_qcd40pp_tempNom = "hpp_qcd40pp_temp_EB";
        else if (part == 1) hpp_qcd40pp_tempNom = "hpp_qcd40pp_temp_EE";
        TH1F *hpp_qcd40pp_temp = new TH1F(hpp_qcd40pp_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_qcd40pp = VariableToDraw + ">>" + hpp_qcd40pp_tempNom;
        if (part == 0)  qcd40pp_Tree->Draw(variable_qcd40pp, *TaoCutMCEB);
        else if (part == 1 )qcd40pp_Tree->Draw(variable_qcd40pp, *TaoCutMCEE);
        TH1F *hpp_qcd40pp = (TH1F*)gDirectory->Get(hpp_qcd40pp_tempNom);
        c1->Clear();

        hpp_qcd30pp->Scale(Kpp);
        hpp_qcd40pp->Scale(Kpp);
        hpp->Add(hpp_qcd30pp, 1.0);
        hpp->Add(hpp_qcd40pp, 1.0);
     

         //===PF Pythia: Gamma Jet=====
        TString hpf_gjpf_tempNom;
        if (part ==0 )   hpf_gjpf_tempNom= "hpf_gjpf_temp_EB";
        else if (part == 1) hpf_gjpf_tempNom= "hpf_gjpf_temp_EE";
        TH1F *hpf_gjpf_temp = new TH1F(hpf_gjpf_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_gjpf = VariableToDraw + ">>" + hpf_gjpf_tempNom;
        if (part == 0)  gjpf_Tree->Draw(variable_gjpf, *TaoCutMCEB);
        else if (part == 1 )  gjpf_Tree->Draw(variable_gjpf, *TaoCutMCEE);
        TH1F *hpf = (TH1F*)gDirectory->Get(hpf_gjpf_tempNom);
        c1->Clear();
        hpf->Scale(Kpf);

        //===PF Pythia: QCD EM=====
        TString hpf_qcd30pf_tempNom;
        if (part ==0 ) hpf_qcd30pf_tempNom = "hpf_qcd30pf_temp_EB";
        else if (part == 1) hpf_qcd30pf_tempNom = "hpf_qcd30pf_temp_EE"; 
        TH1F *hpf_qcd30pf_temp = new TH1F(hpf_qcd30pf_tempNom,"",BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_qcd30pf = VariableToDraw + ">>" + hpf_qcd30pf_tempNom;
        if (part == 0 )   qcd30pf_Tree->Draw(variable_qcd30pf, *TaoCutMCEB);
        else if (part == 1 )  qcd30pf_Tree->Draw(variable_qcd30pf, *TaoCutMCEE);
        TH1F *hpf_qcd30pf = (TH1F*)gDirectory->Get(hpf_qcd30pf_tempNom);
        c1->Clear();

        TString hpf_qcd40pf_tempNom;
        if (part ==0 )  hpf_qcd40pf_tempNom= "hpf_qcd40pf_temp_EB";
        else if (part == 1) hpf_qcd40pf_tempNom= "hpf_qcd40pf_temp_EE";
        TH1F *hpf_qcd40pf_temp = new TH1F(hpf_qcd40pf_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_qcd40pf = VariableToDraw + ">>" + hpf_qcd40pf_tempNom;
        if (part == 0)  qcd40pf_Tree->Draw(variable_qcd40pf, *TaoCutMCEB);
        else if (part ==1 )  qcd40pf_Tree->Draw(variable_qcd40pf, *TaoCutMCEE);
        TH1F *hpf_qcd40pf = (TH1F*)gDirectory->Get(hpf_qcd40pf_tempNom);
        c1->Clear();

        hpf_qcd30pf->Scale(Kpf);
        hpf_qcd40pf->Scale(Kpf);
        hpf->Add(hpf_qcd30pf, 1.0);
        hpf->Add(hpf_qcd40pf, 1.0);
     
 
         //===FF Pythia: QCD EM=====
        TString hpp_qcd30ff_tempNom;
        if (part ==0 )  hpp_qcd30ff_tempNom = "hpp_qcd30ff_temp_EB";
        else if (part == 1 ) hpp_qcd30ff_tempNom = "hpp_qcd30ff_temp_EE";
        TH1F *hpp_qcd30ff_temp = new TH1F(hpp_qcd30ff_tempNom,"",BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_qcd30ff = VariableToDraw + ">>" + hpp_qcd30ff_tempNom;
        if (part ==0 )   qcd30ff_Tree->Draw(variable_qcd30ff, *TaoCutMCEB);
        else if (part ==1 )   qcd30ff_Tree->Draw(variable_qcd30ff, *TaoCutMCEE);
        TH1F *hff = (TH1F*)gDirectory->Get(hpp_qcd30ff_tempNom);
        c1->Clear();

        TString hpp_qcd40ff_tempNom;
        if (part ==0 ) hpp_qcd40ff_tempNom = "hpp_qcd40ff_temp_EB";
        else if (part == 1 ) hpp_qcd40ff_tempNom = "hpp_qcd40ff_temp_EE";
        TH1F *hpp_qcd40ff_temp = new TH1F(hpp_qcd40ff_tempNom,"", BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_qcd40ff = VariableToDraw + ">>" + hpp_qcd40ff_tempNom;
        if (part == 0 )   qcd40ff_Tree->Draw(variable_qcd40ff, *TaoCutMCEB);
        else if (part == 1)   qcd40ff_Tree->Draw(variable_qcd40ff, *TaoCutMCEE);
        TH1F *hff_qcd40ff = (TH1F*)gDirectory->Get(hpp_qcd40ff_tempNom);
        c1->Clear();

        hff->Scale(Kff);
        hff_qcd40ff->Scale(Kff);
        hff->Add(hff_qcd40ff, 1.0);
      

        //=====DY ee=====
        TString hdy_tempNom;
        if (part ==0 )  hdy_tempNom = "hdy_temp_EB";
        else if (part == 1 )  hdy_tempNom = "hdy_temp_EE"; 
        TH1F *hdy_temp = new TH1F(hdy_tempNom,"",BinTotal[i],BinXLow[i],BinXHig[i]);
        TString variable_dy = VariableToDraw + ">>" + hdy_tempNom;
        if (part == 0)  dy_Tree->Draw(variable_dy, *TaoCutMCEB);
        else if (part == 1)  dy_Tree->Draw(variable_dy, *TaoCutMCEE);
        TH1F *hdy = (TH1F*)gDirectory->Get(hdy_tempNom);
        c1->Clear();
        hdy->Scale(Kdy);


        hmc->Add(hpp,1.0);
        hmc->Add(hpf,1.0);
        hmc->Add(hff,1.0);
        hmc->Add(hdy,1.0);



        RooDataHist mcHist("mcHist","the mc histo ",x, hmc);
        nll = theHistoPdf.createNLL(mcHist);
        hLL->Fill(1.0+0.001*j, -1.0*nll->getVal());

        if (minNLL > nll->getVal()) 
                {   
                        minNLL =  nll->getVal();
                        //theGoodI = 1.0+0.01*i;
                        theGoodI = 1.0+0.001*j;
                } 

       


        cout<< "par = " << 1.0+0.001*j << "     nll = " << -1.0*nll->getVal() <<endl;


        double NTotBkgMC= double(hmc->Integral());
        cout<<"Tot bkg (MC): "<<NTotBkgMC<<endl;
        
        delete hmc;
        delete hpp_dpj_temp;
        delete hpp_box10_temp;
        delete hpp_box25_temp;
        delete hpp_box250_temp;
        delete hpp_gjpp_temp;
        delete hpp_qcd30pp_temp;
        delete hpp_qcd40pp_temp;
        delete hpf_gjpf_temp;
        delete hpf_qcd30pf_temp;
        delete hpf_qcd40pf_temp;
        delete hpp_qcd30ff_temp;
        delete hpp_qcd40ff_temp;
        delete hdy_temp;
  }
        hLL->Write();
       
    
 }
 
}



void do_NLLFit() {



       TFile *myFile   = new TFile("histo_file_NLLFit.root","update");
       Fun_do_NLLFit(0, myFile);
       Fun_do_NLLFit(1, myFile);
       myFile->Close();
       


}


    
