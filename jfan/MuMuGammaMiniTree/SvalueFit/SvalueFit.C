//#include "TROOT.h "
#include "SvalueFit.h" 


using namespace RooFit;
using namespace std;


//int main(){
void SvalueFit(int EndCaps, int r9sup, string Category){
//void SvalueFit(int EndCaps, int r9sup){


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gROOT->ProcessLine(".x setTDRStyle.C");



  //***************** initialization *********************
  int isMC = 1; //0 1
  int FitPourcentage = 70;
  double FitPourcentageD = FitPourcentage / 100.0;
  //int EndCaps = 1;//0 1 (0 = Barrel, 1 = Endcaps)
  //int r9sup = 2;//0 1 2 (0 = low r9, 1 = high r9, 2 = no r9 cuts)
  string nomFitMethode = "RooVoigtian2";
  //string Category = "OneBin";  //"Vgamma24" "Vgamma8" "OneBin"
  string SurfaceMethod = "ProfileSurface"; //"ProfileSurface" "FittedSurface"
  string variableX = "Photon_Et";
  bool phiCracks = true;
  bool etaCracks = true;

  string isMCChain = "";
  if(isMC == 0) isMCChain = "Data";
  if(isMC == 1) isMCChain = "MC";


  double MeanTab[100] = {0.0};
  double MeanErrorTab[100] = {0.0};
  double SigmaRTab[100] = {0.0};
  double SigmaLTab[100] = {0.0};
  double SigmaTab[100] = {0.0};
  double SigmaErrorTab[100] = {0.0};
  double SigmaEffTab[100] = {0.0};
  double MinVar[100] = {0.0};
  double MaxVar[100] = {0.0};
  double xValue[100] = {0.0};
  double xErrorR[100] = {0.0};
  double xErrorL[100] = {0.0};
  //Chi2
  double yValue[100] = {0};
  double xErrorRight[100] = {0};
  double xErrorLeft[100] = {0};

  double ChiSquareTab[100] = {0.0};
  double JanChiSquareTab[100] = {0.0};
  double PValueTab[100] = {0.0};
  double PValueErrorTab[100] = {0.0};
  double DegreesOfFreedomTab[100] = {0.0};

  double minLogLikelihoodTab[100] = {0.0};
  double differenceLogLikelihoodTab[100] = {0.0};

  int entriesTot = 0;
  int anti2 = 0;
  int n = 0;
  double Mean = 0;
  double Sigma = 0;

  int fitRoo = 1;
  double pourcentage = 0.78;



  double xminChi2 , xminSigmaLSigmaR , xminSigma , xminSigmaTg , xminSigmaEff , xminSigmaEffTg , xmin1overKrecoClassical , xminJanChi2 , xminPValue;
  xminChi2 = xminSigmaLSigmaR = xminSigma = xminSigmaTg = xminSigmaEff = xminSigmaEffTg = xmin1overKrecoClassical = xminJanChi2 = xminPValue = 0.0;

  if(variableX == "Photon_SC_Eta") xminSigmaLSigmaR = xminSigma = xminSigmaTg = xminSigmaEff = xminSigmaEffTg = xmin1overKrecoClassical = xminJanChi2 = xminPValue = -3.0;

  double xmaxChi2 , xmaxSigmaLSigmaR , xmaxSigma , xmaxSigmaTg , xmaxSigmaEff , xmaxSigmaEffTg , xmax1overKrecoClassical , xmaxJanChi2 , xmaxPValue;

  if(variableX == "Photon_SC_rawEt") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmax1overKrecoClassical = xmaxJanChi2 = xmaxPValue = 250.0;
  if(variableX == "Photon_Et") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmax1overKrecoClassical = xmaxJanChi2 = xmaxPValue = 250.0;
  if(variableX == "Photon_E") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmax1overKrecoClassical = xmaxJanChi2 = xmaxPValue = 1000.0;
  if(variableX == "Photon_SC_Eta") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmax1overKrecoClassical = xmaxJanChi2 = xmaxPValue = 3.0;
  if(variableX == "Photon_SC_brem") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmax1overKrecoClassical = xmaxJanChi2 = xmaxPValue = 15.0;



 
  double yminChi2 , yminSigmaLSigmaR , yminSigma , yminSigmaTg , yminSigmaEff , yminSigmaEffTg , ymin1overKrecoClassical, yminJanChi2, yminPValue;
  yminChi2 = 0.0;
  yminSigmaLSigmaR = -0.15;
  yminSigma = -15;
  yminSigmaTg = 0.0;
  yminSigmaEff = -0.15;
  yminSigmaEffTg = 0.0;
  yminJanChi2 = 0.0;
  yminPValue = 0.0;
  if(r9sup == 2) ymin1overKrecoClassical = -15.0;
  if(r9sup == 1) ymin1overKrecoClassical = -15.0;
  if(r9sup == 0) ymin1overKrecoClassical = -15.0;
 


  double ymaxChi2 , ymaxSigmaLSigmaR , ymaxSigma , ymaxSigmaTg , ymaxSigmaEff , ymaxSigmaEffTg , ymax1overKrecoClassical, ymaxJanChi2, ymaxPValue;
  ymaxChi2 = 15.0;
  ymaxSigmaLSigmaR = 0.15;
  ymaxSigma = 15;
  ymaxSigmaTg = 10.0;
  ymaxSigmaEff = 0.15;
  ymaxSigmaEffTg = 0.1;
  ymaxJanChi2 = 15.0;
  ymaxPValue = 1.0;
  if(r9sup == 2) ymax1overKrecoClassical = 15.0;
  if(r9sup == 1) ymax1overKrecoClassical = 15.0;
  if(r9sup == 0) ymax1overKrecoClassical = 15.0;
  





  double maxDistri = 0;
  double lastX = 0;
  double Entries = 0;
  double mean_value = 0;
  double mean_error = 0;
  double sigma_value = 0;
  double sigma_value_error = 0;
  Double_t sigmaEff_value = 0;
  double sigmaR_value = 0;
  double sigmaL_value = 0;
  double ChiSquare = 0;
  double params[5] = {0};
  double JanChi2 = 0;
  double DegreesOfFreedom = 0;
  double pValue = 0;
  double minLogLikelihood = 0;
  double differenceLogLikelihood = 0;

  double mean = 0;
  double rms = 0;
  double meanMC = 0.0;

  double phiCrackSize = (double)(21.5) / (double) (1290.0);
  double phiCrackPosition = (double)(TMath::Pi()) / (double)(9.0);
  double phiOffset = -(double)(10.0 * TMath::Pi()) / (double)(180.0);

  string DossierCracks = "";
  if(phiCracks == true && etaCracks == true) DossierCracks = "WithCracks";
  if(phiCracks == false && etaCracks == false) DossierCracks = "WithoutCracks";
  if(phiCracks == true && etaCracks == false) DossierCracks = "WithoutEtaCracks";
  if(phiCracks == false && etaCracks == true) DossierCracks = "WithoutPhiCracks";


 
  string EndcapsChain;
  if(variableX == "Photon_SC_rawEt") EndcapsChain = "LimitesAllEtRaw.txt";
  if(variableX == "Photon_Et" && Category == "Vgamma24") EndcapsChain = "LimitesAllPtVgamma24.txt";
  if(variableX == "Photon_Et" && Category == "Vgamma8") EndcapsChain = "LimitesAllPtVgamma8.txt";
  if(variableX == "Photon_Et" && Category == "OneBin") EndcapsChain = "LimitesAllPtOneBin.txt";
  if(variableX == "Photon_E" && EndCaps == 0) EndcapsChain = "LimitesEnergyBARREL.txt";
  if(variableX == "Photon_E" && EndCaps == 1) EndcapsChain = "LimitesEnergyENDCAPS.txt";
  if(variableX == "Photon_SC_Eta" && EndCaps == 0) EndcapsChain = "LimitesAllEtaBARREL.txt";
  if(variableX == "Photon_SC_Eta" && EndCaps == 1) EndcapsChain = "LimitesAllEtaENDCAPS.txt";
  if(variableX == "Photon_SC_brem") EndcapsChain = "LimitesAllBrem.txt";

  n = NbLignesFichier(EndcapsChain.c_str()) - 1;
  ifstream monFlux(EndcapsChain.c_str());
  //**************** initialization ****************************



  //******************* add tree ******************************
  TChain *DataChain = new TChain("miniTree");
  DataChain->Add("/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree/OutputMiniTree/miniTree_53X_Run2012AB.root");
  //******************* add tree *******************************
 


  //******************* directory set **************************
  string nomDossier = Form("/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree/FitSvalue/Results/%s_%s_%s_%s_%dPourcents/",variableX.c_str(),DossierCracks.c_str(),isMCChain.c_str(),SurfaceMethod.c_str(),FitPourcentage);

  string nomDossierFits = "";
  string nomDossierTGraph = "";
  string nomFichier = "";
  //******************* directory set **************************


  //******************* cut set ********************************
  TString temp = "";
  TString tempVarChain = "";


  for(int j = 0; j < n; j++){    //mettre j<n

     TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600); 
     
     TH1D *MC = new TH1D("MC","MC", 1000, 0.8, 1.2);   

      
     temp.Clear();
     tempVarChain.Clear();


     temp += "isLooseMMG == 1";                  
     if(EndCaps == 0 && r9sup == 1 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 > 0.94";
     if(EndCaps == 0 && r9sup == 0 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 < 0.94";
     if(EndCaps == 1 && r9sup == 1 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 > 0.95";
     if(EndCaps == 1 && r9sup == 0 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 < 0.95";

     if(EndCaps == 0 && r9sup == 2 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442";
     if(EndCaps == 1 && r9sup == 2 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566";

     if(EndCaps == 0 && r9sup == 1 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 > 0.94";
     if(EndCaps == 0 && r9sup == 0 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 < 0.94";
     if(EndCaps == 1 && r9sup == 1 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 > 0.95";
     if(EndCaps == 1 && r9sup == 0 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 < 0.95";

     if(EndCaps == 0 && etaCracks == false) temp += " && ( (abs(Photon_SC_Eta)) > 0.018 || ( (abs(Photon_SC_Eta)) < 0.423 && (abs(Photon_SC_Eta)) > 0.461 ) || ( (abs(Photon_SC_Eta)) < 0.770 && (abs(Photon_SC_Eta)) > 0.806 ) || ( (abs(Photon_SC_Eta)) < 1.127 && (abs(Photon_SC_Eta)) > 1.163 ) )";         

     if(EndCaps == 0 && phiCracks == false){
                        temp += Form(" && ((abs(Photon_SC_Phi) + (%f)) < (%f) ||",phiOffset,phiCrackSize);
                        temp += Form(" ( (1.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (1.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (2.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (2.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (3.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (3.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (4.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (4.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (5.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (5.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (6.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (6.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (7.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (7.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (8.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (8.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (9.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) ) )",phiCrackPosition,phiCrackSize,phiOffset);

     }


     temp +=" && " + variableX + " > ";
     double nombre;
     if(j == 0)
     {
             monFlux >> nombre;
             MinVar[j] = nombre;
             monFlux >> nombre;
             MaxVar[j] = nombre;
     }

     if(j > 0)
     {
             MinVar[j] = MaxVar[j-1];
             monFlux >> nombre;
             MaxVar[j] = nombre;

     }

     temp += MinVar[j];
     temp += " && " + variableX + " <= ";
     //temp += " && Photon_Et < ";
     temp += MaxVar[j];
     //**************************** cut set ***************************



     //********************* mmg_s and variable draw ******************
     DataChain->Draw("mmg_s>>MC",temp);
     cout<<"MC->GetMaximum= "<<MC->GetMaximumBin() * 0.001 + 0.5<<endl;
     maxDistri = MC->GetMaximumBin() * 0.001 - 0.5;
     lastX = 0.5;
     Entries = MC->GetEntries();
     meanMC = MC->GetMean();
     cout<<"MC->GetMean() = "<<MC->GetMean()<<endl;


     TH1D *VarXvalue = new TH1D();
     if(variableX == "Photon_SC_rawEt") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 250.0);
     if(variableX == "Photon_Et") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 250.0);
     if(variableX == "Photon_E") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 1000.0);
     if(variableX == "Photon_SC_Eta") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, -3.0, 3.0);
     if(variableX == "Photon_SC_brem") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 15.0);
     tempVarChain = variableX + ">>VarXvalue";
     DataChain->Draw(tempVarChain,temp);
     xValue[j] = VarXvalue->GetMean(1);
     cout<<endl<<"xValue = "<<xValue[j]<<endl;
     c1->Clear();
     //****************** mmg_s and variable draw ************************    



 
     //****************** set fit range *************************** 
     pourcentage = FitPourcentageD;
     double MinRange = 0;
     double MaxRange = 0;
     RangeEstimator3(pourcentage, DataChain, temp, EndCaps, &MinRange, &MaxRange);
     // ***************** set fit range ***************************




     // ***************** start to fit!!! *************************     
     if(fitRoo == 1){
          mean = MC->GetMean();
          rms = MC->GetRMS();
          TF1 * f = new TF1();
          RooPlot* Erawframe = new RooPlot(-1.0,1.0);
          nomFichier = "Svalue";
          nomDossierFits = nomDossier; 


          if(nomFitMethode == "RooVoigtian2"){
              if(r9sup == 2) nomDossierFits += "Fits/Voigtian_r9All";
              if(r9sup == 1) nomDossierFits += "Fits/Voigtian_r9sup";
              if(r9sup == 0) nomDossierFits += "Fits/Voigtian_r9inf";
              RooVoigtian2(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, &minLogLikelihood, &differenceLogLikelihood, MC, j, EndCaps, temp, DataChain, mean, rms, c1, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
          }


          MeanTab[j] = mean_value; //Mean
          MeanErrorTab[j] = mean_error; //MeanError 
          ChiSquareTab[j] = ChiSquare; //Chi2
          SigmaEffTab[j] = sigmaEff_value;
          SigmaTab[j] = sigma_value;
          SigmaErrorTab[j] = sigma_value_error;
          SigmaRTab[j] = sigmaR_value;
          SigmaLTab[j] = sigmaL_value;
          JanChiSquareTab[j] = JanChi2;
          PValueTab[j] = pValue;
          DegreesOfFreedomTab[j] = DegreesOfFreedom;
          minLogLikelihoodTab[j] = minLogLikelihood;
          differenceLogLikelihoodTab[j] = differenceLogLikelihood;  

          f->Delete();
          Erawframe->Delete();
      }
      //******************* start to fit!!! *****************************



      //************** erase all *********************************
      MC->Delete();
      VarXvalue->Delete(); 
      delete c1;

  }  //end fit loop
  
  DataChain->Delete();
  



  TCanvas* c2 = new TCanvas("c2", "c2",0,0,600,600);


//Generation des bins

  for(int k = 0; k < n ; k++)
  {
          xErrorR[k] = MaxVar[k] - xValue[k];
          xErrorL[k] = xValue[k] - MinVar[k];

          ////////// Modifs mean et sigma plots de Jan a enlever si non comparaison  //////////
          MeanTab[k] = MeanTab[k] * 100.0;
          MeanErrorTab[k] = MeanErrorTab[k] * 100.0;
          SigmaTab[k] = SigmaTab[k] * 100.0;
          SigmaErrorTab[k] = SigmaErrorTab[k] * 100.0;

  }






  //********************** txt output ****************************     
  if(EndCaps == 0) nomDossierFits += "_EB/";
  if(EndCaps == 1) nomDossierFits += "_EE/";

  for(int m = 0; m <n ; m++) 
        {    

                ofstream monFlux2(Form("%sMeanTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux2 << MeanTab[m] <<endl;
                monFlux2.close();

                ofstream monFlux3(Form("%sMeanErrorTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux3 << MeanErrorTab[m] <<endl;
                monFlux3.close();

                ofstream monFlux4(Form("%sSigmaRTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux4 << SigmaRTab[m] <<endl;
                monFlux4.close();

                ofstream monFlux5(Form("%sSigmaLTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux5 << SigmaLTab[m] <<endl;
                monFlux5.close();

                ofstream monFlux6(Form("%sSigmaTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux6 << SigmaTab[m] <<endl;
                monFlux6.close();

                ofstream monFlux7(Form("%sMinVar.txt",nomDossierFits.c_str()), ios::app);
                monFlux7 << MinVar[m] <<endl;
                monFlux7.close();

                ofstream monFlux8(Form("%sMaxVar.txt",nomDossierFits.c_str()), ios::app);
                monFlux8 << MaxVar[m] <<endl;
                monFlux8.close();

                ofstream monFlux9(Form("%sxValue.txt",nomDossierFits.c_str()), ios::app);
                monFlux9 << xValue[m] <<endl;
                monFlux9.close();

                ofstream monFlux10(Form("%sxErrorR.txt",nomDossierFits.c_str()), ios::app);
                monFlux10 << xErrorR[m] <<endl;
                monFlux10.close();

                ofstream monFlux11(Form("%sxErrorL.txt",nomDossierFits.c_str()), ios::app);
                monFlux11 << xErrorL[m] <<endl;
                monFlux11.close();

                ofstream monFlux12(Form("%syValue.txt",nomDossierFits.c_str()), ios::app);
                monFlux12 << yValue[m] <<endl;
                monFlux12.close();

                ofstream monFlux13(Form("%sxErrorRight.txt",nomDossierFits.c_str()), ios::app);
                monFlux13 << xErrorRight[m] <<endl;
                monFlux13.close();
     
                ofstream monFlux14(Form("%sxErrorLeft.txt",nomDossierFits.c_str()), ios::app);
                monFlux14 << xErrorLeft[m] <<endl;
                monFlux14.close();

                ofstream monFlux15(Form("%sChiSquareTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux15 << ChiSquareTab[m] <<endl;
                monFlux15.close();

                ofstream monFlux16(Form("%sSigmaEffTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux16 << SigmaEffTab[m] <<endl;
                monFlux16.close();

                ofstream monFlux17(Form("%sSigmaErrorTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux17 << SigmaErrorTab[m] <<endl;
                monFlux17.close();

                ofstream monFlux18(Form("%sJanChiSquareTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux18 << JanChiSquareTab[m] <<endl;
                monFlux18.close();

                ofstream monFlux19(Form("%sPValueTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux19 << PValueTab[m] <<endl;
                monFlux19.close();

                ofstream monFlux20(Form("%sDegreesOfFreedomTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux20 << DegreesOfFreedomTab[m] <<endl;
                monFlux20.close();

                ofstream monFlux21(Form("%sminLogLikelihoodTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux21 << minLogLikelihoodTab[m] <<endl;
                monFlux21.close();
     
                ofstream monFlux22(Form("%sdifferenceLogLikelihoodTab.txt",nomDossierFits.c_str()), ios::app);
                monFlux22 << differenceLogLikelihoodTab[m] <<endl;
                monFlux22.close();


        }
     //******************** txt output *************************





//Graphe E_{RAW*Ceta*PtCor}/E_{TRUE} Converted (MC) vs Var MC 2gammas
//******************************** graph directory and preparation***********************
        nomDossierTGraph = nomDossier;
        if(r9sup == 2)
        {    
                if(nomFitMethode == "RooCrystalBall") nomDossierTGraph += "TGraph/CB_r9All";
                if(nomFitMethode == "RooLogNormal") nomDossierTGraph += "TGraph/LN_r9All";
                if(nomFitMethode == "RooGauss") nomDossierTGraph += "TGraph/G_r9All";         
                if(nomFitMethode == "RooGamma") nomDossierTGraph += "TGraph/Gamma_r9All";             
                if(nomFitMethode == "RooBifurcatedGauss") nomDossierTGraph += "TGraph/BFG_r9All";
                if(nomFitMethode == "RooSumGauss") nomDossierTGraph += "TGraph/SumG_r9All";
                if(nomFitMethode == "RooGenericPDF") nomDossierTGraph += "TGraph/CBG_r9All";
                if(nomFitMethode == "RooLandau2") nomDossierTGraph += "TGraph/RooLandau_r9All";
                if(nomFitMethode == "RooLandauConvGaussian") nomDossierTGraph += "TGraph/LandauConvGaus_r9All";
                if(nomFitMethode == "RooKernel") nomDossierTGraph += "TGraph/Kernel_r9All";
                if(nomFitMethode == "RooVoigtian2") nomDossierTGraph += "TGraph/Voigtian_r9All";
        }    

        if(r9sup == 1)
        {    
                if(nomFitMethode == "RooCrystalBall") nomDossierTGraph += "TGraph/CB_r9sup";
                if(nomFitMethode == "RooLogNormal") nomDossierTGraph += "TGraph/LN_r9sup";
                if(nomFitMethode == "RooGauss") nomDossierTGraph += "TGraph/G_r9sup";                if(nomFitMethode == "RooGamma") nomDossierTGraph += "TGraph/Gamma_r9sup";
                if(nomFitMethode == "RooBifurcatedGauss") nomDossierTGraph += "TGraph/BFG_r9sup";
                if(nomFitMethode == "RooSumGauss") nomDossierTGraph += "TGraph/SumG_r9sup";
                if(nomFitMethode == "RooGenericPDF") nomDossierTGraph += "TGraph/CBG_r9sup";
                if(nomFitMethode == "RooLandau2") nomDossierTGraph += "TGraph/RooLandau_r9sup";
                if(nomFitMethode == "RooLandauConvGaussian") nomDossierTGraph += "TGraph/LandauConvGaus_r9sup";
                if(nomFitMethode == "RooKernel") nomDossierTGraph += "TGraph/Kernel_r9sup";
                if(nomFitMethode == "RooVoigtian2") nomDossierTGraph += "TGraph/Voigtian_r9sup";
        }            

        if(r9sup == 0)
        {
                if(nomFitMethode == "RooCrystalBall") nomDossierTGraph += "TGraph/CB_r9inf";
                if(nomFitMethode == "RooLogNormal") nomDossierTGraph += "TGraph/LN_r9inf";
                if(nomFitMethode == "RooGauss") nomDossierTGraph += "TGraph/G_r9inf";                if(nomFitMethode == "RooGamma") nomDossierTGraph += "TGraph/Gamma_r9inf";
                if(nomFitMethode == "RooBifurcatedGauss") nomDossierTGraph += "TGraph/BFG_r9inf";
                if(nomFitMethode == "RooSumGauss") nomDossierTGraph += "TGraph/SumG_r9inf";
                if(nomFitMethode == "RooGenericPDF") nomDossierTGraph += "TGraph/CBG_r9inf";  
                if(nomFitMethode == "RooLandau2") nomDossierTGraph += "TGraph/RooLandau_r9inf";
                if(nomFitMethode == "RooLandauConvGaussian") nomDossierTGraph += "TGraph/LandauConvGaus_r9inf";
                if(nomFitMethode == "RooKernel") nomDossierTGraph += "TGraph/Kernel_r9inf";
                if(nomFitMethode == "RooVoigtian2") nomDossierTGraph += "TGraph/Voigtian_r9inf";
        }            

  //************************** graph directory *********************************************




  // ************************** draw graph X is variable while Y is fitted mean value *******************************
    TGraphAsymmErrors * MCtg = new TGraphAsymmErrors(n,xValue, MeanTab, xErrorL, xErrorR, MeanErrorTab, MeanErrorTab);
    c2->ToggleEventStatus();
    MCtg->SetTitle("");
    MCtg->SetLineColor(4);
    if(variableX == "Photon_SC_rawEt") MCtg->GetXaxis()->SetTitle("E_{T RAW}");
    if(variableX == "Photon_Et") MCtg->GetXaxis()->SetTitle("P_{T}");
    if(variableX == "Photon_E") MCtg->GetXaxis()->SetTitle("E");
    if(variableX == "Photon_SC_Eta") MCtg->GetXaxis()->SetTitle("#eta");
    if(variableX == "Photon_SC_brem") MCtg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
    //MCtg->GetXaxis()->SetTitle("Var");
    MCtg->GetXaxis()->SetLabelFont(42);
    MCtg->GetXaxis()->SetTitleFont(42);
    MCtg->GetXaxis()->SetLabelSize(0.03);
    MCtg->GetYaxis()->SetTitle("s_{RECO} (%)");
    MCtg->GetYaxis()->SetLabelFont(42);
    MCtg->GetYaxis()->SetTitleOffset(1.24);
    MCtg->GetYaxis()->SetTitleFont(42);
    MCtg->GetYaxis()->SetLabelSize(0.03);
    //MCtg->SetMarkerColor(4);
    //MCtg->SetMarkerStyle(21);
    //MCtg->SetMarkerSize(0.6);
    MCtg->Draw("AP");


    TPaveText *pt = new TPaveText(0.1073826,0.9125874,0.897651,0.9912587,"blNDC");
    pt->SetName("title");
    pt->SetBorderSize(2);
    pt->SetFillColor(kWhite);
    pt->SetTextFont(42);
    TText * text2;
    text2 = pt->AddText("");
    pt->Draw();

    
    TLatex *textF = new TLatex();
    textF->SetNDC();
    textF->SetTextAlign(11);
    //textF->SetTextFont(42);
    textF->SetTextSizePixels(17);
    textF->SetTextSize(0.038);
    //textF->DrawLatex(0.135, 0.93, "s_{RECO} vs Var, MC Z#rightarrow#mu#mu#gamma");
    gStyle->SetPadBorderMode(0);


    TLatex *textL = new TLatex();
    textL = new TLatex();
    textL->SetNDC();
    textL->SetTextAlign(11);
    textL->SetTextFont(42);
    textL->SetTextSizePixels(17);
    textL->SetTextSize(0.028);
    //if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
    //if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
    //if(EndCaps == 2) textL->DrawLatex(0.75, 0.83, "");
    textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
    if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
    if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
    if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
    if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
    if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
    if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");
    gStyle->SetPadBorderMode(0);


    if(EndCaps == 0) MCtg->GetYaxis()->SetRangeUser(ymin1overKrecoClassical,ymax1overKrecoClassical);//a changer
    if(EndCaps == 1) MCtg->GetYaxis()->SetRangeUser(ymin1overKrecoClassical,ymax1overKrecoClassical);//a changer
    if(EndCaps == 2) MCtg->GetYaxis()->SetRangeUser(ymin1overKrecoClassical,ymax1overKrecoClassical);
    MCtg->GetXaxis()->SetLimits(xmin1overKrecoClassical,xmax1overKrecoClassical);
    c2->SetTickx(1);
    c2->SetTicky(1);
    c2->SetGridx(1);
    c2->SetGridy(1);
    c2->Modified();
    c2->cd();
    c2->SetSelected(c2);
    c2->ToggleToolBar();

    nomFichier = "VarVsMeanValue";
    enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c2);
    c2->Clear();
 



     //********************* PValue VS PT TGraph *************************
     TGraphAsymmErrors * PValueTg = new TGraphAsymmErrors(n,xValue, PValueTab, xErrorL, xErrorR, PValueErrorTab, PValueErrorTab);

     c2->ToggleEventStatus();
     PValueTg->SetFillColor(1);
     PValueTg->SetLineColor(4);
     PValueTg->SetMarkerColor(4);
     PValueTg->SetMarkerStyle(21);
     PValueTg->SetMarkerSize(0.6);
     PValueTg->SetLineColor(4);
     PValueTg->SetTitle("");
     PValueTg->SetLineColor(4);
     if(variableX == "Photon_SC_rawEt") PValueTg->GetXaxis()->SetTitle("E_{T RAW}");
     if(variableX == "Photon_Et") PValueTg->GetXaxis()->SetTitle("P_{T}");
     if(variableX == "Photon_E") PValueTg->GetXaxis()->SetTitle("E");
     if(variableX == "Photon_SC_Eta") PValueTg->GetXaxis()->SetTitle("#eta");
     if(variableX == "Photon_SC_brem") PValueTg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
     //PValueTg->GetXaxis()->SetTitle("Var");
     PValueTg->GetXaxis()->SetLabelFont(42);
     PValueTg->GetXaxis()->SetTitleFont(42);
     PValueTg->GetXaxis()->SetLabelSize(0.03);
     PValueTg->GetYaxis()->SetTitle("p-value");
     PValueTg->GetYaxis()->SetLabelFont(42);
     PValueTg->GetYaxis()->SetTitleOffset(1.24);
     PValueTg->GetYaxis()->SetTitleFont(42);
     PValueTg->GetYaxis()->SetLabelSize(0.03);
     //PValueTg->SetMarkerColor(4);
     //PValueTg->SetMarkerStyle(21);
     //PValueTg->SetMarkerSize(0.6);
     PValueTg->Draw("AP");


     pt = new TPaveText(0.1073826,0.9125874,0.897651,0.9912587,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     pt->SetTextFont(42);
     //TText * text2;
     text2 = pt->AddText("");
     pt->Draw();

     //TLatex *textF = new TLatex();
     textF->SetNDC();
     textF->SetTextAlign(11);
     //textF->SetTextFont(42);
     textF->SetTextSizePixels(17);
     textF->SetTextSize(0.038);

     
     gStyle->SetPadBorderMode(0);


     textL = new TLatex();
     textL = new TLatex();
     textL->SetNDC();
     textL->SetTextAlign(11);
     textL->SetTextFont(42);
     textL->SetTextSizePixels(17);
     textL->SetTextSize(0.028);
     //if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
     //if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
     //if(EndCaps == 2) textL->DrawLatex(0.75, 0.83, "");
     textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
     if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
     if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
     if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
     if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
     if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
     if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");
     gStyle->SetPadBorderMode(0);

     PValueTg->GetYaxis()->SetRangeUser(yminPValue,ymaxPValue);//a changer
     
     PValueTg->GetXaxis()->SetLimits(xminPValue,xmaxPValue);

     c2->SetTickx(1);
     c2->SetTicky(1);
     c2->SetGridx(1);
     c2->SetGridy(1);
     c2->Modified();
     c2->cd();
     c2->SetSelected(c2);
     c2->ToggleToolBar();

     nomFichier = "PvalueVsVar";
     enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c2); 

     c2->Clear();





     //*****************Graph chi2 et sigmaRL*****************//

     gROOT->SetStyle("Plain");
     gStyle->SetOptTitle(0);


     //c2 = new TCanvas("c2", "c2",0,0,600,600);
     c2->Range(0,0,1,1);
     c2->SetBorderSize(2);
     c2->SetFrameFillColor(0);

     TPad *c2_1 = new TPad("c2_1", "c2_1",0.01,0.51,0.99,0.99);
     c2_1->Draw();
     c2_1->cd();
     c2_1->Range(-12.5,-0.4455269,12.5,1.229671);
     c2_1->SetBorderSize(2);
     c2_1->SetFrameFillColor(0);
     //c2_1->SetLogy();
     c2_1->SetTickx(1);
     c2_1->SetTicky(1);
     c2_1->SetGridx(1);
     c2_1->SetGridy(1);

     TGraph* chi2 = new TGraph(n,xValue,ChiSquareTab);
     chi2->SetFillColor(1);
     chi2->SetLineColor(4);
     chi2->SetMarkerColor(4);
     chi2->SetMarkerStyle(21);
     chi2->SetMarkerSize(0.6);

     chi2->GetYaxis()->SetTitle("#chi^{2} / ndf");
     chi2->GetYaxis()->CenterTitle(true);
     chi2->GetYaxis()->SetLabelSize(0.06);
     chi2->GetXaxis()->SetLabelSize(0.06);
     chi2->GetYaxis()->SetTitleSize(0.07);
     chi2->GetYaxis()->SetTitleOffset(0.5);
     chi2->GetXaxis()->SetLabelFont(42);
     chi2->GetYaxis()->SetTitleFont(42);
     chi2->GetYaxis()->SetLabelFont(42);
     chi2->Draw("AP");
     chi2->SetTitle("");
     chi2->GetYaxis()->SetRangeUser(yminChi2,ymaxChi2);
     chi2->GetXaxis()->SetLimits(xminChi2,xmaxChi2);
     //chi2->GetXaxis()->SetRangeUser(0,17);

     pt = new TPaveText(0,0,0,0,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     TText * text4 = pt->AddText("");
     pt->Draw();


     TLine *line = new TLine(0,1,xmaxChi2,1);
     line->SetLineStyle(3);
     line->Draw();

     c2_1->Modified();
     c2->cd();


     TPad *c2_2 = new TPad("c2_2", "c2_2",0.01,0.01,0.99,0.49);
     c2_2->Draw();
     c2_2->cd();
     c2_2->Range(-12.5,-1.375519,12.5,1.380519);
     c2_2->SetBorderSize(2);
     c2_2->SetFrameFillColor(0);
     c2_2->SetTickx(1);
     c2_2->SetTicky(1);
     c2_2->SetGridx(1);
     c2_2->SetGridy(1);






     //********************** not used **************************
/*
     TGraphAsymmErrors * sigmaRLGraph = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaLTab, SigmaRTab);
     sigmaRLGraph->SetFillColor(4);
     sigmaRLGraph->SetLineColor(4);
     sigmaRLGraph->GetYaxis()->SetTitle("_{L} - #sigma_{R}");
     sigmaRLGraph->GetYaxis()->CenterTitle(true);
     sigmaRLGraph->GetYaxis()->SetLabelSize(0.06);
     sigmaRLGraph->GetXaxis()->SetLabelSize(0.06);
     sigmaRLGraph->GetYaxis()->SetTitleSize(0.07);
     sigmaRLGraph->GetYaxis()->SetTitleOffset(0.5);
     sigmaRLGraph->GetXaxis()->SetLabelFont(42);
     sigmaRLGraph->GetYaxis()->SetTitleFont(42);
     sigmaRLGraph->GetYaxis()->SetLabelFont(42);
     sigmaRLGraph->SetLineWidth(1);
     sigmaRLGraph->Draw("AP");

     sigmaRLGraph->SetTitle("");

     

     pt = new TPaveText(0,0,0,0,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     TText *text5 = pt->AddText("");
     pt->Draw();

     sigmaRLGraph->GetYaxis()->SetRangeUser(yminSigmaLSigmaR,ymaxSigmaLSigmaR);
     sigmaRLGraph->GetXaxis()->SetLimits(xminSigmaLSigmaR,xmaxSigmaLSigmaR);
     TLine *line2 = new TLine(0,0,xmaxSigmaLSigmaR,0);
     line2->SetLineStyle(3);
     line2->Draw();


     c2_2->Modified();
     c2->cd();
     c2->Modified();
     c2->cd();
     c2->SetSelected(c2);

     nomFichier = "Chi2SigmaLRVsVar";
     enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c2);
     

     c2->Clear();
     //c2->Delete();
     chi2->Delete();
     sigmaRLGraph->Delete();
*/




     //*****************Graph sigma et sigmaEff*****************//
/*
     gStyle->SetOptTitle(0);


     //c2 = new TCanvas("c2", "c2",0,0,600,600);
     c2->Range(0,0,1,1);
     c2->SetBorderSize(2);
     c2->SetFrameFillColor(0);

     c2_1 = new TPad("c2_1", "c2_1",0.01,0.51,0.99,0.99);
     c2_1->Draw();
     c2_1->cd();
     c2_1->Range(-12.5,-0.4455269,12.5,1.229671);
     c2_1->SetBorderSize(2);
     c2_1->SetFrameFillColor(0);
     //c2_1->SetLogy();
     c2_1->SetTickx(1);
     c2_1->SetTicky(1);
     c2_1->SetGridx(1);
     c2_1->SetGridy(1);

     TGraphAsymmErrors * sigmaGraph = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaTab, SigmaTab);
     sigmaGraph->SetFillColor(4);
     sigmaGraph->SetLineColor(4);
     sigmaGraph->GetYaxis()->SetTitle(" (%)");
     sigmaGraph->GetYaxis()->CenterTitle(true);
     sigmaGraph->GetYaxis()->SetLabelSize(0.06);
     sigmaGraph->GetXaxis()->SetLabelSize(0.06);
     sigmaGraph->GetYaxis()->SetTitleSize(0.07);
     sigmaGraph->GetYaxis()->SetTitleOffset(0.5);

     sigmaGraph->GetXaxis()->SetLabelFont(42);
     sigmaGraph->GetYaxis()->SetTitleFont(42);
     sigmaGraph->GetYaxis()->SetLabelFont(42);
     sigmaGraph->SetLineWidth(1);
     sigmaGraph->Draw("AP");

     sigmaGraph->SetTitle("");

     pt = new TPaveText(0,0,0,0,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     TText *text6 = pt->AddText("");
     pt->Draw();

     sigmaGraph->GetYaxis()->SetRangeUser(yminSigma,ymaxSigma);


     sigmaGraph->GetXaxis()->SetLimits(xminSigma,xmaxSigma);
     TLine *line3 = new TLine(0,0,xmaxSigma,0);
     line3->SetLineStyle(3);
     line3->Draw();


     c2_1->Modified();
     c2->cd();


     c2_2 = new TPad("c2_2", "c2_2",0.01,0.01,0.99,0.49);
     c2_2->Draw();
     c2_2->cd();
     c2_2->Range(-12.5,-1.375519,12.5,1.380519);
     c2_2->SetBorderSize(2);
     c2_2->SetFrameFillColor(0);
     c2_2->SetTickx(1);
     c2_2->SetTicky(1);
     c2_2->SetGridx(1);
     c2_2->SetGridy(1);
*/




     //******************** not used *********************************
/*
     TGraphAsymmErrors * sigmaEffGraph = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaEffTab, SigmaEffTab);
     sigmaEffGraph->SetFillColor(4);
     sigmaEffGraph->SetLineColor(4);
     sigmaEffGraph->GetYaxis()->SetTitle("_{eff}");
     sigmaEffGraph->GetYaxis()->CenterTitle(true);
     sigmaEffGraph->GetYaxis()->SetLabelSize(0.06);
     sigmaEffGraph->GetXaxis()->SetLabelSize(0.06);
     sigmaEffGraph->GetYaxis()->SetTitleSize(0.07);
     sigmaEffGraph->GetYaxis()->SetTitleOffset(0.5);
     sigmaEffGraph->GetXaxis()->SetLabelFont(42);
     sigmaEffGraph->GetYaxis()->SetTitleFont(42);
     sigmaEffGraph->GetYaxis()->SetLabelFont(42);
     sigmaEffGraph->SetLineWidth(1);
     sigmaEffGraph->Draw("AP");

     sigmaEffGraph->SetTitle("");

     pt = new TPaveText(0,0,0,0,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     TText *text7 = pt->AddText("");
     pt->Draw();

     sigmaEffGraph->GetYaxis()->SetRangeUser(yminSigmaEff,ymaxSigmaEff);


     sigmaEffGraph->GetXaxis()->SetLimits(xminSigmaEff,xmaxSigmaEff);
     TLine *line4 = new TLine(0,0,xmaxSigmaEff,0);
     line4->SetLineStyle(3);
     line4->Draw();


     c2_2->Modified();
     c2->cd();
     c2->Modified();
     c2->cd();
     c2->SetSelected(c2);

     nomFichier = "SigmaSigmaEffVsVar";
     enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c2);
     c2->Clear();
     //c2->Delete();
     sigmaGraph->Delete();
     sigmaEffGraph->Delete();
*/




     //********************* TH1D sigmaEff ****************************
     TGraph* sigmaEffTg = new TGraph(n,xValue,SigmaEffTab);

     pt = new TPaveText(0.1073826,0.9125874,0.897651,0.9912587,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     pt->SetTextFont(42);
     TText * text = pt->AddText("");
     pt->Draw();



     sigmaEffTg->SetFillColor(1);
     sigmaEffTg->SetLineColor(4);
     sigmaEffTg->SetMarkerColor(4);
     sigmaEffTg->SetMarkerStyle(21);
     sigmaEffTg->SetMarkerSize(0.6);

     sigmaEffTg->GetYaxis()->SetTitle("_{eff}");
     sigmaEffTg->GetYaxis()->SetLabelSize(0.03);
     sigmaEffTg->GetXaxis()->SetLabelSize(0.03);
     //sigmaEffTg->GetYaxis()->SetTitleSize(0.07);
     sigmaEffTg->GetYaxis()->SetTitleOffset(1.24);
     if(variableX == "Photon_SC_rawEt") sigmaEffTg->GetXaxis()->SetTitle("E_{T RAW}");
     if(variableX == "Photon_Et") sigmaEffTg->GetXaxis()->SetTitle("P_{T}");
     if(variableX == "Photon_E") sigmaEffTg->GetXaxis()->SetTitle("E");
     if(variableX == "Photon_SC_Eta") sigmaEffTg->GetXaxis()->SetTitle("#eta");
     if(variableX == "Photon_SC_brem") sigmaEffTg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
     //sigmaEffTg->GetXaxis()->SetTitle("Var");
     sigmaEffTg->GetXaxis()->SetLabelFont(42);
     sigmaEffTg->GetYaxis()->SetTitleFont(42);
     sigmaEffTg->GetYaxis()->SetLabelFont(42);
     sigmaEffTg->Draw("AP");
     sigmaEffTg->SetTitle("");

     sigmaEffTg->GetYaxis()->SetRangeUser(yminSigmaEffTg,ymaxSigmaEffTg);
     sigmaEffTg->GetXaxis()->SetLimits(xminSigmaEffTg,xmaxSigmaEffTg);
     //sigmaEffTg->GetXaxis()->SetRangeUser(0,17);

     pt = new TPaveText(0,0,0,0,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     TText *text8 = pt->AddText("");
     pt->Draw();
     //if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
     //if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
     textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
     if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
     if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
     if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
     if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
     if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
     if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");

     c2->SetTickx(1);
     c2->SetTicky(1);
     c2->SetGridx(1);
     c2->SetGridy(1);
     c2->Modified();
     c2->cd();
     c2->SetSelected(c2);
     c2->ToggleToolBar();

     nomFichier = "SigmaEffVsVarTg";
     enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c2);

     c2->Clear();


      



     //***************************** TGraphAsym sigma CB ****************
     TGraphAsymmErrors * SigmaCBtg = new TGraphAsymmErrors(n,xValue, SigmaTab, xErrorL, xErrorR, SigmaErrorTab, SigmaErrorTab);

     SigmaCBtg->SetFillColor(1);
     SigmaCBtg->SetLineColor(4);
     //SigmaCBtg->SetMarkerColor(4);
     //SigmaCBtg->SetMarkerStyle(21);
     //SigmaCBtg->SetMarkerSize(0.6);

     SigmaCBtg->GetYaxis()->SetTitle("_{CB} (%)");
     SigmaCBtg->GetYaxis()->SetLabelSize(0.03);
     SigmaCBtg->GetXaxis()->SetLabelSize(0.03);
     //SigmaCBtg->GetYaxis()->SetTitleSize(0.07);
     SigmaCBtg->GetYaxis()->SetTitleOffset(1.24);
     if(variableX == "Photon_SC_rawEt") SigmaCBtg->GetXaxis()->SetTitle("E_{T RAW}");
     if(variableX == "Photon_Et") SigmaCBtg->GetXaxis()->SetTitle("P_{T}");
     if(variableX == "Photon_E") SigmaCBtg->GetXaxis()->SetTitle("E");
     if(variableX == "Photon_SC_Eta") SigmaCBtg->GetXaxis()->SetTitle("#eta");
     if(variableX == "Photon_SC_brem") SigmaCBtg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
     //SigmaCBtg->GetXaxis()->SetTitle("Var");
     SigmaCBtg->GetXaxis()->SetLabelFont(42);
     SigmaCBtg->GetYaxis()->SetTitleFont(42);
     SigmaCBtg->GetYaxis()->SetLabelFont(42);
     SigmaCBtg->Draw("AP");


     SigmaCBtg->GetYaxis()->SetRangeUser(yminSigmaTg,ymaxSigmaTg);
     SigmaCBtg->GetXaxis()->SetLimits(xminSigmaTg,xmaxSigmaTg);
     //SigmaCBtg->GetXaxis()->SetRangeUser(0,17);
     //if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
     //if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
     textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
     if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
     if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
     if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
     if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
     if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
     if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");

     c2->SetTickx(1);
     c2->SetTicky(1);
     c2->SetGridx(1);
     c2->SetGridy(1);
     c2->Modified();
     c2->cd();
     c2->SetSelected(c2);
     c2->ToggleToolBar();

     nomFichier = "SigmaCBVsVarTg";
     enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c2);

     c2->Clear();

     



     //***************************** Graph Janchi2 et p-value ****************
     c2->Range(0,0,1,1);
     c2->SetBorderSize(2);
     c2->SetFrameFillColor(0);

     c2_1 = new TPad("c2_1", "c2_1",0.01,0.51,0.99,0.99);
     c2_1->Draw();
     c2_1->cd();
     c2_1->Range(-12.5,-0.4455269,12.5,1.229671);
     c2_1->SetBorderSize(2);
     c2_1->SetFrameFillColor(0);
     //c2_1->SetLogy();
     c2_1->SetTickx(1);
     c2_1->SetTicky(1);
     c2_1->SetGridx(1);
     c2_1->SetGridy(1);


     TGraph* Janchi2 = new TGraph(n,xValue,JanChiSquareTab);
     Janchi2->SetFillColor(1);
     Janchi2->SetLineColor(4);
     Janchi2->SetMarkerColor(4);
     Janchi2->SetMarkerStyle(21);
     Janchi2->SetMarkerSize(0.6);

     Janchi2->GetYaxis()->SetTitle("Jan #chi^{2} / ndf");
     Janchi2->GetYaxis()->CenterTitle(true);
     Janchi2->GetYaxis()->SetLabelSize(0.06);
     Janchi2->GetXaxis()->SetLabelSize(0.06);
     Janchi2->GetYaxis()->SetTitleSize(0.07);
     Janchi2->GetYaxis()->SetTitleOffset(0.5);
     Janchi2->GetXaxis()->SetLabelFont(42);
     Janchi2->GetYaxis()->SetTitleFont(42);
     Janchi2->GetYaxis()->SetLabelFont(42);
     Janchi2->Draw("AP");
     Janchi2->SetTitle("");

     Janchi2->GetYaxis()->SetRangeUser(yminJanChi2,ymaxJanChi2);
     Janchi2->GetXaxis()->SetLimits(xminJanChi2,xmaxJanChi2);
     //Janchi2->GetXaxis()->SetRangeUser(0,17);


     pt = new TPaveText(0,0,0,0,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     //TText * text11 = pt->AddText("");
     pt->Draw();


     line = new TLine(0,1,xmaxJanChi2,1);
     line->SetLineStyle(3);
     line->Draw();

     c2_1->Modified();
     c2->cd();


     c2_2 = new TPad("c2_2", "c2_2",0.01,0.01,0.99,0.49);
     c2_2->Draw();
     c2_2->cd();
     c2_2->Range(-12.5,-1.375519,12.5,1.380519);
     c2_2->SetBorderSize(2);
     c2_2->SetFrameFillColor(0);
     c2_2->SetTickx(1);
     c2_2->SetTicky(1);
     c2_2->SetGridx(1);
     c2_2->SetGridy(1);


     TGraph* PValueGraph = new TGraph(n,xValue,PValueTab);
     PValueGraph->SetFillColor(1);
     PValueGraph->SetLineColor(4);
     PValueGraph->SetMarkerColor(4);
     PValueGraph->SetMarkerStyle(21);
     PValueGraph->SetMarkerSize(0.6);
     PValueGraph->SetLineColor(4);
     PValueGraph->GetYaxis()->SetTitle("p-value");
     PValueGraph->GetYaxis()->CenterTitle(true);
     PValueGraph->GetYaxis()->SetLabelSize(0.06);
     PValueGraph->GetXaxis()->SetLabelSize(0.06);
     PValueGraph->GetYaxis()->SetTitleSize(0.07);
     PValueGraph->GetYaxis()->SetTitleOffset(0.7);
     PValueGraph->GetXaxis()->SetLabelFont(42);
     PValueGraph->GetYaxis()->SetTitleFont(42);
     PValueGraph->GetYaxis()->SetLabelFont(42);
     PValueGraph->SetLineWidth(1);
     PValueGraph->Draw("AP");


     PValueGraph->SetTitle("");

     pt = new TPaveText(0,0,0,0,"blNDC");
     pt->SetName("title");
     pt->SetBorderSize(2);
     pt->SetFillColor(kWhite);
     //TText *text12 = pt->AddText("");
     pt->Draw();

     PValueGraph->GetYaxis()->SetRangeUser(yminPValue,ymaxPValue);
     PValueGraph->GetXaxis()->SetLimits(xminPValue,xmaxPValue);

     c2_2->Modified();
     c2->cd();
     c2->Modified();
     c2->cd();
     c2->SetSelected(c2);

     nomFichier = "JanChi2PValueVsVar";
     enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c2);

     c2->Clear();
     //c2->Delete();
     Janchi2->Delete();
     PValueGraph->Delete();




}
 
