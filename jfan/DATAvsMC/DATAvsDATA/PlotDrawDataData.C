#include <TLorentzVector.h>
#include <TVector3.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFormula.h>
#include <TF1.h>
#include <TF2.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TObject.h>
#include <string.h>
#include <algorithm>
#include <TBranch.h>
#include <TString.h>
#include <TBits.h>
#include <TMath.h>
#include "TROOT.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <TLatex.h>
#include <THStack.h>
#include <TLegendEntry.h>
#include <TMinuit.h>
#include <TPaveStats.h>


#include "CMSStyle.C"



using namespace std;

double fonction_affine(double *x, double *par){
	return x[0]*par[0] - x[1]*par[1];
}

void DrawDataDataplot(TTree *Data_miniTree, TTree *OLD_Data_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c0){

	CMSstyle();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

  // Get Histo_Data from eventTree
  TH1F *Histo_Data_temp = new TH1F();
  string variable_Data = var + ">>Histo_Data_temp"+limits;
  Data_miniTree->Draw(variable_Data.c_str(), cut.c_str());
  TH1F *Histo_Data = (TH1F*)gDirectory->Get("Histo_Data_temp");
  //c1->Clear();
/*
  // Get Histo_DYToMuMu from eventTree
  TH1F *Histo_DYToMuMu_temp = new TH1F();
  string variable_DYToMuMu = var + ">>Histo_DYToMuMu_temp";
  DYToMuMu_miniTree->Draw(variable_DYToMuMu.c_str(), cut.c_str());
  TH1F *Histo_DYToMuMu = (TH1F*)gDirectory->Get("Histo_DYToMuMu_temp");
  c1->Clear();


  // Get Histo_OLD_DYToMuMu from eventTree
  TH1F *Histo_OLD_DYToMuMu_temp = new TH1F();
  string variable_OLD_DYToMuMu = var + ">>Histo_OLD_DYToMuMu_temp";
  OLD_DYToMuMu_miniTree->Draw(variable_OLD_DYToMuMu.c_str(), cut.c_str());
  TH1F *Histo_OLD_DYToMuMu = (TH1F*)gDirectory->Get("Histo_OLD_DYToMuMu_temp");
  c1->Clear();
*/

  // Get Histo_OLD_Data from eventTree
  TH1F *Histo_OLD_Data_temp = new TH1F();
  string variable_OLD_Data = var + ">>Histo_OLD_Data_temp"+limits;
  OLD_Data_miniTree->Draw(variable_OLD_Data.c_str(), cut.c_str());
  TH1F *Histo_OLD_Data = (TH1F*)gDirectory->Get("Histo_OLD_Data_temp");
  //c1->Clear();
  // Get Histo_WJetsToLNu from eventTree
//  TH1F *Histo_WJetsToLNu_temp = new TH1F();
//  string variable_WJetsToLNu = var + ">>Histo_WJetsToLNu_temp";
//  WJetsToLNu_miniTree->Draw(variable_WJetsToLNu.c_str(), cut.c_str());
//  TH1F *Histo_WJetsToLNu = (TH1F*)gDirectory->Get("Histo_WJetsToLNu_temp");
//  c1->Clear();

  // Get Histo_QCDMu from eventTree
//  TH1F *Histo_QCDMu_temp = new TH1F();
//  string variable_QCDMu = var + ">>Histo_QCDMu_temp";
//  QCDMu_miniTree->Draw(variable_QCDMu.c_str(), cut.c_str());
//  TH1F *Histo_QCDMu = (TH1F*)gDirectory->Get("Histo_QCDMu_temp");
//  c1->Clear();

  // Get the number of entries for further normalization
//  double a = Histo_Data->Integral();
/*
  double b_DYToMuMu = Histo_DYToMuMu->Integral();
  if( (a==0.0) || (b_DYToMuMu==0.0) ){
    cout << "no entries to plots" <<endl;
    return; 
  }*/
  // Normalize
  Histo_Data->Sumw2(); // In order to have the correct error bars on data after renormalization
  Histo_OLD_Data->Sumw2(); // In order to have the correct error bars on data after renormalization
  // // Normalize MC and Data to 1
  //Histo_Data->Scale((double)((double)1.0/(double)a));
  //Histo_MC->Scale((double)((double)1.0/(double)b));
  // // Normalize MC to Data number of entries
//  double integratedLuminosity = 191.09326;

//  double XSectionDYToMuMu = 1614.0;
//  double XSectionOLD_DYToMuMu = 1614.0;
//  double XSectionOLD_Data = 121.0;
//	double XSectionWJetsToLNu = 24640.0;
//	double XSectionQCDMu = 84679.30;

//  double InitialNumberDYToMuMu = 1995236.0;
//  double InitialNumberOLD_DYToMuMu = 1995236.0;
//  double InitialNumberOLD_Data = 1164208.0;
//	double InitialNumberWJetsToLNu = 15110974.0;
//	double InitialNumberQCDMu = 29434562.0;

// Normalize everything to 1
	double N_Data = Histo_Data->Integral();
	double N_OLD_Data = Histo_OLD_Data->Integral();
//	double N_DYToMuMu = Histo_DYToMuMu->Integral();
//	double N_OLD_DYToMuMu = Histo_OLD_DYToMuMu->Integral();

	Histo_Data->Scale((double)((double)1.0/(double)N_Data));
	Histo_OLD_Data->Scale((double)((double)1.0/(double)N_OLD_Data));
//	Histo_DYToMuMu->Scale((double)((double)1.0/(double)N_DYToMuMu));
//	Histo_OLD_DYToMuMu->Scale((double)((double)1.0/(double)N_OLD_DYToMuMu));


//  Histo_DYToMuMu->Scale((double)(  (double)((double)(XSectionDYToMuMu) / (double)(InitialNumberDYToMuMu)) * (double)integratedLuminosity));
//  Histo_OLD_DYToMuMu->Scale((double)(  (double)((double)(XSectionOLD_DYToMuMu) / (double)(InitialNumberOLD_DYToMuMu)) * (double)integratedLuminosity));
//  Histo_OLD_Data->Scale((double)(  (double)((double)(XSectionOLD_Data) / (double)(InitialNumberOLD_Data)) * (double)integratedLuminosity));
//  Histo_WJetsToLNu->Scale((double)(  (double)((double)(XSectionWJetsToLNu) / (double)(InitialNumberWJetsToLNu)) * (double)integratedLuminosity));
//  Histo_QCDMu->Scale((double)(  (double)((double)(XSectionQCDMu) / (double)(InitialNumberQCDMu)) * (double)integratedLuminosity));
  // Adding histograms for binned samples
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt30);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt80);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt170);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt300);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt470);

//	Histo_WJetsToLNu->Add(Histo_QCDMu);
//	Histo_OLD_Data->Add(Histo_WJetsToLNu);
//	Histo_OLD_DYToMuMu->Add(Histo_OLD_Data);
//	Histo_DYToMuMu->Add(Histo_OLD_DYToMuMu);

	// Total MC histo for comupting min/max
//	TH1F *Histo_allMC = new TH1F(*Histo_QCD_Mu_Pt20to30);
//	Histo_allMC->Add(Histo_QCD_Pt15);
//	Histo_allMC->Add(Histo_InclusiveMu15);
//	Histo_allMC->Add(Histo_ZmumuJet_Pt0to15);
//	Histo_allMC->Add(Histo_ZJets_7TeV);
//	Histo_allMC->Add(Histo_WJets_7TeV);
//	Histo_allMC->Add(Histo_TTbarJets_Tauola);
//	Histo_allMC->Add(Histo_DYToMuMu);


  // Get the maxs and the mins to further correct the Y-axis
  double dataMax = Histo_Data->GetMaximum();
  double YMax = dataMax;

//  double DYToMuMuMax = Histo_DYToMuMu->GetMaximum();
//  YMax = max(YMax, DYToMuMuMax);

	double OLD_dataMax = Histo_OLD_Data->GetMaximum();
	YMax = max(YMax, OLD_dataMax);

//	double OLD_DYToMuMuMax = Histo_OLD_DYToMuMu->GetMaximum();
//  YMax = max(YMax, OLD_DYToMuMuMax);

//	double allMCMax = Histo_allMC->GetMaximum();
//	YMax = max(YMax, allMCMax);

  double dataMin = YMax;
//  double OLD_DYToMuMuMin = YMax;
//  double DYToMuMuMin = YMax;
  double OLD_DataMin = YMax;
//  double WJetsToLNuMin = YMax;
//  double QCDMuMin = YMax;

	double allMCMin = YMax;

  double YMin = YMax;

  // Gets the actual minimum for each histogram, and not the unfilled bin if any

  for( int ibin=1 ; ibin<Histo_Data->GetNbinsX() ; ibin++ ){
		if( ((Histo_Data->GetBinContent(ibin))!=0) ){
			YMax = max(YMax, (Histo_Data->GetBinContent(ibin) + Histo_Data->GetBinError(ibin)));
//			cout << "YMax= " << YMax << endl;
		}
//		cout << "ibin= " << ibin << "\tcontent= " << Histo_Data->GetBinContent(ibin) << "\terror= " << Histo_Data->GetBinError(ibin) << endl;
    if( ((Histo_Data->GetBinContent(ibin))!=0) && ((Histo_Data->GetBinContent(ibin))<dataMin) ){
      dataMin = Histo_Data->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, dataMin);
/*
  for( int ibin=1 ; ibin<Histo_DYToMuMu->GetNbinsX() ; ibin++ ){
    if( ((Histo_DYToMuMu->GetBinContent(ibin))!=0) && ((Histo_DYToMuMu->GetBinContent(ibin))<DYToMuMuMin) ){
      DYToMuMuMin = Histo_DYToMuMu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, DYToMuMuMin);

 for( int ibin=1 ; ibin<Histo_OLD_DYToMuMu->GetNbinsX() ; ibin++ ){
    if( ((Histo_OLD_DYToMuMu->GetBinContent(ibin))!=0) && ((Histo_OLD_DYToMuMu->GetBinContent(ibin))<OLD_DYToMuMuMin) ){
      OLD_DYToMuMuMin = Histo_OLD_DYToMuMu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, OLD_DYToMuMuMin);
*/
  for( int ibin=1 ; ibin<Histo_OLD_Data->GetNbinsX() ; ibin++ ){
    if( ((Histo_OLD_Data->GetBinContent(ibin))!=0) && ((Histo_OLD_Data->GetBinContent(ibin))<OLD_DataMin) ){
      OLD_DataMin = Histo_OLD_Data->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, OLD_DataMin);

//  for( int ibin=1 ; ibin<Histo_WJetsToLNu->GetNbinsX() ; ibin++ ){
//    if( ((Histo_WJetsToLNu->GetBinContent(ibin))!=0) && ((Histo_WJetsToLNu->GetBinContent(ibin))<WJetsToLNuMin) ){
//      WJetsToLNuMin = Histo_WJetsToLNu->GetBinContent(ibin);
//    }
//  }
//  YMin = min(YMin, WJetsToLNuMin);

//  for( int ibin=1 ; ibin<Histo_QCDMu->GetNbinsX() ; ibin++ ){
//    if( ((Histo_QCDMu->GetBinContent(ibin))!=0) && ((Histo_QCDMu->GetBinContent(ibin))<QCDMuMin) ){
//      QCDMuMin = Histo_QCDMu->GetBinContent(ibin);
//    }
//  }
//  YMin = min(YMin, QCDMuMin);


//  cout << "YMax= "<< YMax << "\t\tYMin= " << YMin << endl;
  double YMin_lin = (double)YMin / (double)10.0;
//  double Range_lin = ((double)(YMax - YMin_lin)) / ((double)(0.8));
  double Range_lin = ((double)(YMax - YMin_lin)) / ((double)(1.0));
  double YMax_lin = 0.2*Range_lin + YMax;
/*
  double Range_lin = ((double)(YMax - YMin)) / ((double)(0.77));
  double YMax_lin = 0.2*Range_lin + YMax;
  double YMin_lin = max(YMin - 0.03*Range_lin, (double)YMin / (double)10.0);
*/
  double Range_log = ((double)(log10(YMax) - log10(YMin))) / ((double)(0.77));
//  cout << "Range_lin= " << Range_lin << "\t\tRange_log= " << Range_log << endl;
  double YMax_log = pow(10.0, 0.2*Range_log + log10(YMax));
  double YMin_log = pow(10.0, log10(YMin) - 0.03*Range_log);
//  cout << "YMin_lin= " << YMin_lin << "\t\tYMax_lin= " << YMax_lin << endl;
//  cout << "YMin_log= " << YMin_log << "\t\tYMax_log= " << YMax_log << endl;




	c0->Divide(1,2);
	c0->cd(1);
//	gPad->SetNumber(1);
	gPad->SetPad(0,0.2,1,1);
//	gPad->SetBottomMargin(0);
	gPad->Draw();                                                             

	c0->cd(2);
//	gPad->SetNumber(2);
	gPad->SetPad(0,0.,1,0.2);
//	gPad->SetTopMargin(0);
//	gPad->SetBottomMargin(0.3);
	gPad->Draw();                                                   


	

/*
TPad *pad =new TPad("haut","haut",0,0.4,1,1);
	    pad->SetNumber(1);
//	    cout << pad->GetBottomMargin() << endl;
//	    pad->SetBottomMargin(0);
	    pad->Draw();
	    
	    TPad *pad2 =new TPad("milieu","milieu",0,0.2,1,0.4);
	    pad2->SetNumber(2);
//	    pad2->SetTopMargin(0);
//	   pad2->SetBottomMargin(0.3);
	    pad2->Draw();

		TPad *pad3=new TPad("bas", "bas", 0, 0, 1, 0.2);
		pad3->SetNumber(3);
		pad3->Draw();
*/

  c0->cd(1);
  // Setup the histo and canvas names and title
  string data_name = "Data_" + pic + "_" + name;
  string mc_name = "MC_" + pic + "_" + name;
  string canvas_name = "DataData_" + pic + "_" + name;
  std::ostringstream binWidthOSS;
  binWidthOSS << (double)Histo_Data->GetBinWidth(1);
  string binWidth = binWidthOSS.str();
  string YaxisTitle = "";
  if((Title.rfind("[") < Title.size()) && (Title.rfind("]") < Title.size())){
//    string unit = Title.substr(Title.rfind("[")+1, Title.size()-Title.rfind("]")-2);
    string unit = Title.substr(Title.rfind("[")+1, Title.rfind("]")-Title.rfind("[")-1);
    YaxisTitle = "a.u. / " + binWidth + " " + unit;
  } else {
    YaxisTitle = "a.u. / " + binWidth;
  }
  Histo_Data->SetName(data_name.c_str());
//  Histo_QCDMu->SetName(mc_name.c_str());
//	Histo_DYToMuMu->SetName(mc_name.c_str());
	
  c0->SetName(canvas_name.c_str());
  c0->SetTitle(canvas_name.c_str());

  // Draw the comparison plots
//	// Template empty histo
//	TH1F *Histo_template = new TH1F("Histo_template", "Histo_template", Histo_Data->GetNbinsX(), Histo_Data->GetXaxis()->GetXmin(),  Histo_Data->GetXaxis()->GetXmax());
//	Histo_template->SetAxisRange(Histo_Data->GetXaxis()->GetXmin(),  Histo_Data->GetXaxis()->GetXmax(), "X");
//	Histo_template->SetAxisRange(YMin_lin, YMax_lin,"Y");
//	Histo_template->SetMaximum(YMax_lin);
//	Histo_template->SetMinimum(YMin_lin);
//	Histo_template->Draw();
//	c0->Update();
//	c0->Draw();

  TLegend *legend = new TLegend(0.55, 0.82, 0.8, 0.93, "");
  legend->SetTextSize(0.05);
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetShadowColor(kWhite);
  legend->AddEntry(Histo_Data->GetName(), "Data Run2012 52X", "lp");
  legend->AddEntry(Histo_OLD_Data->GetName(), "Data Run2012 53X", "l");
//  legend->AddEntry(Histo_DYToMuMu->GetName(), "Z#mu#muJet 41X", "f");
//  legend->AddEntry(Histo_OLD_DYToMuMu->GetName(), "Z#mu#muJets 39X", "f");
//  legend->AddEntry(Histo_WJetsToLNu->GetName(), "WJets", "f");
//  legend->AddEntry(Histo_QCDMu->GetName(), "QCD #mu", "f");
//  legend->AddEntry(Histo_DYToMuMu->GetName(), "PhotonJet", "f");
//  legend->AddEntry(Histo_QCD_Mu_Pt20to30->GetName(), "QCD Mu", "f");
//  legend->AddEntry(Histo_ZJets_7TeV->GetName(), "ZJets madgraph", "f");



  // // First: draw the data to get correct Y-axis scale
  gPad->Update();
  Histo_Data->GetXaxis()->SetTitle(Title.c_str());
  Histo_Data->GetYaxis()->SetTitle(YaxisTitle.c_str());
  Histo_Data->SetTitleSize(0.05,"XY");;
  Histo_Data->SetLineColor(kBlack);
  Histo_Data->SetMarkerColor(kBlack);
  Histo_Data->SetMarkerSize(0.7);
  Histo_Data->SetMarkerStyle(20);
  Histo_Data->SetMaximum(YMax_lin);
  Histo_Data->SetMinimum(YMin_lin);
//  Histo_Data->GetYaxis()->SetRangeUser(YMin_lin, YMax_lin);
//  Histo_Data->Draw("E1sames");

  // // Second: draw MC on the same canvas
//  Histo_InclusiveMu15->SetLineColor(kBlack);
//  Histo_InclusiveMu15->SetFillColor(kGreen-6);
//  Histo_InclusiveMu15->SetFillStyle(3001);
//  Histo_InclusiveMu15->SetMaximum(YMax_lin);
//  Histo_InclusiveMu15->SetMinimum(YMin_lin);
//  Histo_InclusiveMu15->Draw("same");  
//  OLD_DataAddEntry(Histo_InclusiveMu15->GetName(), "InclusiveMu15", "f");

//  Histo_QCDMu->SetLineColor(kBlack);
//  Histo_QCDMu->SetFillColor(kGreen-6);
//  Histo_QCDMu->SetFillStyle(3001);
//  Histo_QCDMu->SetMaximum(YMax_lin);
//  Histo_QCDMu->SetMinimum(YMin_lin);

//  Histo_WJetsToLNu->SetLineColor(kBlack);
//  Histo_WJetsToLNu->SetFillColor(kMagenta+3);
//  Histo_WJetsToLNu->SetFillStyle(3001);
//  Histo_WJetsToLNu->SetMaximum(YMax_lin);
//  Histo_WJetsToLNu->SetMinimum(YMin_lin);
  Histo_Data->Draw("E1");                                             

  Histo_OLD_Data->SetLineColor(kRed);
  Histo_OLD_Data->SetMarkerColor(kRed);
  Histo_OLD_Data->SetMarkerSize(0.7);
  Histo_OLD_Data->SetMarkerStyle(20);
//  Histo_OLD_Data->SetFillColor(kBlue);
//  Histo_OLD_Data->SetFillStyle(3001);
  Histo_OLD_Data->SetMaximum(YMax_lin);
  Histo_OLD_Data->SetMinimum(YMin_lin);
  Histo_OLD_Data->Draw("same");                                           

//  Histo_DYToMuMu->SetLineColor(kBlack);
//  Histo_DYToMuMu->SetFillColor(kMagenta);
//  Histo_DYToMuMu->SetFillStyle(3001);
//  Histo_DYToMuMu->SetMaximum(YMax_lin);
//  Histo_DYToMuMu->SetMinimum(YMin_lin);

/*
  Histo_DYToMuMu->SetLineColor(kBlack);
  Histo_DYToMuMu->SetFillColor(kRed);
  Histo_DYToMuMu->SetFillStyle(3001);
  Histo_DYToMuMu->SetMaximum(YMax_lin);
  Histo_DYToMuMu->SetMinimum(YMin_lin);
  Histo_DYToMuMu->GetXaxis()->SetTitle(Title.c_str());
  Histo_DYToMuMu->GetYaxis()->SetTitle(YaxisTitle.c_str());
  Histo_DYToMuMu->Draw("HISTsame");

	Histo_OLD_DYToMuMu->SetLineColor(kBlack);
  Histo_OLD_DYToMuMu->SetFillColor(kOrange);
  Histo_OLD_DYToMuMu->SetFillStyle(3001);
  Histo_OLD_DYToMuMu->SetMaximum(YMax_lin);
  Histo_OLD_DYToMuMu->SetMinimum(YMin_lin);
  Histo_OLD_DYToMuMu->Draw("HISTsame");
*/
//	Histo_QCD_Mu_Pt20to30->Draw("same");
//  Histo_WJetsToLNu->Draw("same");
//  Histo_QCDMu->Draw("same");



  // // Third: re-draw Data so that data appears in front of MC
  //Histo_Data->Draw("PE1 same");                                          
  //Histo_OLD_Data->Draw("same");

  // // Fourth: redraw axis so that axis appears in front of everything
  gPad->RedrawAxis();

  // // Fifth: draw legend
  legend->Draw("same");
	c0->Update();

	c0->cd(2);
	TH1F *Histo_DataRatio = (TH1F*) Histo_Data->Clone("Dataratio");
	Histo_DataRatio->Reset();
//	Histo_DataRatio->Sumw2();
	Histo_DataRatio->Divide(Histo_Data, Histo_OLD_Data, 1.0, 1.0);
	Histo_DataRatio->SetMinimum(0.0);
	Histo_DataRatio->SetMaximum(2);
	Histo_DataRatio->GetYaxis()->SetTitle("DATA / OLD DATA");
	Histo_DataRatio->GetXaxis()->SetTitle("");
	Histo_DataRatio->SetTitleOffset(0.35,"Y");
	Histo_DataRatio->SetTitleSize(0.1,"Y");
	Histo_DataRatio->SetLabelSize(0.15,"XY"); 
	Histo_DataRatio->GetYaxis()->CenterTitle();
	Histo_DataRatio->SetNdivisions(509 ,"Y");
	Histo_DataRatio->Draw("EP");                                     
	TLine *l = new TLine(Histo_DataRatio->GetXaxis()->GetXmin(),1.,Histo_DataRatio->GetXaxis()->GetXmax(),1.);
	l->SetLineWidth(1.5); 
	l->Draw("same");
	c0->Update();
/*
	c0->cd(3);
 	TH1F *Histo_DYToMuMuRatio = (TH1F*) Histo_DYToMuMu->Clone("DYToMuMuratio");
	Histo_DYToMuMuRatio->Reset();
//	Histo_DYToMuMuRatio->Sumw2();
	Histo_DYToMuMuRatio->Divide(Histo_DYToMuMu, Histo_OLD_DYToMuMu, 1.0, 1.0);
	Histo_DYToMuMuRatio->SetMinimum(0.0);
	Histo_DYToMuMuRatio->SetMaximum(2);
	Histo_DYToMuMuRatio->GetYaxis()->SetTitle("MC / OLD MC");
	Histo_DYToMuMuRatio->SetTitleOffset(0.35,"Y");
	Histo_DYToMuMuRatio->SetTitleSize(0.11,"Y");
	Histo_DYToMuMuRatio->SetLabelSize(0.1,"Y"); 
	Histo_DYToMuMuRatio->GetYaxis()->CenterTitle();
	Histo_DYToMuMuRatio->SetNdivisions(509 ,"Y");
	Histo_DYToMuMuRatio->Draw("EP");
	TLine *l1 = new TLine(Histo_DYToMuMuRatio->GetXaxis()->GetXmin(),1.,Histo_DYToMuMuRatio->GetXaxis()->GetXmax(),1.);
	l1->SetLineWidth(1.5); 
	l1->Draw("same");
	c0->Update();
*/
	c0->cd(1);
  TLatex latexLabel;
//  std::ostringstream intLumiString;
//  intLumiString << setprecision (2) << fixed << integratedLuminosity;
//  string intLumiText = "#intL= " + intLumiString.str() + " pb^{-1}";
  latexLabel.SetTextSize(0.04);
  latexLabel.SetNDC();
  latexLabel.DrawLatex(0.15, 0.96, "CMS 2012");
  latexLabel.DrawLatex(0.42, 0.96, "#sqrt{s} = 8 TeV");
//  latexLabel.DrawLatex(0.57, 0.96, intLumiText.c_str());

  // // Sixth: update canvas
  c0->Update();
  //c0->Draw();                            

  // Print the canvas
  string PicName="PlotGif/DataData_" + pic + "_" + name + ".gif";
  c0->Print(PicName.c_str());
  PicName="PlotEps/DataData_" + pic + "_" + name + ".eps";
  c0->Print(PicName.c_str());
  string convert = "convert PlotEps/DataData_" + pic + "_" + name + ".eps" + " PlotPdf/DataData_" + pic + "_" + name + ".pdf";
  system(convert.c_str());


  if (inlog==true) {
    c0->cd(1);
    Histo_Data->SetMaximum(YMax_log);
    Histo_Data->SetMinimum(YMin_log);
    Histo_Data->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
//    Histo_DYToMuMu->SetMaximum(YMax_log);
//    Histo_DYToMuMu->SetMinimum(YMin_log);
//    Histo_DYToMuMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
//    Histo_QCD_Pt15->SetMaximum(YMax_log);
//    Histo_QCD_Pt15->SetMinimum(YMin_log);
//    Histo_QCD_Pt15->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    Histo_OLD_Data->SetMaximum(YMax_log);
    Histo_OLD_Data->SetMinimum(YMin_log);
    Histo_OLD_Data->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_WJetsToLNu->SetMaximum(YMax_log);
//    Histo_WJetsToLNu->SetMinimum(YMin_log);
//    Histo_WJetsToLNu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_QCDMu->SetMaximum(YMax_log);
//    Histo_QCDMu->SetMinimum(YMin_log);
//    Histo_QCDMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_InclusiveMu15->SetMaximum(YMax_log);
//    Histo_InclusiveMu15->SetMinimum(YMin_log);
//    Histo_InclusiveMu15->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_OLD_DYToMuMu->SetMaximum(YMax_log);
//    Histo_OLD_DYToMuMu->SetMinimum(YMin_log);
//    Histo_OLD_DYToMuMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    c0->cd(1)->SetLogy(1);
    c0->Update();
    //c0->Draw();                     fan
    string PicName_log="PlotGif/DataData_" + pic + "_" + name + "_log.gif";
    c0->Print(PicName_log.c_str());
    PicName="PlotEps/DataData_" + pic + "_" + name + "_log.eps";
    c0->Print(PicName.c_str());
    string convert = "convert PlotEps/DataData_" + pic + "_" + name + "_log.eps" + " PlotPdf/DataData_" + pic + "_" + name + "_log.pdf";
    system(convert.c_str());
    c0->cd(1)->SetLogy(0);
    c0->Update();
  }


  // Clean the memory
  c0->Clear();       
  legend->Clear();
//	Histo_template->Delete();
    Histo_Data_temp->Delete();    
    Histo_Data->Delete();        
//  Histo_DYToMuMu_temp->Delete();
//  Histo_DYToMuMu->Delete();

//  Histo_OLD_DYToMuMu_temp->Delete();
//  Histo_OLD_DYToMuMu->Delete();

  Histo_OLD_Data_temp->Delete();
  Histo_OLD_Data->Delete();

//  Histo_WJetsToLNu_temp->Delete();
//  Histo_WJetsToLNu->Delete();

//  Histo_QCDMu_temp->Delete();
//  Histo_QCDMu->Delete();

}



void PlotDrawDataData(){
   
        gROOT->ProcessLine(".x setTDRStyle.C");
        string Data = "/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree/OutputMiniTree/miniTree_53X_Run2012AB.root";    
        //string OLD_Data = "/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree/OutputMiniTree/miniTree_52X_Run2012AB.root";
        string OLD_Data = "/sps/cms/jfan/CMSSW_5_3_2_patch4/src/Morgan/IpnTreeProducer/MuGammaMiniTree/OutputMiniTree/miniTree_53X_Run2012PromptC.root";

        TFile *Data_File = new TFile(Data.c_str());
        TTree* Data_miniTree = (TTree*) Data_File->Get("miniTree");
        TTree* Data_miniTree_allmuons = (TTree*) Data_File->Get("miniTree_allmuons");
        TTree* Data_miniTree_allphotons = (TTree*) Data_File->Get("miniTree_allphotons");
/*
        TFile *DYToMuMu_File = new TFile(DYToMuMu.c_str());
        TTree* DYToMuMu_miniTree = (TTree*) DYToMuMu_File->Get("miniTree");
        TTree* DYToMuMu_miniTree_allmuons = (TTree*) DYToMuMu_File->Get("miniTree_allmuons");
        TTree* DYToMuMu_miniTree_allphotons = (TTree*) DYToMuMu_File->Get("miniTree_allphotons");
        TFile *OLD_DYToMuMu_File = new TFile(OLD_DYToMuMu.c_str());
        TTree* OLD_DYToMuMu_miniTree = (TTree*) OLD_DYToMuMu_File->Get("miniTree");
        TTree* OLD_DYToMuMu_miniTree_allmuons = (TTree*) OLD_DYToMuMu_File->Get("miniTree_allmuons");
        TTree* OLD_DYToMuMu_miniTree_allphotons = (TTree*) OLD_DYToMuMu_File->Get("miniTree_allphotons");
*/
        TFile *OLD_Data_File = new TFile(OLD_Data.c_str());
        TTree* OLD_Data_miniTree = (TTree*) OLD_Data_File->Get("miniTree");
        TTree* OLD_Data_miniTree_allmuons = (TTree*) OLD_Data_File->Get("miniTree_allmuons");
        TTree* OLD_Data_miniTree_allphotons = (TTree*) OLD_Data_File->Get("miniTree_allphotons");
//      TFile *WJetsToLNu_File = new TFile(WJetsToLNu.c_str());
//      TTree* WJetsToLNu_miniTree = (TTree*) WJetsToLNu_File->Get("miniTree");
//      TTree* WJetsToLNu_miniTree_allmuons = (TTree*) WJetsToLNu_File->Get("miniTree_allmuons");
//      TTree* WJetsToLNu_miniTree_allphotons = (TTree*) WJetsToLNu_File->Get("miniTree_allphotons");

//      TFile *QCDMu_File = new TFile(QCDMu.c_str());
//      TTree* QCDMu_miniTree = (TTree*) QCDMu_File->Get("miniTree");
//      TTree* QCDMu_miniTree_allmuons = (TTree*) QCDMu_File->Get("miniTree_allmuons");
//      TTree* QCDMu_miniTree_allphotons = (TTree*) QCDMu_File->Get("miniTree_allphotons");


        TCanvas *c0 = new TCanvas("Default", "Default");

//////  DrawDataDataplot(Data_miniTree_allmuons, DYToMuMu_miniTree_allmuons, OLD_DYToMuMu_miniTree_allmuons, QCDMu_miniTree_allmuons, OLD_Data_miniTree_allmuons, WJetsToLNu_miniTree_allmuons, "Ptmumu", "Ptmumu", "(100,0,200)", "isMM", "dimuon", "p_{T}^{#mu#mu} [GeV]", true, false, c0);

        vector<string> set_of_cuts;
        vector<string> name;

//      set_of_cuts.push_back("isMMGCandidate");
//      name.push_back("selected-00-beforeFSRcuts");
/*
        set_of_cuts.push_back("isAfterFSRCut1");
        name.push_back("selected-01");
        set_of_cuts.push_back("isAfterFSRCut2");
        name.push_back("selected-02");
        set_of_cuts.push_back("isAfterFSRCut3");
        name.push_back("selected-03");
        set_of_cuts.push_back("isAfterFSRCut4");
        name.push_back("selected-04");
        set_of_cuts.push_back("isVeryLooseMMG");
        name.push_back("selected-veryloose");
*/
  set_of_cuts.push_back("isLooseMMG");
  name.push_back("selected-loose");
  set_of_cuts.push_back("isLooseMMG && Photon_isEB");
  name.push_back("selected-loose-EB");
  set_of_cuts.push_back("isLooseMMG && Photon_isEB && Photon_r9 < .94");
  name.push_back("selected-loose-EB-lowR9");
  set_of_cuts.push_back("isLooseMMG && Photon_isEB && Photon_r9 > .94");
  name.push_back("selected-loose-EB-highR9");
  set_of_cuts.push_back("isLooseMMG && Photon_isEE");
  name.push_back("selected-loose-EE");
  set_of_cuts.push_back("isLooseMMG && Photon_isEE && Photon_r9 < .95");
  name.push_back("selected-loose-EE-lowR9");
  set_of_cuts.push_back("isLooseMMG && Photon_isEE && Photon_r9 > .95");
  name.push_back("selected-loose-EE-highR9");

  set_of_cuts.push_back("isTightMMG");
  name.push_back("selected-tight");
  set_of_cuts.push_back("isTightMMG && Photon_isEB");
  name.push_back("selected-tight-EB");
  set_of_cuts.push_back("isTightMMG && Photon_isEB && Photon_r9 < .94");
  name.push_back("selected-tight-EB-lowR9");
  set_of_cuts.push_back("isTightMMG && Photon_isEB && Photon_r9 > .94");
  name.push_back("selected-tight-EB-highR9");
  set_of_cuts.push_back("isTightMMG && Photon_isEE");
  name.push_back("selected-tight-EE");
  set_of_cuts.push_back("isTightMMG && Photon_isEE && Photon_r9 < .95");
  name.push_back("selected-tight-EE-lowR9");
  set_of_cuts.push_back("isTightMMG && Photon_isEE && Photon_r9 > .95");
  name.push_back("selected-tight-EE-highR9");

/*
        set_of_cuts.push_back("isLooseMMG && isMultipleCandidate==0");
  name.push_back("selected-loose-nomultiple");
        set_of_cuts.push_back("isLooseMMG && Photon_isEB && isMultipleCandidate==0");
  name.push_back("selected-loose-EB-nomultiple");
        set_of_cuts.push_back("isLooseMMG && Photon_isEE && isMultipleCandidate==0");
  name.push_back("selected-loose-EE-nomultiple");
        set_of_cuts.push_back("isTightMMG && isMultipleCandidate==0");
  name.push_back("selected-tight-nomultiple");
        set_of_cuts.push_back("isTightMMG && Photon_isEB && isMultipleCandidate==0");
  name.push_back("selected-tight-EB-nomultiple");
        set_of_cuts.push_back("isTightMMG && Photon_isEE && isMultipleCandidate==0");
  name.push_back("selected-tight-EE-nomultiple");
*/

        //for(int i=0; i<set_of_cuts.size() ; i++)
        for(int i=0; i<5 ; i++)
        {   
//      DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Ptmumu", "Ptmumu", "(100,0,200)", set_of_cuts[i], name[i], "p_{T}^{#mu#mu} [GeV]", true, false, c0);
/*

//              DrawDataDataplot_TH1I(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "NbMuons", "NbMuons", "(10,0,10)", set_of_cuts[i], name[i], "# of muons", true, false, c0);
//              DrawDataDataplot_TH1I(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "NbPhotons", "NbPhotons", "(10,0,10)", set_of_cuts[i], name[i], "# of photons", true, false, c0);

*/
/*

//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Pt", "MuonL_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{leading} [GeV]", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Pt", "MuonS_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{trailing} [GeV]", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Eta", "MuonL_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{leading}", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Eta", "MuonS_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{trailing}", true, false, c0);

//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Phi", "MuonL_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{leading}", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Phi", "MuonS_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{trailing}", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Charge", "MuonL_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{leading}", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Charge", "MuonS_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{trailing}", true, false, c0);

//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonF_Pt", "MuonF_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{far} [GeV]", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonN_Pt", "MuonN_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{near} [GeV]", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonF_Eta", "MuonF_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{far}", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonN_Eta", "MuonN_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{near}", true, false, c0);

                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "MuonF_Phi", "MuonF_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{far}", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "MuonN_Phi", "MuonN_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{near}", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "MuonF_Charge", "MuonF_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{far}", true, false, c0);

                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "MuonN_Charge", "MuonN_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{near}", true, false, c0);
*/
              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Mmumu", "Mmumu", "(150,0,300)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Mmumu", "Mmumu_extended", "(30,30,90)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", true, false, c0);








                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Mmumugamma", "Mmumugamma", "(120,60,120)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma} [GeV]", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Mmumugamma", "Mmumugamma_extended", "(60,60,120)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma} [GeV]", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Mmumugamma", "Mmumugamma_zoom", "(48,85,97)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma} [GeV]", true, false, c0);

                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Mmumugamma_SCraw", "Mmumugamma_SCraw", "(120,60,120)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma}^{SCraw} [GeV]", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Mmumugamma_SCraw_fEta", "Mmumugamma_SCraw_fEta", "(50,0,400)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma}^{SCraw x fEta} [GeV]", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "deltaRNear", "deltaRNear", "(100,0,1)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{near})", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "deltaRFar", "deltaRFar", "(100,0,5)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{far})", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "deltaRPlus", "deltaRPlus", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{plus})", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "deltaRMinus", "deltaRMinus", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{minus})", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "deltaRLeading", "deltaRLeading", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{leading})", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "deltaRSubleading", "deltaRSubleading", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{trailing})", true, false, c0);




//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_Eta", "Photon_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta^{#gamma}", true, false, c0);







                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_Eta", "Photon_Eta", "(16,-3.2,3.2)", set_of_cuts[i], name[i], "#eta^{#gamma}", true, false, c0);

                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_Phi", "Photon_Phi", "(21,-3.15,3.15)", set_of_cuts[i], name[i], "#phi^{#gamma}", true, false, c0);







/*
*/
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_E", "Photon_E", "(100, 0, 100)", set_of_cuts[i], name[i], "E^{#gamma} [GeV]", true, false, c0);
//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_E", "Photon_E_extended", "(100, 0, 150)", set_of_cuts[i], name[i], "E^{#gamma} [GeV]", true, false, c0);

     
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_SC_rawE", "Photon_SC_rawE", "(50, 0, 400)", set_of_cuts[i], name[i], "E^{SC raw} [GeV]", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_SC_rawEt", "Photon_SC_rawEt", "(40, 0, 160)", set_of_cuts[i], name[i], "E^{SC raw}_{T} [GeV]", true, false, c0);

                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_E", "Photon_E", "(50, 0, 250)", set_of_cuts[i], name[i], "E^{#gamma} [GeV]", true, false, c0);



//              DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_Et", "Photon_Et", "(100, 0, 100)", set_of_cuts[i], name[i], "E_{T}^{#gamma} [GeV]", true, false, c0);
                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "Photon_Et", "Photon_Et", "(50, 0, 100)", set_of_cuts[i], name[i], "E_{T}^{#gamma} [GeV]", true, false, c0);  









/*
                DrawDataDataplot_TH2F(Data_miniTree, OLD_Data_miniTree, "Mmumu", "Mmumugamma", "(900,0,300,900,0,300)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", "m_{#mu#mu#gamma} [GeV]", "Mmumu_VS_Mmumugamma", false, false, c0);
                DrawDataDataplot_TH2F(Data_miniTree, OLD_Data_miniTree, "Mmumu", "Mmumugamma", "(450,20,100,450,80,110)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", "m_{#mu#mu#gamma} [GeV]", "Mmumu_VS_Mmumugamma_extended", false, false, c0);

DrawDataDataplot_TH2F(Data_miniTree, OLD_Data_miniTree, "((91.1876**2 - Mmumu**2 ) / (Mmumugamma**2 - Mmumu**2))*Photon_E", "Photon_E", "(400,0,100,400,0,100)", set_of_cuts[i], name[i], "E_{true} = k*E_{reco} [GeV]", "E_{reco} [GeV]", "Etrue_VS_Ereco", false, false, c0, true);DrawDataDataplot_TH2F(Data_miniTree, OLD_Data_miniTree, "((91.1876**2 - Mmumu**2 ) / (Mmumugamma**2 - Mmumu**2))*Photon_E", "Photon_E", "(600,0,150,600,0,150)", set_of_cuts[i], name[i], "E_{true} = k*E_{reco} [GeV]", "E_{reco} [GeV]", "Etrue_VS_Ereco_extended", false, false, c0, true);*/



                DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_k", "mmg_k", "(40,0,2.0)", set_of_cuts[i], name[i], "k = E_{muons} / E_{reco}", true, false, c0);
        DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_s", "mmg_s", "(40,-1.0,1.0)", set_of_cuts[i], name[i], "s = E_{reco} / E_{muons} - 1", true, false, c0);

        DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_k", "mmg_k_extended", "(20,0,2.0)", set_of_cuts[i], name[i], "k = E_{muons} / E_{reco}", true, false, c0);
        DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_s", "mmg_s_extended", "(20,-1.0,1.0)", set_of_cuts[i], name[i], "s = E_{reco} / E_{muons} - 1", true, false, c0);
        DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_k", "mmg_k_zoom", "(20,0.6,1.6)", set_of_cuts[i], name[i], "k = E_{muons} / E_{reco}", true, false, c0);
        DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_s", "mmg_s_zoom", "(20,-0.5,0.5)", set_of_cuts[i], name[i], "s = E_{reco} / E_{muons} - 1", true, false, c0);

        DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_k_SCraw", "mmg_k_SCraw", "(60,-1.0,3.0)", set_of_cuts[i], name[i], "k_{SCraw} = E_{muons} / E_{SCraw}", true, false, c0);
        DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_s_SCraw", "mmg_s_SCraw", "(60,-2.0,2.0)", set_of_cuts[i], name[i], "s_{SCraw} = E_{SCraw} / E_{muons} - 1", true, false, c0);




//      DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_k_SCraw_fEta", "mmg_k_SCraw_fEta", "(60,-1.0,3.0)", set_of_cuts[i], name[i], "k_{SCraw x fEta} = E_{muons} / E_{SCraw x fEta}", true, false, c0);
//      DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "mmg_s_SCraw_fEta", "mmg_s_SCraw_fEta", "(60,-2.0,2.0)", set_of_cuts[i], name[i], "s_{SCraw x fEta} = E_{SCraw x fEta} /  E_{muons} -1", true, false, c0);









/*
//DrawDataDataplot_TH2F(Data_miniTree, OLD_Data_miniTree, "Photon_convNTracks", "MuonN_isoR03_sumPt", "(16,0,4,200,0,10)", set_of_cuts[i], name[i], "Photon_convNTracks", "MuonN_isoR03_sumPt", "Photon_convNTracks_VS_MuonN_isoR03_sumPt", false, false, c0, false);
//DrawDataDataplot_TH2F(Data_miniTree, OLD_Data_miniTree, "Photon_convNTracks", "MuonN_isoR03_sumPt", "(12,0,3,200,0,3)", set_of_cuts[i], name[i], "Photon_convNTracks", "MuonN_isoR03_sumPt", "Photon_convNTracks_VS_MuonN_isoR03_sumPt_extended", false, false, c0, false);
        //      DrawDataDataplot(Data_miniTree, OLD_Data_miniTree, "MuonN_isoR03_sumPt", "MuonN_isoR03_sumPt", "(100,0,10)", set_of_cuts[i], name[i], "MuonN_isoR03_sumPt", true, false, c0);
*/

        }   
}
 
