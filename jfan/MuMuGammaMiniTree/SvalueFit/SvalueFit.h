//#include "TROOT.h "
#ifndef __CINT__
#endif
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
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TGraph.h"
#include "TText.h"
#include "TLine.h"
#include "TGaxis.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include "TROOT.h"
#include "TRint.h"
#include "TDirectory.h"
#include "CrystalBall.C"
#include "setTDRStyle.C"
#include "CMSStyle.C"
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
//#include "PValue.h"
#pragma optimize 0

using namespace RooFit;
using namespace std;

string DoubleToString(double x);

Double_t myfunction(Double_t *x, Double_t *par);
Double_t myfunction2(Double_t *x, Double_t *par);
void myfunc(int EndCaps);
void myfit(TH1D *h1);

float PtCor(float brem, int EndCaps, int ParametersVersion); //ParametersVersion = 1,2,3,4 (=electrons old, new, photons old, new);
Double_t effSigma(TH1 * hist);
void ChaineChi2(char * buffer, double chi2);
double SigmaR(TF1* crystalBall, double Xmin, double Xmax);
double SigmaL(TF1* crystalBall, double Xmin, double Xmax);
float fEta(float eta);
void CrystalBallMethode(double * mean_value, double * mean_error, double * sigma_value, double * ChiSquare, double * DegreesOfFreedom, double param[5], TH1D Data);

void RooCrystalBall(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooLogNormal(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooGamma2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooBifurcatedGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooSumGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooGenericPDF(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier);

void RooLandau2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooLandauConvGaussian(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooKernel(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooVoigtian2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1);

void RangeEstimator(double pourcentage, double centralValue, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax);

void RangeEstimator2(double pourcentage, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax);

void RangeEstimator3(double pourcentage, TChain * chain, TString temp, int Endcaps, double * MinRange, double * MaxRange);

void SymetricRangeEstimator(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp);

void SymetricRangeEstimator2(TChain * chain, double lastX, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp);

void SymetricRangeEstimator3(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp);

Double_t chiSquare(RooPlot* plot_, char* pdfname, char* histname, int nFitParam, double* JanChi2, double* DegreesOfFreedom, double* pValue, int* fewBins);

RooHist* residHist(RooPlot* plot_, char *histname, char* curvename, bool normalize, string dossierSauvegardePull, int iteration);

RooHist* pullHist(RooPlot* plot_, char* histname, char* pdfname, string dossierSauvegardePull, int iteration) { return residHist(plot_, histname, pdfname, true, dossierSauvegardePull, iteration); }


int NbLignesFichier(string fichier);

string DoubleToString(double x)
{

        std::string s;
        {
                std::ostringstream oss; 
                oss << x;
                s = oss.str();
        }
        std::cout << "x = " << x << " s = " << s << std::endl;

        return s;
}


Double_t myfunction(Double_t *x, Double_t *par)
{
        float bremLowThr  = 1.10000002384185791015625;
        float bremHighThr = 8.0;
        if ( x[0] < bremLowThr  ) x[0] = bremLowThr;
        if ( x[0] > bremHighThr ) x[0] = bremHighThr;

        float threshold = par[4];

        float y = par[0]*threshold*threshold + par[1]*threshold + par[2];
        float yprime = 2*par[0]*threshold + par[1];
        float a = par[3];
        float b = yprime - 2*a*threshold;
        float c = y - a*threshold*threshold - b*threshold;

        float bremCorrection = 1.0;
        if ( x[0] < threshold ) bremCorrection = par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
        else bremCorrection = a*x[0]*x[0] + b*x[0] + c;

        return bremCorrection;

}

Double_t myfunction2(Double_t *x, Double_t *par)
{
        float bremLowThr  = 0.89999997615814208984375;
  	float bremHighThr = 6.5;
	if ( x[0] < bremLowThr  ) x[0] = bremLowThr;
        if ( x[0] > bremHighThr ) x[0] = bremHighThr;

        float threshold = par[4];

        float y = par[0]*threshold*threshold + par[1]*threshold + par[2];
        float yprime = 2*par[0]*threshold + par[1];
        float a = par[3];
        float b = yprime - 2*a*threshold;
        float c = y - a*threshold*threshold - b*threshold;

        float bremCorrection = 1.0;
        if ( x[0] < threshold ) bremCorrection = par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
        else bremCorrection = a*x[0]*x[0] + b*x[0] + c;

        return bremCorrection;

}


void myfunc(int Endcaps)
{

	TF1 *f1;
      if(Endcaps == 0) 
      {	
		TF1 *f1 = new TF1("myfunc",myfunction,0,100,5);
		//f1->Draw();
      }
      if(Endcaps == 1) 
      {		
		TF1 *f1 = new TF1("myfunc",myfunction2,0,100,5);
      		//f1->Draw();
      }
	f1->Draw();

	//f1->SetParameters(-0.002362,0.004314,1.001,0.0003413,3.124);
      //f1->SetParNames("constant","coefficient");
}

void myfit(TH1D *h1)
{


      //TH1F *h1=new TH1F("h1","test",100,0,12);
      //h1->FillRandom("myfunc",20000);
      TF1 *f1=(TF1*)gROOT->GetFunction("myfunc");
      //f1->SetParameters(800,1);
      f1->SetFillColor(19);
      f1->SetFillStyle(0);
      f1->SetLineWidth(3);


	TPaveStats *ptstats = new TPaveStats(0.6845638,0.7226277,0.9865772,0.9912587,"brNDC");
        ptstats->SetName("stats");
        ptstats->SetBorderSize(2);
        ptstats->SetFillColor(kWhite);
        ptstats->SetTextColor(19);
        ptstats->SetTextAlign(12);
        ptstats->SetTextFont(42);
        ptstats->SetOptStat(0);
        ptstats->SetOptFit(11111111);
        ptstats->Draw();


      h1->Fit("myfunc");
}




float PtCor(float brem, int EndCaps, int ParametersVersion) 
{
  // brem == phiWidth/etaWidth of the SuperCluster 
  // e  == energy of the SuperCluster 
  // first parabola (for br > threshold) 
  // p0 + p1*x + p2*x^2 
  // second parabola (for br <= threshold) 
  // ax^2 + bx + c, make y and y' the same in threshold 
  // y = p0 + p1*threshold + p2*threshold^2  
  // yprime = p1 + 2*p2*threshold 
  // a = p3 
  // b = yprime - 2*a*threshold 
  // c = y - a*threshold^2 - b*threshold 
  float p0 = 0;
  float p1 = 0;
  float p2 = 0;
  float p3 = 0;
  float p4 = 0;

/*
  int offset;
  if ( EndCaps == 0 ) offset = 0; //Barrel
  else if ( EndCaps == 1 ) offset = 20; //End caps

  else {
    // not supported, produce no correction
    return 1.0;
  }
*/

  //Make No Corrections if brem is invalid! 
	//cout<<endl<<setprecision( 10 )<<"brem in function= "<<brem<<endl;
  if ( brem == 0 ) return 1.0;

 if(EndCaps == 1) //End caps
 {
  float bremLowThr  = 0.89999997615814208984375;
  float bremHighThr = 6.5;
  if ( brem < bremLowThr  ) brem = bremLowThr;
  if ( brem > bremHighThr ) brem = bremHighThr;

  if(ParametersVersion == 1)
  {
  	p0 = -0.07945;
	p1 = 0.1298;
	p2 = 0.9147;
	p3 = -0.001565;
	p4 = 0.9;
  }

  if(ParametersVersion == 2)
  {
  	p0 = -0.12139999866485595703125;
	p1 = 0.2362000048160552978515625;
	p2 = 0.884700000286102294921875;
	p3 = -0.00193000002764165401458740234375;
	p4 = 1.05700004100799560546875;

  }

 if(ParametersVersion == 3)
  {
	p0 = -0.0728;
	p1 = 0.1323;
	p2 = 0.9201;
	p3 = 0.002501;
	p4 = 1.118;
  }

 if(ParametersVersion == 4)
  {
	p0 = -0.07667;
	p1 = 0.1407;
	p2 = 0.9157;
	p3 = 0.00251;
	p4 = 1.117;
  }



 }

if(EndCaps == 0) //Barrel
 {
  float bremLowThr  = 1.10000002384185791015625;
  float bremHighThr = 8.0;
  if ( brem < bremLowThr  ) brem = bremLowThr;
  if ( brem > bremHighThr ) brem = bremHighThr;

	
  if(ParametersVersion == 1)
  {
	p0 = -0.05185;
	p1 = 0.1354;
	p2 = 0.9165;
	p3 = -0.0005626;
	p4 = 1.385;    
  }

  if(ParametersVersion == 2)
  {
	p0 = -0.05289;
	p1 = 0.1374;
	p2 = 0.9141;
	p3 = -0.000669;
	p4 = 1.38;
  }

 if(ParametersVersion == 3)
  {
	p0 = -0.0625;
	p1 = 0.1331;
	p2 = 0.9329;
	p3 = -0.0007823;
	p4 = 1.1;
  }

 if(ParametersVersion == 4)
  {
	p0 = -0.002362;
	p1 = 0.004314;
	p2 = 1.001;
	p3 = 0.0003413;
	p4 = 3.124;
  }

if(ParametersVersion == 5)
  {
        p0 = -0.0004081;
        p1 = -0.005385;
        p2 = 1.008;
        p3 = -0.02381;
        p4 = 123.4;
  }

 }

  float threshold = p4;

  float y = p0*threshold*threshold + p1*threshold + p2;
  float yprime = 2*p0*threshold + p1;
  float a = p3;
  float b = yprime - 2*a*threshold;
  float c = y - a*threshold*threshold - b*threshold;

  float bremCorrection = 1.0;
  if ( brem < threshold ) bremCorrection = p0*brem*brem + p1*brem + p2;
  else bremCorrection = a*brem*brem + b*brem + c;

  return bremCorrection;
}




Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}




void ChaineChi2(char * buffer, double chi2)
{
	
	string chaine1("#chi^{2} / ndf = ");
        string chi2Chaine;
        string chaineFinale;
	
	ostringstream DoubleToString;
        DoubleToString << chi2;
        chi2Chaine = DoubleToString.str();
        chaineFinale = chaine1 + chi2Chaine;

        //size_t size = chaineFinale.size() + 1; 
        //char * buffer = new char[ size ];
        //strncpy( buffer, chaineFinale.c_str(), size );
	strncpy( buffer, chaineFinale.c_str(), 25 );

}

double SigmaR(TF1* crystalBall, double Xmin, double Xmax)
{

	double y = crystalBall->GetMaximum(Xmin, Xmax) * 1.0/exp(1.0);
	double MaxX = crystalBall->GetMaximumX(Xmin, Xmax);
	double sigmaR = crystalBall->GetX(y, MaxX, Xmax) - MaxX;

	return sigmaR;
}

double SigmaL(TF1* crystalBall, double Xmin, double Xmax)
{

        double y = crystalBall->GetMaximum(Xmin, Xmax) * 1.0/exp(1.0);
        double MaxX = crystalBall->GetMaximumX(Xmin, Xmax);
	double sigmaL = MaxX - crystalBall->GetX(y, Xmin, MaxX);

        return sigmaL;
}

float fEta(float eta)
{

  
  float ieta = fabs(eta)*(5/0.087);
  float p0 = 40.2198;  // should be 40.2198
  float p1 = -3.03103e-6;  // should be -3.03103e-6
  float feta = 1;

  if ( ieta < p0 || fabs(eta) > 1.4442 ) feta = 1.0;
  else feta = 1.0/(1.0 + p1*(ieta-p0)*(ieta-p0));

  //correctedEnergy = energy/(1.0 + p1*(ieta-p0)*(ieta-p0));
  return feta;

}


void CrystalBallMethode(double * mean_value, double * mean_error, double * sigma_value, double * ChiSquare, double * DegreesOfFreedom, double param[5], TH1D Data)
{
	TCanvas *cf = new TCanvas("cf", "cf",0,0,600,600);

	int b0 = 0;
        int b1 = 0;
        int b3 = 0;
        int b4 = 0;
        double fit0 = 10;
        double fit1 = 0.6;
        double fit2 = 5;
        double fit3 = 1;
        double fit4 = 0.02;
	
	double ChiSquareBest = 10000000;
	double mean_errorBest = 10000000;
	*mean_error = 10000000;	

	double fit0Temp = 0;
	double fit1Temp = 0;
	double fit2Temp = 0;
	double fit3Temp = 0;
	double fit4Temp = 0;

	double fit0Best = 0;
	double fit1Best = 0;
	double fit2Best = 0;
	double fit3Best = 0;
	double fit4Best = 0;
	
	ChiSquareBest = 10000000;
        mean_errorBest = 10000000;

	for(int r = 0; r<400; r++)
       	{
               	cout<<endl<<"r = "<<r<<endl;
               	//fit0 = r;
		fit0Temp = r;
               	Data.Draw("");
		TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0Temp, fit1, fit2, fit3, fit4);
               	func->SetLineColor(3);
         	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
               	*ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
               	if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
			if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                	{
				//if(mean_error <= mean_errorBest)
				if(*mean_error <0.01)
				{
                       			ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                       			fit0Best = r;
					mean_errorBest = *mean_error;
                       			cout<<endl<<"fit0Best = "<<fit0Best<<endl;
                       			cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
					b0 = 1;
				}
                	}
		}
               	func->Delete();
        }

        //cout<<endl<<"fit0Best = "<<fit0Best<<endl;
        //cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
	if(b0 == 1) fit0 = fit0Best;
        ChiSquareBest = 10000000;	
	mean_errorBest = 10000000;

	for(double s = -1.0; s<1; s+=0.05)
       	{
            	cout<<endl<<"s = "<<s<<endl;
               	//fit1 = s;
		fit1Temp = s;
               	Data.Draw("");
               	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0, fit1Temp, fit2, fit3, fit4);
               	func->SetLineColor(3);
               	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
                *ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
		if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
                        if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                        {
                                //if(mean_error <= mean_errorBest)
                                if(*mean_error <0.01)
				{
                                        ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                                        fit1Best = s;
					mean_errorBest = *mean_error;
                                        cout<<endl<<"fit1Best = "<<fit1Best<<endl;
                                        cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
                                      	b1 = 1;
				}
                        }
		}
               	func->Delete();
        }
        if(b1 == 1) fit1 = fit1Best;
        ChiSquareBest = 10000000;
	mean_errorBest = 10000000;

	for(double t = 0.900; t<1.010; t+=0.001)
        {
               	cout<<endl<<"t = "<<t<<endl;
               	//fit3 = t;
		fit3Temp = t;
               	Data.Draw("");
               	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0, fit1, fit2, fit3Temp, fit4);
               	func->SetLineColor(3);
               	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
                *ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
                if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
			if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                        {
                                //if(mean_error <= mean_errorBest)
                                if(*mean_error <0.01)
				{
                                        ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                                        fit3Best = t;
					mean_errorBest = *mean_error;
                                        cout<<endl<<"fit3Best = "<<fit3Best<<endl;
                                        cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
                                       	b3 = 1;
				}
                        }
		}				
                func->Delete();

        }
        if(b3 == 1) fit3 = fit3Best;
        ChiSquareBest = 10000000;
	mean_errorBest = 10000000;

	for(double u = 0.001; u<0.1; u+=0.001)
        {
               	cout<<endl<<"u = "<<u<<endl;
               	//fit4 = u;
		fit4Temp = u;
               	Data.Draw("");
               	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0, fit1, fit2, fit3, fit4Temp);
               	func->SetLineColor(3);
               	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
                *ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
		if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
                        if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                        {
                	        //if(mean_error <= mean_errorBest)
                                if(*mean_error <0.01)
				{
                                        ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                                        fit4Best = u;
					mean_errorBest = *mean_error;
                                        cout<<endl<<"fit4Best = "<<fit4Best<<endl;
                                        cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
					b4 = 1;
                                }
                        }
		}
                func->Delete();
        }
        if( b4 == 1 ) fit4 = fit4Best;
	cout<<endl<<"fit0Best = "<<fit0Best<<endl<<"fit1Best = "<<fit1Best<<endl<<"fit3Best = "<<fit3Best<<endl<<"fit4Best = "<<fit4Best<<endl;

	cf->Clear();
       	Data.Draw("");
       	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
       	func->SetParameters(fit0, fit1, fit2, fit3, fit4);
	func->SetLineColor(3);
        func->SetLineWidth(3);
       	Data.Fit(func);
       	*mean_value  = func->GetParameter(3);
       	*mean_error  = func->GetParError(3); 
       	*sigma_value = func->GetParameter(4);
	*ChiSquare = func->GetChisquare();
	*DegreesOfFreedom = func->GetNDF();	
	

	param[0] = fit0;
	param[1] = fit1;
	param[2] = fit2;
	param[3] = fit3;
	param[4] = fit4;

	
	cf->Clear();
	func->Delete();
	//cf->Delete();
	cf->Close();
}

void RooCrystalBall(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
	RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");	
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, ""); 
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, ""); 
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG, *isMultipleCandidate, *weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *isLooseMMG, *isMultipleCandidate, *weight_pileUp);
	//RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isLooseMMG, Photon_Et);
        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars); //non weighté
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
	//RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        
/*	
	Erawframe = mmg_s->frame();
	Data_subset->plotOn(Erawframe,Name("myhist"),Binning(300));
        Erawframe->Draw();
*/

	RooRealVar * Eraw_CB_m0 = new RooRealVar("Eraw_CB_m0", "CB #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_CB_sigma = new RooRealVar("Eraw_CB_sigma", "CB ", rms, 0.0, 2 * rms, "GeV");
        RooRealVar * Eraw_CB_alpha = new RooRealVar("Eraw_CB_alpha", "CB #alpha", 0.9, -5.0, 5.0);
        RooRealVar * Eraw_CB_n = new RooRealVar("Eraw_CB_n", "CB n", 10.0, 1.0, 500.0);


/*
	//Ne pas supprimer !!!!//
	RooRealVar Eraw_CB_m0("Eraw_CB_m0", "CB #Delta m_{0}", 1.0, 0.9, 1.05, "GeV");
        RooRealVar Eraw_CB_sigma("Eraw_CB_sigma", "CB ", 0.0215, 0.0, 0.5, "GeV");
        RooRealVar Eraw_CB_alpha("Eraw_CB_alpha", "CB #alpha", 0.9, 0.0, 3.0);
        RooRealVar Eraw_CB_n("Eraw_CB_n", "CB n", 10.0, 0roo.0, 60.0);
*/
        RooCBShape * Eraw_CrystalBall = new RooCBShape("Eraw_CrystalBall","mmg_s_CrystalBall", *mmg_s, *Eraw_CB_m0, *Eraw_CB_sigma, *Eraw_CB_alpha, *Eraw_CB_n);

	int fewBins = 1;
	Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	{
		cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
		Erawframe->Clear();
		Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);		
		Eraw_CrystalBall = new RooCBShape("Eraw_CrystalBall","mmg_s_CrystalBall", *mmg_s, *Eraw_CB_m0, *Eraw_CB_sigma, *Eraw_CB_alpha, *Eraw_CB_n); 

		Erawframe = mmg_s->frame();
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
		//Erawframe->Draw();
        	//Eraw_CrystalBall.fitTo(*Data_subset, Range(0.6,1.2));
		//Eraw_CrystalBall.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
		//Eraw_CrystalBall.fitTo(Data_subset, Range(mean - 5 * rms,mean + 5 * rms));
		//double maxDistri = hh->GetMaximumBin() * 0.01;
		//Eraw_CrystalBall.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
		//Eraw_CrystalBall.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
		Eraw_CrystalBall->fitTo(*Data_subset, Range(RangeMin, RangeMax));
		//RooArgSet* Eraw_CrystalBall_param = Eraw_CrystalBall.getVariables();
		//Eraw_CrystalBall_param->Print("v");
		//Eraw_CrystalBall.plotOn(Erawframe);
        	Eraw_CrystalBall->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


		Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",4, JanChi2, DegreesOfFreedom, pValue, &fewBins);
		if(fewBins == 0) break; 


	}


	//RooDataHist data("data", "dataset with mmg_s",*mmg_s,hh);
	//TF1 * f = Eraw_CrystalBall.asTF( RooArgList(mmg_s) );
	f = Eraw_CrystalBall->asTF( RooArgList(*mmg_s) );
	//ftest = Eraw_CrystalBall.asTF( RooArgList(*mmg_s) );
	//cout<<endl<<"ftest->GetXmax() = "<<ftest->GetXmax()<<endl;
	//cout<<endl<<"ftest->GetMaximum = "<<ftest->GetMaximum(0.0, 2.0,1.E-10,100,false)<<endl;
	*mean_value  = Eraw_CB_m0->getVal();
        *mean_error  = Eraw_CB_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_CB_sigma->getVal();
	*sigma_value_error = Eraw_CB_sigma->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

 
        double entries = hh->GetEntries();
        
        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");	
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
	if(r9sup == 2) *EndCapsR9Chain += "All r9";
        string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);     
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",Eraw_CB_m0->getVal(), Eraw_CB_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma = %f +/- %f}",Eraw_CB_sigma->getVal(), Eraw_CB_sigma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#alpha = %f +/- %f}",Eraw_CB_alpha->getVal(), Eraw_CB_alpha->getError()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{n = %f +/- %f}",Eraw_CB_n->getVal(), Eraw_CB_n->getError()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
	latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
	latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{p-value = %0.6e}",*pValue));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);




	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);	

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
	residuals->Delete();
	pulls->Delete();	
	textResid->Delete();
	lineResid->Delete();
	Data_subset->Delete();
	Data->Delete();
	ntplVars->Delete();
	mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
	Photon_r9->Delete();
        isLooseMMG->Delete();
        isMultipleCandidate->Delete();
	Photon_Et->Delete();
	Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_CB_m0->Delete();
        Eraw_CB_sigma->Delete();
        Eraw_CB_alpha->Delete();
        Eraw_CB_n->Delete();
	Eraw_CrystalBall->Delete();



}


void RooLogNormal(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);

	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);


	RooRealVar * Eraw_LN_m0 = new RooRealVar("mmg_s_LN_m0", "CB #Delta m_{0}", 1.0, -5.0, 5.0, "GeV");
        RooRealVar * Eraw_LN_k = new RooRealVar("mmg_s_LN_sigma", "CB ", 0.45, 0.01, 2.0, "GeV");

        RooLognormal * Eraw_LogNormal = new RooLognormal("Eraw_LogNormal", "mmg_s_LogNormal", *mmg_s, *Eraw_LN_m0, *Eraw_LN_k);


	int fewBins = 1;
        Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        {
                cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
                Erawframe->Clear();
                Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Eraw_LogNormal = new RooLognormal("Eraw_LogNormal", "mmg_s_LogNormal", *mmg_s, *Eraw_LN_m0, *Eraw_LN_k);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
                Eraw_LogNormal->fitTo(*Data_subset, Range(RangeMin, RangeMax));
                Eraw_LogNormal->plotOn(Erawframe,Name("mycurve"));
                Erawframe->Draw();

                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",2, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break;
        }

	//RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = Eraw_LogNormal.asTF( RooArgList(mmg_s) );
        f = Eraw_LogNormal->asTF( RooArgList(*mmg_s) );
        *mean_value  = Eraw_LN_m0->getVal();
        *mean_error  = Eraw_LN_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_LN_k->getVal();
        *sigma_value_error = Eraw_LN_k->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",Eraw_LN_m0->getVal(), Eraw_LN_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{k = %f +/- %f}",Eraw_LN_k->getVal(), Eraw_LN_k->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{p-value = %0.6e}",*pValue));

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

        string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

        residuals->Draw("AP");
        nomFichier = "Chi2Residuals";
      
        residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

        TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3); 
        lineResid->Draw();

        TLatex *textResid = new TLatex();
	textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        residuals->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();

        pulls->Draw("AP");
        nomFichier = "Chi2Pulls";

	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

        lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");        
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        pulls->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();
	
	delete EndCapsR9Chain;
        delete tempLegChain;
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_LN_m0->Delete();
        Eraw_LN_k->Delete();
        Eraw_LogNormal->Delete();

}


void RooGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);

	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);


	RooRealVar * sigmean = new RooRealVar("sigmean","B^{#pm} mass",1.0,0.7,1.2);
        RooRealVar * sigwidth = new RooRealVar("sigwidth","B^{#pm} width",0.5,0.0,1.0) ;

        // --- Build Gaussian PDF ---
        RooGaussian * Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth) ;



	int fewBins = 1;
        Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        {
                cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
                Erawframe->Clear();
                Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
                Eraw_Gaussian->fitTo(*Data_subset, Range(RangeMin, RangeMax));
                Eraw_Gaussian->plotOn(Erawframe,Name("mycurve"));
                Erawframe->Draw();

                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",2, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break;
        }

	//RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = Eraw_Gaussian.asTF( RooArgList(mmg_s) );
        f = Eraw_Gaussian->asTF( RooArgList(*mmg_s) );
        *mean_value  = sigmean->getVal();
        *mean_error  = sigmean->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = sigwidth->getVal();
        *sigma_value_error = sigwidth->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	

        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",sigmean->getVal(), sigmean->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma = %f +/- %f}",sigwidth->getVal(), sigwidth->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{p-value = %0.6e}",*pValue));

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

        string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

        residuals->Draw("AP");
        nomFichier = "Chi2Residuals";
      
        residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

        TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3); 
        lineResid->Draw();

        TLatex *textResid = new TLatex();
	textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        residuals->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();

        pulls->Draw("AP");
        nomFichier = "Chi2Pulls";

	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

        lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");        
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        pulls->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();
	
	delete EndCapsR9Chain;
        delete tempLegChain;
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        sigmean->Delete();
        sigwidth->Delete();
        Eraw_Gaussian->Delete();

}


void RooGamma2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);

	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);



 
	// --- Parameters ---
/*	RooRealVar beta("beta","beta",20.0,0.0,200.0);
	RooRealVar gamma("gamma","gamma",1.0,0.0,5.0) ;
	RooRealVar mu("mu","mu",0.5,0.0,1.0) ; 
*/
	RooRealVar * beta = new RooRealVar("beta","beta",20.0,0.0,200.0);
        RooRealVar * gamma = new RooRealVar("gamma","gamma",1.0,0.0,5.0) ;
        RooRealVar * mu = new RooRealVar("mu","mu",0.5,0.0,1.0) ;


	// --- Build Gaussian PDF ---
	RooGamma * Eraw_Gamma = new RooGamma("Eraw_Gamma","mmg_s_Gamma",*mmg_s,*beta,*gamma,*mu);

	int fewBins = 1;
        Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        {
                cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
                Erawframe->Clear();
                Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Eraw_Gamma = new RooGamma("Eraw_Gamma","mmg_s_Gamma",*mmg_s,*beta,*gamma,*mu);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
                Eraw_Gamma->fitTo(*Data_subset, Range(RangeMin, RangeMax));
                Eraw_Gamma->plotOn(Erawframe,Name("mycurve"));
                Erawframe->Draw();

                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break;
        }



        //RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = Eraw_Gamma.asTF( RooArgList(mmg_s) );
       	f = Eraw_Gamma->asTF( RooArgList(*mmg_s) ); 
	//*mean_value  = sigmean.getVal();
        //*mean_error  = sigmean.getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{#beta = %f +/- %f}",beta->getVal(), beta->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#gamma = %f +/- %f}",gamma->getVal(), gamma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#mu = %f +/- %f}",mu->getVal(), mu->getError()));
	latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{p-value = %0.6e}",*pValue));

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

        string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

        residuals->Draw("AP");
        nomFichier = "Chi2Residuals";
      
        residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

        TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3); 
        lineResid->Draw();

        TLatex *textResid = new TLatex();
	textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        residuals->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();

        pulls->Draw("AP");
        nomFichier = "Chi2Pulls";

	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

        lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");        
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        pulls->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();
	
	delete EndCapsR9Chain;
        delete tempLegChain;
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        beta->Delete();
        gamma->Delete();
        mu->Delete();
        Eraw_Gamma->Delete();

}

void RooBifurcatedGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
	RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
	RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, ""); 
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, ""); 
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG, *isMultipleCandidate, *weight_pileUp);


        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *isLooseMMG, *isMultipleCandidate);
        //RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isLooseMMG, Photon_Et);
        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars);
        RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
        //RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        //Erawframe = mmg_s->frame();
        //Data_subset->plotOn(Erawframe);
        //Data_subset.plotOn(Erawframe);

        //Erawframe->Draw();

 
	// --- Parameters ---
	/*RooRealVar BifurGauss_mean("BifurGauss_mean","BifurGauss_mean",1.0,0.7,1.2);
	RooRealVar BifurGauss_sigmaR("BifurGauss_sigmaR","BifurGauss_sigmaR",0.5,0.0,1.0);
	RooRealVar BifurGauss_sigmaL("BifurGauss_sigmaL","BifurGauss_sigmaL",0.5,0.0,1.0);
 */

	RooRealVar * BifurGauss_mean = new RooRealVar("BifurGauss_mean","BifurGauss_mean",0.0,-0.3,0.2);
        RooRealVar * BifurGauss_sigmaR = new RooRealVar("BifurGauss_sigmaR","BifurGauss_sigmaR",0.5,0.0,1.0);
        RooRealVar * BifurGauss_sigmaL = new RooRealVar("BifurGauss_sigmaL","BifurGauss_sigmaL",0.5,0.0,1.0);


	// --- Build Bifurcated Gaussian PDF ---
	RooBifurGauss * Eraw_BifurGauss = new RooBifurGauss("Eraw_BifurGauss","mmg_s_BifurGauss",*mmg_s,*BifurGauss_mean,*BifurGauss_sigmaR,*BifurGauss_sigmaL);


	int fewBins = 1;
        Double_t Chi2J;

	RooFitResult *res;
	Double_t minNll;
		

        //for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        int RightBinning = 60;
	for(int i = 0; i<4; i++)
	{
                Erawframe->Clear();
                Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Eraw_BifurGauss = new RooBifurGauss("Eraw_BifurGauss","mmg_s_BifurGauss",*mmg_s,*BifurGauss_mean,*BifurGauss_sigmaR,*BifurGauss_sigmaL);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
                //Eraw_BifurGauss->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true));
                res = Eraw_BifurGauss->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true),Save());
		res->Print();
		minNll = res->minNll();
		Eraw_BifurGauss->plotOn(Erawframe,Name("mycurve"));
                Erawframe->Draw();


                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break; 


        }



        //Eraw_CrystalBall.fitTo(*Data_subset, Range(0.6,1.2));
        //Eraw_BifurGauss.fitTo(Data_subset, Range(0.8,1.2)); //ATTENTION NE PAS SUPPRIMER!!!!
       
	//double maxDistri = hh->GetMaximumBin() * 0.01;
        //Eraw_BifurGauss.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
        //Eraw_BifurGauss.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms)); 
        //Eraw_BifurGauss.fitTo(*Data_subset, Range(RangeMin, RangeMax));
	//Eraw_BifurGauss.plotOn(Erawframe);
        //Erawframe->Draw();



        //RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = Eraw_BifurGauss.asTF( RooArgList(mmg_s) );
        f = Eraw_BifurGauss->asTF( RooArgList(*mmg_s) );
	*mean_value  = BifurGauss_mean->getVal();
        *mean_error  = BifurGauss_mean->getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	*minLogLikelihood = minNll;


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	//latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean = %f +/- %f}",BifurGauss_mean->getVal(), BifurGauss_mean->getError()));
        //latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{R} = %f +/- %f}",BifurGauss_sigmaR->getVal(), BifurGauss_sigmaR->getError()));
	//latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#sigma_{L} = %f +/- %f}",BifurGauss_sigmaL->getVal(), BifurGauss_sigmaL->getError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean = %f^{+ %f}_{- %f}}",BifurGauss_mean->getVal(), BifurGauss_mean->getErrorHi(),BifurGauss_mean->getErrorLo()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{R} = %f^{+ %f}_{- %f}}",BifurGauss_sigmaR->getVal(), BifurGauss_sigmaR->getErrorHi(),BifurGauss_sigmaR->getErrorLo()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#sigma_{L} = %f^{+ %f}_{- %f}}",BifurGauss_sigmaL->getVal(), BifurGauss_sigmaL->getErrorHi(),BifurGauss_sigmaL->getErrorLo()));	
	latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{p-value = %0.6e}",*pValue));


	cout<<endl<<"HERE"<<endl<<"BifurGauss_mean->getErrorHi() = "<<BifurGauss_mean->getErrorHi()<<endl<<"BifurGauss_mean->getErrorLo() = "<<BifurGauss_mean->getErrorLo()<<endl<<endl;


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// RooMCStudy //////////

	
	RooMCStudy* mcs = new RooMCStudy(*Eraw_BifurGauss,*mmg_s,FitModel(*Eraw_BifurGauss),Silence(),Binned()) ;
        RooChi2MCSModule chi2mod ;
        mcs->addModule(chi2mod) ;

        mcs->generateAndFit(2000,entriesInt) ;

        //RooPlot* frame = mcs->plotNLL(Bins(50));
        RooPlot* frame = mcs->plotNLL(Name("myhistTest"));
	frame->Draw();

	//RooRealVar * VArTest = (RooRealVar*) frame->getPlotVar(); //artefact

	
	RooHist* histTest = frame->getHist("myhistTest");
	cout<<endl<<"histTest->GetMean() = "<<histTest->GetMean()<<endl;

	*differenceLogLikelihood = fabs(histTest->GetMean() - minNll);

	TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());

	textResid->DrawLatex(0.61, 0.90, Form("#color[4]{Fit_minll = %f}",minNll));
	textResid->DrawLatex(0.61, 0.85, Form("#color[4]{Distrib_mean = %f}",histTest->GetMean()));
	textResid->DrawLatex(0.61, 0.80, Form("#color[3]{Difference = %f}",*differenceLogLikelihood));
	
	//RooAbsRealLValue VArTest = frame->getPlotVar();


	enregistrementPlots(nomDossier, "LogLikelihood", EndCaps, j, c2);

	c2->Clear();


	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";


        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);
	

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);


	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();


        //TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);	

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
	res->Delete();
	mcs->Delete();
	frame->Delete();
	residuals->Delete();
	pulls->Delete();	
	textResid->Delete();
	lineResid->Delete();
	Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
	Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
	isMultipleCandidate->Delete();
        Photon_Et->Delete();
	Photon_SC_rawEt->Delete();
	Photon_E->Delete();
	Photon_SC_brem->Delete();
	weight_pileUp->Delete();
	BifurGauss_mean->Delete();
	BifurGauss_sigmaR->Delete();
	BifurGauss_sigmaL->Delete();
	Eraw_BifurGauss->Delete();



}


void RooSumGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, ""); 
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, ""); 
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, ""); 
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, ""); 
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, ""); 
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");     
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);
	
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);

 
	// --- Parameters ---
	//RooRealVar sigmean1("sigmean1","sigmean1",1.0,0.6,1.4);
	//RooRealVar sigwidth1("sigwidth1","sigwidth1",0.5,0.0,1.0) ;
 
	RooRealVar * sigmean1 = new RooRealVar("sigmean1","sigmean1",1.0,0.6,1.4);
        RooRealVar * sigwidth1 = new RooRealVar("sigwidth1","sigwidth1",0.5,0.0,1.0) ;

	// --- Build Gaussian PDF ---
	RooGaussian Eraw_Gaussian1("Eraw_Gaussian1","mmg_s_Gaussian1",*mmg_s,*sigmean1,*sigwidth1) ;
	
	// --- Parameters ---
        //RooRealVar sigmean2("sigmean2","sigmean2",0.9,0.6,1.4);
	//RooRealVar sigwidth2("sigwidth2","sigwidth2",0.5,0.0,1.0) ;
        RooRealVar * sigmean2 = new RooRealVar("sigmean2","sigmean2",0.9,0.6,1.4);
	RooRealVar * sigwidth2 = new RooRealVar("sigwidth2","sigwidth2",0.5,0.0,1.0) ;


	// --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian2("Eraw_Gaussian2","mmg_s_Gaussian2",*mmg_s,*sigmean2,*sigwidth2) ;


	// --- Parameters ---
        //RooRealVar sigmean3("sigmean3","sigmean3",1.1,0.6,1.4);
	//RooRealVar sigwidth3("sigwidth3","sigwidth3",0.5,0.0,1.0) ;
	RooRealVar * sigmean3 = new RooRealVar("sigmean3","sigmean3",1.1,0.6,1.4);
	RooRealVar * sigwidth3 = new RooRealVar("sigwidth3","sigwidth3",0.5,0.0,1.0) ;


        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian3("Eraw_Gaussian3","mmg_s_Gaussian3",*mmg_s,*sigmean3,*sigwidth3) ;


	// --- Parameters ---
        //RooRealVar sigmean4("sigmean4","sigmean4",0.8,0.6,1.4);
        //RooRealVar sigwidth4("sigwidth4","sigwidth4",0.5,0.0,1.0);
	RooRealVar * sigmean4 = new RooRealVar("sigmean4","sigmean4",0.8,0.6,1.4);
	RooRealVar * sigwidth4 = new RooRealVar("sigwidth4","sigwidth4",0.5,0.0,1.0);


        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian4("Eraw_Gaussian4","mmg_s_Gaussian4",*mmg_s,*sigmean4,*sigwidth4) ;	


	// --- Parameters ---
        //RooRealVar sigmean5("sigmean5","sigmean5",1.2,0.6,1.4);
        //RooRealVar sigwidth5("sigwidth5","sigwidth5",0.5,0.0,1.0);
	RooRealVar * sigmean5 = new RooRealVar("sigmean5","sigmean5",1.2,0.6,1.4);
	RooRealVar * sigwidth5 = new RooRealVar("sigwidth5","sigwidth5",0.5,0.0,1.0);

        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian5("Eraw_Gaussian5","mmg_s_Gaussian5",*mmg_s,*sigmean5,*sigwidth5) ;


	// --- Parameters ---
        //RooRealVar sigmean6("sigmean6","sigmean6",0.7,0.6,1.4);
        //RooRealVar sigwidth6("sigwidth6","sigwidth6",0.5,0.0,1.0);
	RooRealVar * sigmean6 = new RooRealVar("sigmean6","sigmean6",0.7,0.6,1.4);
	RooRealVar * sigwidth6 = new RooRealVar("sigwidth6","sigwidth6",0.5,0.0,1.0);

        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian6("Eraw_Gaussian6","mmg_s_Gaussian6",*mmg_s,*sigmean6,*sigwidth6) ;


	// --- Parameters ---
        //RooRealVar sigmean7("sigmean7","sigmean7",1.3,0.6,1.4);
        //RooRealVar sigwidth7("sigwidth7","sigwidth7",0.5,0.0,1.0);
	RooRealVar * sigmean7 = new RooRealVar("sigmean7","sigmean7",1.3,0.6,1.4);
	RooRealVar * sigwidth7 = new RooRealVar("sigwidth7","sigwidth7",0.5,0.0,1.0);

        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian7("Eraw_Gaussian7","mmg_s_Gaussian7",*mmg_s,*sigmean7,*sigwidth7) ;




	// --- Build Sum of Gaussian PDF ---
/*
	RooRealVar * g1Frac = new RooRealVar("g1Frac","g1Frac",0.1428);
	RooRealVar * g2Frac = new RooRealVar("g2Frac","g2Frac",0.1428);
	RooRealVar * g3Frac = new RooRealVar("g3Frac","g3Frac",0.1428);
	RooRealVar * g4Frac = new RooRealVar("g4Frac","g4Frac",0.1428);
	RooRealVar * g5Frac = new RooRealVar("g5Frac","g5Frac",0.1428);
        RooRealVar * g6Frac = new RooRealVar("g6Frac","g6Frac",0.1428);
        RooRealVar * g7Frac = new RooRealVar("g7Frac","g7Frac",0.1428);	
*/
	double princ = 0.99;

	double frac = (1.0 - princ) / 6.0;

	RooRealVar * g1Frac = new RooRealVar("g1Frac","g1Frac",princ);
        RooRealVar * g2Frac = new RooRealVar("g2Frac","g2Frac",frac);
        RooRealVar * g3Frac = new RooRealVar("g3Frac","g3Frac",frac);
        RooRealVar * g4Frac = new RooRealVar("g4Frac","g4Frac",frac);
        RooRealVar * g5Frac = new RooRealVar("g5Frac","g5Frac",frac);
        RooRealVar * g6Frac = new RooRealVar("g6Frac","g6Frac",frac);
        RooRealVar * g7Frac = new RooRealVar("g7Frac","g7Frac",frac);	



	RooAddPdf SumGaussians("SumGaussians","SumGaussians", RooArgList(Eraw_Gaussian1,Eraw_Gaussian2,Eraw_Gaussian3,Eraw_Gaussian4,Eraw_Gaussian5,Eraw_Gaussian6,Eraw_Gaussian7), RooArgList(*g1Frac,*g2Frac,*g3Frac,*g4Frac,*g5Frac,*g6Frac,*g7Frac));


        //Eraw_CrystalBall.fitTo(*Data_subset, Range(0.6,1.2));
        //Eraw_CrystalBall.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
        double maxDistri = hh->GetMaximumBin() * 0.01;
        //SumGaussians.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
        SumGaussians.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
        //SumGaussians.fitTo(Data_subset, Range(0.4,1.6));
	SumGaussians.plotOn(Erawframe);
        Erawframe->Draw();



        //RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = SumGaussians.asTF( RooArgList(mmg_s) );
	f = SumGaussians.asTF( RooArgList(*mmg_s) );        

	//*mean_value  = sigmean.getVal();
        //*mean_error  = sigmean.getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


	cout<<"SigmaR(f, 0.0, 2.0) = "<<SigmaR(f, 0.0, 2.0)<<endl;
	cout<<"SigmaL(f, 0.0, 2.0) = "<<SigmaL(f, 0.0, 2.0)<<endl;

	double fxmax = f->GetMaximumX(0.4,1.6,1.E-10,100,false);

        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	TLatex latexLabel;
        latexLabel.SetTextSize(0.028);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");
        latexLabel.DrawLatex(0.16, 0.80, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.75, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.70, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f}",fxmax));
        //latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{ = %f +/- %f}",sigwidth.getVal(), sigwidth.getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#chi^{2} = %f}",Erawframe->chiSquare()));


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
	isMultipleCandidate->Delete();
        Photon_Et->Delete();
	sigmean1->Delete();
        sigwidth1->Delete();
	sigmean2->Delete();
        sigwidth2->Delete();
	sigmean3->Delete();
        sigwidth3->Delete();
	sigmean4->Delete();
        sigwidth4->Delete();
	sigmean5->Delete();
        sigwidth5->Delete();
	sigmean6->Delete();
        sigwidth6->Delete();
	sigmean7->Delete();
	sigwidth7->Delete();
	g1Frac->Delete();
	g2Frac->Delete();
	g3Frac->Delete();
	g4Frac->Delete();
	g5Frac->Delete();
	g6Frac->Delete();
	g7Frac->Delete();



}




void RooGenericPDF(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
	RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");


        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_r9, *Photon_Et, *isLooseMMG, *isMultipleCandidate);
        //RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isLooseMMG, Photon_Et);
        RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars);
        //RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
        //RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        Erawframe = mmg_s->frame();
        Data_subset->plotOn(Erawframe);
        //Data_subset.plotOn(Erawframe);

        Erawframe->Draw();



 
	// --- Parameters 1/Gaussian ---
	//RooRealVar GenericPDF_mean("GenericPDF_mean","GenericPDF_mean",1.0,0.6,1.4);
	//RooRealVar GenericPDF_mean("GenericPDF_mean","GenericPDF_mean", mean, mean - rms, mean + rms, "GeV");
	//RooRealVar GenericPDF_width("GenericPDF_width","GenericPDF_width",0.5,0.0,1.0);

	RooRealVar * GenericPDF_mean = new RooRealVar("GenericPDF_mean","GenericPDF_mean", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * GenericPDF_width = new RooRealVar("GenericPDF_width","GenericPDF_width",0.5,0.0,1.0);

 
	// --- Build 1/Gaussian PDF ---
	RooGenericPdf Eraw_GenericPdf("Eraw_GenericPdf","mmg_s_GenericPdf","(GenericPDF_width * sqrt(2*3.14159265358979323846))*exp((mmg_s - GenericPDF_mean)*(mmg_s - GenericPDF_mean)/(2*GenericPDF_width*GenericPDF_width))",RooArgSet(*mmg_s, *GenericPDF_mean,*GenericPDF_width));



	// --- Parameters CB ---
/*	RooRealVar Eraw_CB_m0("Eraw_CB_m0", "CB #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar Eraw_CB_sigma("Eraw_CB_sigma", "CB ", rms, 0.0, 2 * rms, "GeV");
        RooRealVar Eraw_CB_alpha("Eraw_CB_alpha", "CB #alpha", 0.9, -5.0, 5.0);
        RooRealVar Eraw_CB_n("Eraw_CB_n", "CB n", 10.0, 0.0, 500.0);    
*/
	RooRealVar * Eraw_CB_m0 = new RooRealVar("Eraw_CB_m0", "CB #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_CB_sigma = new RooRealVar("Eraw_CB_sigma", "CB ", rms, 0.0, 2 * rms, "GeV");
        RooRealVar * Eraw_CB_alpha = new RooRealVar("Eraw_CB_alpha", "CB #alpha", 0.9, -5.0, 5.0);
        RooRealVar * Eraw_CB_n = new RooRealVar("Eraw_CB_n", "CB n", 10.0, 0.0, 500.0); 

/*
        //Ne pas supprimer !!!!//
        RooRealVar Eraw_CB_m0("Eraw_CB_m0", "CB #Delta m_{0}", 1.0, 0.9, 1.05, "GeV");
        RooRealVar Eraw_CB_sigma("Eraw_CB_sigma", "CB ", 0.0215, 0.0, 0.5, "GeV");
        RooRealVar Eraw_CB_alpha("Eraw_CB_alpha", "CB #alpha", 0.9, 0.0, 3.0);
        RooRealVar Eraw_CB_n("Eraw_CB_n", "CB n", 10.0, 0.0, 60.0);
*/

        //RooRealVar Eraw_CB_m0("mmg_s_CB_m0", "CB #Delta m_{0}", 1.0, -5.0, 5.0, "GeV");
        //RooRealVar Eraw_CB_sigma("mmg_s_CB_sigma", "CB ", 0.45, 0.01, 0.5, "GeV");
        //RooRealVar Eraw_CB_alpha("mmg_s_CB_alpha", "CB #alpha", -1.0, -10.01, 10.0);
        //RooRealVar Eraw_CB_n("mmg_s_CB_n", "CB n", 2.0, 0.5, 500.0);


	// --- Build CB PDF ---
        RooCBShape Eraw_CrystalBall("Eraw_CrystalBall","mmg_s_CrystalBall", *mmg_s, *Eraw_CB_m0, *Eraw_CB_sigma, *Eraw_CB_alpha, *Eraw_CB_n);

	RooProdPdf prod("CB/Gaussian","CB/Gaussian",RooArgList(Eraw_CrystalBall,Eraw_GenericPdf));



        //Eraw_CrystalBall.fitTo(*Data_subset, Range(0.6,1.2));
        //prod.fitTo(Data_subset, Range(0.6,1.2)); //ATTENTION NE PAS SUPPRIMER!!!!
	double maxDistri = hh->GetMaximumBin() * 0.01;
        //prod.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
        prod.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
        prod.plotOn(Erawframe);
        Erawframe->Draw();



        //RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = prod.asTF( RooArgList(mmg_s) );
        f = prod.asTF( RooArgList(*mmg_s) );
	*mean_value  = Eraw_CB_m0->getVal();
        *mean_error  = Eraw_CB_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_CB_sigma->getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	TLatex latexLabel;
        latexLabel.SetTextSize(0.028);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");
        latexLabel.DrawLatex(0.16, 0.80, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        
	latexLabel.DrawLatex(0.16, 0.75, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.70, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        latexLabel.SetTextSize(0.023);
        latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean_{1/Gauss} = %f +/- %f}",GenericPDF_mean->getVal(), GenericPDF_mean->getError()));
	latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{width_{1/Gauss} = %f +/- %f}",GenericPDF_width->getVal(), GenericPDF_width->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{m_{0}_{CB} = %f +/- %f}",Eraw_CB_m0->getVal(), Eraw_CB_m0->getError()));
	latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{#sigma_{CB} = %f +/- %f}",Eraw_CB_sigma->getVal(), Eraw_CB_sigma->getError()));
	latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{#alpha_{CB} = %f +/- %f}",Eraw_CB_alpha->getVal(), Eraw_CB_alpha->getError()));
	latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{n_{CB} = %f +/- %f}",Eraw_CB_n->getVal(), Eraw_CB_n->getError()));
        latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{#chi^{2} = %f}",Erawframe->chiSquare()));


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
	isMultipleCandidate->Delete();
        Photon_Et->Delete();
	GenericPDF_mean->Delete();
	GenericPDF_width->Delete();
	Eraw_CB_m0->Delete();
        Eraw_CB_sigma->Delete();
        Eraw_CB_alpha->Delete();
        Eraw_CB_n->Delete();


}

void RooLandau2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
	RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");	
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, ""); 
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, ""); 
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG, *isMultipleCandidate, *weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *isLooseMMG, *isMultipleCandidate, *weight_pileUp);
	//RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isLooseMMG, Photon_Et);
        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars); //non weighté
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
	//RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        
/*	
	Erawframe = mmg_s->frame();
	Data_subset->plotOn(Erawframe,Name("myhist"),Binning(300));
        Erawframe->Draw();
*/

	RooRealVar * Eraw_Landau_m0 = new RooRealVar("Eraw_Landau_m0", "Landau #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_Landau_sigma = new RooRealVar("Eraw_Landau_sigma", "Landau ", rms, 0.0, 2 * rms, "GeV");


/*
	//Ne pas supprimer !!!!//
	RooRealVar Eraw_Landau_m0("Eraw_Landau_m0", "Landau #Delta m_{0}", 1.0, 0.9, 1.05, "GeV");
        RooRealVar Eraw_Landau_sigma("Eraw_Landau_sigma", "Landau ", 0.0215, 0.0, 0.5, "GeV");
        RooRealVar Eraw_Landau_alpha("Eraw_Landau_alpha", "Landau #alpha", 0.9, 0.0, 3.0);
        RooRealVar Eraw_Landau_n("Eraw_Landau_n", "Landau n", 10.0, 0roo.0, 60.0);
*/
        RooLandau * Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma);

	int fewBins = 1;
	Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	{
		cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
		Erawframe->Clear();
		Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);		
		Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma); 

		Erawframe = mmg_s->frame();
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
		//Erawframe->Draw();
        	//Eraw_Landau.fitTo(*Data_subset, Range(0.6,1.2));
		//Eraw_Landau.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
		//Eraw_Landau.fitTo(Data_subset, Range(mean - 5 * rms,mean + 5 * rms));
		//double maxDistri = hh->GetMaximumBin() * 0.01;
		//Eraw_Landau.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
		//Eraw_Landau.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
		Eraw_Landau->fitTo(*Data_subset, Range(RangeMin, RangeMax));
		//RooArgSet* Eraw_Landau_param = Eraw_Landau.getVariables();
		//Eraw_Landau_param->Print("v");
		//Eraw_Landau.plotOn(Erawframe);
        	Eraw_Landau->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


		Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",2, JanChi2, DegreesOfFreedom, pValue, &fewBins);
		if(fewBins == 0) break; 


	}


	//RooDataHist data("data", "dataset with mmg_s",*mmg_s,hh);
	//TF1 * f = Eraw_Landau.asTF( RooArgList(mmg_s) );
	f = Eraw_Landau->asTF( RooArgList(*mmg_s) );
	//ftest = Eraw_Landau.asTF( RooArgList(*mmg_s) );
	//cout<<endl<<"ftest->GetXmax() = "<<ftest->GetXmax()<<endl;
	//cout<<endl<<"ftest->GetMaximum = "<<ftest->GetMaximum(0.0, 2.0,1.E-10,100,false)<<endl;
	*mean_value  = Eraw_Landau_m0->getVal();
        *mean_error  = Eraw_Landau_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_Landau_sigma->getVal();
	*sigma_value_error = Eraw_Landau_sigma->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

 
        double entries = hh->GetEntries();
        
        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);     
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",Eraw_Landau_m0->getVal(), Eraw_Landau_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{landau} = %f +/- %f}",Eraw_Landau_sigma->getVal(), Eraw_Landau_sigma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
	latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
	latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{p-value = %0.6e}",*pValue));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);




	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);	

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
	residuals->Delete();
	pulls->Delete();	
	textResid->Delete();
	lineResid->Delete();
	Data_subset->Delete();
	Data->Delete();
	ntplVars->Delete();
	mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
	Photon_r9->Delete();
        isLooseMMG->Delete();
        isMultipleCandidate->Delete();
	Photon_Et->Delete();
	Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_Landau_m0->Delete();
        Eraw_Landau_sigma->Delete();
	Eraw_Landau->Delete();



}

void RooLandauConvGaussian(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
	RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");	
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, ""); 
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, ""); 
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG, *isMultipleCandidate, *weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *isLooseMMG, *isMultipleCandidate, *weight_pileUp);
	//RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isLooseMMG, Photon_Et);
        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars); //non weighté
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
	//RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        
/*	
	Erawframe = mmg_s->frame();
	Data_subset->plotOn(Erawframe,Name("myhist"),Binning(300));
        Erawframe->Draw();
*/

	RooRealVar * Eraw_Landau_m0 = new RooRealVar("Eraw_Landau_m0", "Landau #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_Landau_sigma = new RooRealVar("Eraw_Landau_sigma", "Landau ", rms , 0.0,1.0, "GeV");

	// --- Build Landau PDF ---
        RooLandau * Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma);


	/*
        RooRealVar * sigmean = new RooRealVar("sigmean","B^{#pm} mass",0.0,0.0,0.0);
        RooRealVar * sigwidth = new RooRealVar("sigwidth","B^{#pm} width",rms, 0.0,1.0);


        // --- Build Gaussian PDF ---
        RooGaussian * Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth) ;


        // CONVOLUTION
        RooFFTConvPdf * model = new RooFFTConvPdf("model", "model", *mmg_s, *Eraw_Landau, *Eraw_Gaussian);
*/

        int fewBins = 1; 
        Double_t Chi2J;

        Erawframe = mmg_s->frame();
        Data_subset->plotOn(Erawframe,Name("myhist"));
        Eraw_Landau->fitTo(*Data_subset, Range(RangeMin, RangeMax));
        Eraw_Landau->plotOn(Erawframe,Name("mycurve"));
        Erawframe->Draw();

        RooRealVar * sigmean = new RooRealVar("sigmean","B^{#pm} mass",0.0,0.0,0.0);
        RooRealVar * sigwidth = new RooRealVar("sigwidth","B^{#pm} width",(Eraw_Landau_sigma->getVal() / 4.0), 0.0,(Eraw_Landau_sigma->getVal() / 1.5));

        // --- Build Gaussian PDF ---
        RooGaussian * Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth) ;


        // CONVOLUTION
        RooFFTConvPdf * model = new RooFFTConvPdf("model", "model", *mmg_s, *Eraw_Landau, *Eraw_Gaussian);




	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	{
		cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
		Erawframe->Clear();
		Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);		
		Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma);
		Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth);
		model = new RooFFTConvPdf("model", "model", *mmg_s, *Eraw_Landau, *Eraw_Gaussian); 

		Erawframe = mmg_s->frame();
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
		//Erawframe->Draw();
        	//model.fitTo(*Data_subset, Range(0.6,1.2));
		//model.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
		//model.fitTo(Data_subset, Range(mean - 5 * rms,mean + 5 * rms));
		//double maxDistri = hh->GetMaximumBin() * 0.01;
		//model.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
		//model.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
		model->fitTo(*Data_subset, Range(RangeMin, RangeMax));
		//RooArgSet* model_param = model.getVariables();
		//model_param->Print("v");
		//model.plotOn(Erawframe);
        	model->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


		Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",4, JanChi2, DegreesOfFreedom, pValue, &fewBins);
		if(fewBins == 0) break; 


	}

	// Create function that represents derivative
        //RooAbsReal * deriv = (RooAbsReal*) model->derivative(*mmg_s,1) ;
        RooDerivative * deriv = model->derivative(*mmg_s,1) ;

        // Find point where derivative is zero
        Double_t x_max = deriv->findRoot(*mmg_s,RangeMin, RangeMax,0) ;


	//RooDataHist data("data", "dataset with mmg_s",*mmg_s,hh);
	//TF1 * f = model.asTF( RooArgList(mmg_s) );
	f = model->asTF( RooArgList(*mmg_s) );
	//ftest = model.asTF( RooArgList(*mmg_s) );
	//cout<<endl<<"ftest->GetXmax() = "<<ftest->GetXmax()<<endl;
	//cout<<endl<<"ftest->GetMaximum = "<<ftest->GetMaximum(0.0, 2.0,1.E-10,100,false)<<endl;
	*mean_value  = x_max;
        *mean_error  = Eraw_Landau_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_Landau_sigma->getVal();
	*sigma_value_error = Eraw_Landau_sigma->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

 
        double entries = hh->GetEntries();
        
        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);     
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0 Landau} = %f +/- %f}",Eraw_Landau_m0->getVal(), Eraw_Landau_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{Landau} = %f +/- %f}",Eraw_Landau_sigma->getVal(), Eraw_Landau_sigma->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{m_{0 Gaus} = %f +/- %f}",sigmean->getVal(), sigmean->getError()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{#sigma_{Gaus} = %f +/- %f}",sigwidth->getVal(), sigwidth->getError()));
	latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{m_{0 Conv} = %f}",*mean_value));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.55, Form("#color[4]{p-value = %0.6e}",*pValue));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);




	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);	

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
	residuals->Delete();
	pulls->Delete();	
	textResid->Delete();
	lineResid->Delete();
	Data_subset->Delete();
	Data->Delete();
	ntplVars->Delete();
	mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
	Photon_r9->Delete();
        isLooseMMG->Delete();
        isMultipleCandidate->Delete();
	Photon_Et->Delete();
	Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_Landau_m0->Delete();
        Eraw_Landau_sigma->Delete();
	Eraw_Landau->Delete();
	sigmean->Delete();
	sigwidth->Delete();
	Eraw_Gaussian->Delete();
	model->Delete();



}


void RooKernel(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{
	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s_{RECO}", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
	RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);



        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);


        RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);

	RooKeysPdf * K = new RooKeysPdf("K","K",*mmg_s,*Data_subset);
	Erawframe = mmg_s->frame();
	
	//Data_subset->plotOn(Erawframe,Name("myhist"),Binning(60));
	Data_subset->plotOn(Erawframe,Name("myhist"));
	RooFitResult* r = K->fitTo(*Data_subset,Minos(true),Save());
        K->plotOn(Erawframe,Name("mycurve"));
        Erawframe->Draw();


	double entries = hh->GetEntries();

        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
        string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        /*latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0 Landau} = %f +/- %f}",Eraw_Landau_m0->getVal(), Eraw_Landau_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{Landau} = %f +/- %f}",Eraw_Landau_sigma->getVal(), Eraw_Landau_sigma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{m_{0 Gaus} = %f +/- %f}",sigmean->getVal(), sigmean->getError()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{#sigma_{Gaus} = %f +/- %f}",sigwidth->getVal(), sigwidth->getError()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{m_{0 Conv} = %f}",*mean_value));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.55, Form("#color[4]{p-value = %0.6e}",*pValue));
	*/

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
        //residuals->Delete();
        //pulls->Delete();     
        //textResid->Delete();
        //lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
	isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        K->Delete();

}

void RooVoigtian2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{
	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 100);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);



        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);


	//RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "");
	RooDataSet *Data_subset = new RooDataSet("Data_subset", "Data_subset", Tree_Data, *ntplVars, temp, "weight_pileUp");


	RooRealVar * meanV = new RooRealVar("meanV","meanV",0.0,-0.1,0.1);
        RooRealVar * sigma = new RooRealVar("sigma","sigma",0.5,0.0,1.0);
        RooRealVar * width = new RooRealVar("width","width",0.5,0.0,1.0);


	// --- Build Bifurcated Gaussian PDF ---
	RooVoigtian * Voigtian = new RooVoigtian("Voigtian","Voigtian",*mmg_s,*meanV,*sigma,*width);


	int fewBins = 1;
        Double_t Chi2J;

	RooFitResult *res;
	Double_t minNll;
		

        //for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        int RightBinning = 25;
	//for(int i = 0; i<4; i++)
	for(int i = 0; i<1; i++)
	{
                Erawframe->Clear();
                //Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Voigtian = new RooVoigtian("Voigtian","Voigtian",*mmg_s,*meanV,*sigma,*width);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning),DataError(RooAbsData::SumW2));
                //Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true));
                //res = Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true),Save(),SumW2Error(kFALSE));
		res = Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Save(),SumW2Error(kTRUE));
		res->Print();
		minNll = res->minNll();
		Voigtian->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break; 


        }
	cout << endl << "ICI !!!!!!!!!!!!!!!!!! "<<Data_subset->weight() << endl ;

        f = Voigtian->asTF( RooArgList(*mmg_s) );
	*mean_value  = meanV->getVal();
        *mean_error  = meanV->getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	*minLogLikelihood = minNll;


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
	latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011, #sqrt{s} = 7 TeV");
        if(isMC == 1) latexLabel.DrawLatex(0.17, 0.88, "Simulation");
        if(isMC == 0) latexLabel.DrawLatex(0.17, 0.88, "Data, #int L = 4,89 fb^{-1}"); //lumi a changer !!!

	if(EndCaps == 0) latexLabel.DrawLatex(0.17, 0.83,"ECAL Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.17, 0.83,"ECAL Endcaps");
        if(r9sup == 0 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
        if(r9sup == 0 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
        if(r9sup == 1 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,94");
        if(r9sup == 1 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,95");
        if(r9sup == 2) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");
	
	latexLabel.DrawLatex(0.17, 0.73, Form("Entries = %d",entriesInt));


	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean = %f^{+ %f}_{ %f}}",meanV->getVal(), meanV->getErrorHi(),meanV->getErrorLo()));
	latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma = %f^{+ %f}_{ %f}}",sigma->getVal(), sigma->getErrorHi(),sigma->getErrorLo()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{width = %f^{+ %f}_{ %f}}",width->getVal(), width->getErrorHi(),width->getErrorLo()));	


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// RooMCStudy //////////

	
	RooMCStudy* mcs = new RooMCStudy(*Voigtian,*mmg_s,FitModel(*Voigtian),Silence(),Binned()) ;
        RooChi2MCSModule chi2mod ;
        mcs->addModule(chi2mod) ;

        mcs->generateAndFit(2000,entriesInt) ;

        //RooPlot* frame = mcs->plotNLL(Bins(50));
        RooPlot* frame = mcs->plotNLL(Name("myhistTest"));
	frame->Draw();

	//RooRealVar * VArTest = (RooRealVar*) frame->getPlotVar(); //artefact

	
	RooHist* histTest = frame->getHist("myhistTest");
	cout<<endl<<"histTest->GetMean() = "<<histTest->GetMean()<<endl;

	*differenceLogLikelihood = fabs(histTest->GetMean() - minNll);

	TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());

	textResid->DrawLatex(0.61, 0.90, Form("#color[4]{Fit_minll = %f}",minNll));
	textResid->DrawLatex(0.61, 0.85, Form("#color[4]{Distrib_mean = %f}",histTest->GetMean()));
	textResid->DrawLatex(0.61, 0.80, Form("#color[3]{Difference = %f}",*differenceLogLikelihood));
	
	//RooAbsRealLValue VArTest = frame->getPlotVar();


	enregistrementPlots(nomDossier, "LogLikelihood", EndCaps, j, c2);

	c2->Clear();


	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";


        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);
	

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);


	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s_{RECO}");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();


        //TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s_{RECO}");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	



	c2->Clear();

        //delete EndCapsR9Chain;
        //delete tempLegChain;
        res->Delete();
        mcs->Delete();
        frame->Delete();
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        //Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        meanV->Delete();
        sigma->Delete();
        width->Delete();
        Voigtian->Delete();

}


void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1)
{    
     
	if(EndCaps == 0) nomDossier += "_EB/";
	if(EndCaps == 1) nomDossier += "_EE/";
	//mkdir(nomDossier.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	system(Form("mkdir -p %s", nomDossier.c_str()));
	if(iteration != 10000) nomFichier += Form("%d",iteration);
        c1->Print(Form("%s%s.root",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.C",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.pdf",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.ps",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.png",nomDossier.c_str(),nomFichier.c_str()));
        return;
}

void RangeEstimator(double pourcentage, double centralValue, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax)
{

	float mmg_s;
        float Photon_SC_Eta;
        float Photon_r9;
        float isLooseMMG;
        float isMultipleCandidate;
        float Photon_Et;
        chain->SetBranchAddress("mmg_s",&mmg_s);
        chain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
        chain->SetBranchAddress("Photon_r9",&Photon_r9);
        chain->SetBranchAddress("isLooseMMG",&isLooseMMG);
        chain->SetBranchAddress("isMultipleCandidate",&isMultipleCandidate);
        chain->SetBranchAddress("Photon_Et",&Photon_Et);
        vector <float> OneOverKrecoClassicalVector;

	for (int ievt = 0 ; ievt < chain->GetEntries() ; ievt++)
        {
                chain->GetEntry(ievt);

                if(Endcaps == 0)
                {

                        if(( (fabs(Photon_SC_Eta)) > 0.018 || ( (fabs(Photon_SC_Eta)) < 0.423 && (fabs(Photon_SC_Eta)) > 0.461 ) || ( (fabs(Photon_SC_Eta)) < 0.770 && (fabs(Photon_SC_Eta)) > 0.806 ) || ( (fabs(Photon_SC_Eta)) < 1.127 && (fabs(Photon_SC_Eta)) > 1.163 ) ) && (fabs(Photon_SC_Eta)) < 1.4442 && isMultipleCandidate == 0 && isLooseMMG != 0 && Photon_r9 > 0.94 && Photon_Et > Varmin && Photon_Et <= Varmax)
                        {

                                OneOverKrecoClassicalVector.push_back(mmg_s);

                        }


                }

                if(Endcaps == 1)
                {

        //cout<<endl<<BremVector[1000]<<endl;
                        if((fabs(Photon_SC_Eta)) > 1.566 && isLooseMMG != 0 && isMultipleCandidate == 0 && Photon_r9 > 0.95 && Photon_Et > Varmin && Photon_Et <= Varmax)  OneOverKrecoClassicalVector.push_back(mmg_s);

                }


        }
        sort(OneOverKrecoClassicalVector.begin(), OneOverKrecoClassicalVector.end());

	cout<<endl<<"OneOverKrecoClassicalVector.size() = "<<OneOverKrecoClassicalVector.size()<<endl;

	int meanVector = 0;
	for(int h = 0; h < OneOverKrecoClassicalVector.size(); h++)
	{
		if(OneOverKrecoClassicalVector[h] >= centralValue)
		{
			meanVector = h;
			h = OneOverKrecoClassicalVector.size();
			break;
		}


	}

        
	cout<<endl<<"meanVector = "<<meanVector<<endl;

        double stopVectorMoins = meanVector - OneOverKrecoClassicalVector.size() * pourcentage / 2.0;
        double stopVectorPlus = meanVector + OneOverKrecoClassicalVector.size() * pourcentage / 2.0;
        cout<<endl<<"stopVectorPlus = "<<stopVectorPlus<<endl;
        int intStopVectorMoins = (int) stopVectorMoins;
        int intStopVectorPlus = (int) stopVectorPlus;
        cout<<endl<<"intStopVectorMoins = "<<intStopVectorMoins<<endl;
        cout<<endl<<"intStopVectorPlus = "<<intStopVectorPlus<<endl;

        cout<<endl<<"OneOverKrecoClassicalVector[meanVector] = "<<OneOverKrecoClassicalVector[meanVector]<<endl;

        cout<<endl<<"Le range est compris entre : "<<OneOverKrecoClassicalVector[intStopVectorMoins]<<" et "<<OneOverKrecoClassicalVector[intStopVectorPlus]<<endl;

	*MinRange = OneOverKrecoClassicalVector[intStopVectorMoins];
	*MaxRange = OneOverKrecoClassicalVector[intStopVectorPlus]; 


}


void RangeEstimator2(double pourcentage, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax)
{

	float mmg_s;
        float Photon_SC_Eta;
        float Photon_r9;
        float isLooseMMG;
        float isMultipleCandidate;
        float Photon_Et;
        chain->SetBranchAddress("mmg_s",&mmg_s);
        chain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
        chain->SetBranchAddress("Photon_r9",&Photon_r9);
        chain->SetBranchAddress("isLooseMMG",&isLooseMMG);
        chain->SetBranchAddress("isMultipleCandidate",&isMultipleCandidate);
        chain->SetBranchAddress("Photon_Et",&Photon_Et);
        vector <float> OneOverKrecoClassicalVector;

        for (int ievt = 0 ; ievt < chain->GetEntries() ; ievt++)
        {
                chain->GetEntry(ievt);

                if(Endcaps == 0)
                {

                        if(( (fabs(Photon_SC_Eta)) > 0.018 || ( (fabs(Photon_SC_Eta)) < 0.423 && (fabs(Photon_SC_Eta)) > 0.461 ) || ( (fabs(Photon_SC_Eta)) < 0.770 && (fabs(Photon_SC_Eta)) > 0.806 ) || ( (fabs(Photon_SC_Eta)) < 1.127 && (fabs(Photon_SC_Eta)) > 1.163 ) ) && (fabs(Photon_SC_Eta)) < 1.4442 && isMultipleCandidate == 0 && isLooseMMG != 0 && Photon_r9 > 0.94 && Photon_Et > Varmin && Photon_Et <= Varmax)
                        {

                                OneOverKrecoClassicalVector.push_back(mmg_s);

                        }


                }

                if(Endcaps == 1)
                {

        //cout<<endl<<BremVector[1000]<<endl;
                        if((fabs(Photon_SC_Eta)) > 1.566 && isLooseMMG != 0 && isMultipleCandidate == 0 && Photon_r9 > 0.95 && Photon_Et > Varmin && Photon_Et <= Varmax)  OneOverKrecoClassicalVector.push_back(mmg_s);

                }


        }
        sort(OneOverKrecoClassicalVector.begin(), OneOverKrecoClassicalVector.end());

        double Min = OneOverKrecoClassicalVector.size() * (pourcentage / 4.0);
        int MinInt = (int) Min;
        double Max = OneOverKrecoClassicalVector.size() * (3.0 * pourcentage / 4.0);
        int MaxInt = (int) Max;

        *MinRange = OneOverKrecoClassicalVector[MinInt];
        *MaxRange = OneOverKrecoClassicalVector[MaxInt];

}



void RangeEstimator3(double pourcentage, TChain * chain, TString temp, int Endcaps, double * MinRange, double * MaxRange)
{
      

        TChain * ReducedChain = (TChain *) chain->CopyTree(temp);
      
        float mmg_s;
        ReducedChain->SetBranchAddress("mmg_s",&mmg_s);

        vector <float> ErecoOverEtrueVector;


        for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
        {
                ReducedChain->GetEntry(ievt);
                ErecoOverEtrueVector.push_back(mmg_s);

        }


        sort(ErecoOverEtrueVector.begin(), ErecoOverEtrueVector.end());

        size_t interval_entries = TMath::Ceil(pourcentage * ErecoOverEtrueVector.size());

        vector<float>::iterator lower = ErecoOverEtrueVector.begin();
        vector<float>::iterator upper = ErecoOverEtrueVector.begin() + interval_entries - 1; 

        double dx = *upper - *lower;

        for(vector<float>::iterator first = lower, last = upper; last < ErecoOverEtrueVector.end(); first++, last++)
        {
                if((*last - *first) < dx)
                {

                        lower = first;
                        upper = last;
                        dx = *upper - *lower;
                }
                      
        }

        *MinRange = *lower;
        *MaxRange = *upper;

}



void SymetricRangeEstimator(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp)
{

	double min = 0.0;
	double max = 0.0;
	string minresult;
	string maxresult;
	std::ostringstream oss;
	std::ostringstream oss2;	
	TString temp2;
	

	for(double sigma = 0.001; sigma < 1.0; sigma += 0.001)
	{
		temp2.Clear();
		oss.str("");
		oss2.str("");
		minresult.clear();
		maxresult.clear();

		min = centralValue - sigma; 
		max = centralValue + sigma;

        	oss << min;
        	minresult = oss.str();
		
                oss2 << max; 
                maxresult = oss2.str();	
		
		TH1D *histo = new TH1D("histo","histo", 200, 0.0, 2.0);
		
		temp2 = temp;
		temp2 += " && mmg_s > ";
		temp2 += minresult;
		temp2 += " && mmg_s < ";
		temp2 += maxresult;

		chain->Draw("mmg_s>>histo", temp2);
		
		if(histo->GetEntries() >= (Entries * pourcentage))
		{
			*MinRange = centralValue - sigma;
			*MaxRange = centralValue + sigma;
			sigma = 1.0;
		}

		histo->Delete();	
	}



}

void SymetricRangeEstimator2(TChain * chain, double lastX, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp)
{

        double min = 0.0;
        double max = 0.0;
        int itermin = 0;
        int itermax = 0;
	int loop1 = 0;
	int loop2 = 0;
        string minresult;
        string maxresult;
        std::ostringstream oss;
        std::ostringstream oss2;
        TString temp2;

	double limite = ((1.0 - pourcentage) / 2.0) * Entries;

        for(double sigma = 0.0; sigma < 2.0; sigma += 0.01)
        {
                temp2.Clear();
                oss.str("");
                oss2.str("");
                minresult.clear();
                maxresult.clear();

                min = sigma;
                max = lastX - sigma;

                oss << min;
                minresult = oss.str();

                oss2 << max;
                maxresult = oss2.str();

                TH1D *histo = new TH1D("histo","histo", 2000, 0.0, 2.0);
                TH1D *histo2 = new TH1D("histo2","histo2", 2000, 0.0, 2.0);

		temp2 = temp;
                temp2 += " && mmg_s < ";
                temp2 += minresult;

                chain->Draw("mmg_s>>histo", temp2);


                if(histo->GetEntries() >= limite && itermin == 0)
		{
                        *MinRange = sigma;
                        itermin = 1;
		}

                temp2.Clear();

                temp2 = temp;
                temp2 += " && mmg_s > ";
                temp2 += maxresult;

		chain->Draw("mmg_s>>histo2", temp2);


                if(histo2->GetEntries() >= limite && itermax == 0)
		{
		        *MaxRange =  lastX - sigma;
                        itermax = 1;
                }

		if (itermin == 1 && itermax == 1) sigma = 2.0;

		

		if((histo->GetEntries() < (limite * 0.1)) && (histo2->GetEntries() < (limite * 0.1)) && loop1 == 0)
		{
			 sigma +=0.05;
		}
		else loop1 = 1;
/*	
                if(loop1 == 1 && (histo->GetEntries() < (limite * 0.9)) && (histo2->GetEntries() < (limite * 0.9)))
		{
                        sigma +=0.01;
                }
		else if(loop1 == 1) loop2 = 1;

		if(loop1 == 1 && loop2 == 1 && (histo->GetEntries() > limite) || (histo2->GetEntries() > limite) && (histo->GetEntries() > (limite * 0.9)) || (histo2->GetEntries() > (limite * 0.9)))
                {
                         sigma +=0.01;
                }

 */
		

                histo->Delete();
                histo2->Delete();
        }


}

void SymetricRangeEstimator3(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp)
{
	

        double min = 0.0;
        double max = 0.0;
	int iterMin = 0;
	int iterMax = 0;
        string minresult;
        string maxresult;
        std::ostringstream oss;
        std::ostringstream oss2;
        TString temp2;


        for(double sigma = 0.001; sigma < 1.0; sigma += 0.001)
        //for(double sigma = 0.01; sigma < 1.0; sigma += 0.01)
	{
		temp2.Clear();
                oss.str("");
                oss2.str("");
                minresult.clear();
                maxresult.clear();

                min = centralValue - sigma; 
                max = centralValue + sigma;

                oss << min; 
                minresult = oss.str();
      
                oss2 << max; 
                maxresult = oss2.str(); 
      
		if(iterMin == 0)
		{
                	TH1D *histo = new TH1D("histo","histo", 200, 0.0, 2.0);
      
                	temp2 = temp;
                	temp2 += " && mmg_s > ";
                	temp2 += minresult;
                	temp2 += " && mmg_s < ";
                	temp2 += centralValue;

                	chain->Draw("mmg_s>>histo", temp2);
			if(histo->GetEntries() >= (Entries * (pourcentage / 2.0)))
                	{
                        	*MinRange = centralValue - sigma;
                        	iterMin = 1;
                	}
			histo->Delete();
		}

		if(iterMax == 0)
                {		

			TH1D *histo2 = new TH1D("histo2","histo2", 200, 0.0, 2.0);
                	
			temp2 = temp;
                	temp2 += " && mmg_s > ";
                	temp2 += centralValue;
                	temp2 += " && mmg_s < ";
                	temp2 += maxresult;
                
			chain->Draw("mmg_s>>histo2", temp2);            
	
			if(histo2->GetEntries() >= (Entries * (pourcentage / 2.0)))
                	{
                        	*MaxRange = centralValue + sigma;
                        	iterMax = 1;
                	}

                	histo2->Delete();
		}
		
		if(iterMin == 1 && iterMax == 1) sigma = 1.0;

        }
		if(iterMin == 0) *MinRange = 0.0;
		if(iterMax == 0) *MaxRange = 1.5;		

}


Double_t chiSquare(RooPlot* plot_, char* pdfname, char* histname, int nFitParam, double* JanChi2, double* DegreesOfFreedom, double* pValue, int* fewBins)
{
  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit

  // Find curve object
  RooCurve* curve = (RooCurve*) plot_->findObject(pdfname, RooCurve::Class());
  //RooCurve* curve = plot_->getCurve(pdfname);  
  //curve->Print();

  if (!curve) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::chiSquare(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot_->findObject(histname, RooHist::Class()) ;
  //RooHist* hist = plot_->getHist(histname);

  if (!hist) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::chiSquare(..) cannot find histogram" << endl ;
    return 0 ;
  }


  Int_t i,np = hist->GetN() ;
  Double_t x,y,/*eyl,eyh,*/ xl,xh ;

  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;

#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif

  Int_t nbin(0) ;

  Double_t chisq(0) ;
  for (i=0 ; i<np ; i++) {   

    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // eyl = hist->GetEYlow()[i] ;
    // eyh = hist->GetEYhigh()[i] ;

    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;

    if(y != 0 && y < 35.0)
    {
    	cout<<endl<<"Trop peu d'entree : "<<y<<" dans le bin : "<<i<<"  >>>Need to reduce the binning for the p-value calculation!"<<endl;
	*fewBins = 1;
	break;
	
    }
    else *fewBins = 0;

    nbin++ ;

    // Integrate function over this bin.
    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    if (avg + avg2 > 0 &&
	(avg2 - avg) / (avg2 + avg) > 0.1) {
      avg = curve->interpolate(x);
    }
    // End of hack around the bug in RooCurve::interpolate

    // JV: Adjust observed and expected number of events for bin width to represent
    // number of events.
    Double_t norm = (xh - xl) / plot_->getFitRangeBinW();
    y *= norm;
    avg *= norm;

    if (avg < 5.) {
      cout << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
			    << ")::chiSquare(..) expectation in bin "
			    << i << " is " << avg << " < 5!" << endl ;
    }

    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf

    // Add pull^2 to chisq
    if (avg != 0) {      
      Double_t resid = y - avg;
      chisq += (resid * resid / avg) ;
    }
  }

  // Return chisq/nDOF 
  *JanChi2 = chisq / (nbin - nFitParam);
  *DegreesOfFreedom = (nbin - nFitParam);
  *pValue =  TMath::Prob(chisq, nbin - nFitParam);

  return chisq / (nbin - nFitParam) ;
}

  
RooHist* residHist(RooPlot* plot_, char *histname, char* curvename, bool normalize, string dossierSauvegardePull, int iteration)
{
  if(normalize == true)
  {
        FILE *fPullsX = fopen(Form("%sPullsX%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsX;
        FILE *fPullsErrorX = fopen(Form("%sPullsErrorX%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsErrorX;
        FILE *fPullsY = fopen(Form("%sPullsY%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsY;
        FILE *fPullsErrorY = fopen(Form("%sPullsErrorY%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsErrorY;
  }

  // Create and return RooHist containing  residuals w.r.t to given curve->
  // If normalize is true, the residuals are normalized by the histogram
  // errors creating a RooHist with pull values


  // Find curve object
  RooCurve* curve = (RooCurve*) plot_->findObject(curvename, RooCurve::Class());
  if (!curve) {
    cout << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::residHist(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot_->findObject(histname, RooHist::Class()) ;
  if (!hist) {
    cout << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::residHist(..) cannot find histogram" << endl ;
    return 0 ;
  }


  // Copy all non-content properties from hist
  RooHist* ret = new RooHist(plot_->getFitRangeBinW()) ;
  if (normalize) {
    ret->SetName(Form("pull_%s_%s", hist->GetName(), curve->GetName())) ;
    ret->SetTitle(Form("Pull of %s and %s", hist->GetTitle(), curve->GetTitle())) ;
  } else {
    ret->SetName(Form("resid_%s_%s", hist->GetName(), curve->GetName())) ;
    ret->SetTitle(Form("Residual of %s and %s", hist->GetTitle(), curve->GetTitle())) ;
  }

  // Determine range of curve
  Double_t xstart, xstop, y ;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0, xstart, y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve&>(curve)->GetPoint(0, xstart, y) ;
  const_cast<RooCurve&>(curve)->GetPoint(curve->GetN()-1, xstop, y) ;
#endif
  // cout << "cit::RooChi2Calculator::residHist dumping curve:\n";
  // for (int i=0; i<curve->GetN(); ++i){
  //   Double_t xi, yi;
  //   curve->GetPoint(i, xi, yi);
  //   printf("i=%d x,y: %.3g, %.3g\n", i, xi, yi);
  // }

  // cout << "cit::RooChi2Calculator::residHist  adding bins with error:\n";

  // Add histograms, calculate Poisson confidence interval on sum value
  for(Int_t i=0 ; i < hist->GetN() ; i++) {
    Double_t x, point;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
    hist->GetPoint(i,x,point) ;
#else
    const_cast<RooHist&>(hist)->GetPoint(i,x,point) ;
#endif
    Double_t xl = x - hist->GetErrorXlow(i);
    Double_t xh = x + hist->GetErrorXhigh(i);

    // Only calculate pull for bins inside curve range
    if (xl < xstart || xstop < xh) continue ;

    Double_t norm = (xh - xl) / plot_->getFitRangeBinW();
    point *= norm;

    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    Double_t yexpected;
    if (avg + avg2 > 0 && (avg2 - avg) / (avg2 + avg) > 0.1) {
      yexpected = curve->interpolate(x);
    } else {
      yexpected = avg;
    }
    // End of hack around the bug in RooCurve::interpolate

    // Correct the expected number of events in this bin for the non-uniform
    // bin width.
    yexpected *= norm;

    Double_t yy = point - yexpected;
    // Normalize to the number of events per bin taking into account
    // variable bin width.
    Double_t dy = TMath::Sqrt(yexpected);
    if (normalize) {
	if (dy==0.) {
	  cout << "cit::RooChi2Calculator::residHist(histname ="
               << hist->GetName() << ", ...) WARNING: point "
               << i << " has zero error, setting residual to zero"
               << endl ;
	  yy=0 ;
	  dy=0 ;
	} else {
	  yy /= dy;
	  dy = 1.;
	}
    }
    // printf("bin=%3d n=%5.3g nu=%5.3g x=%5.3g .. %5.3g y=%5.3g +/- %5.3g "
    //	   "norm=%5.3g\n", i, point, yexpected, xl, xh, yy, dy, norm);
    ret->addBinWithError(x,yy,dy,dy);
  
    if(normalize == true)
    {
        ofstream monFluxPullsX(Form("%sPullsX%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsX << yy <<endl;
        monFluxPullsX.close();

        ofstream monFluxPullsErrorX(Form("%sPullsErrorX%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsErrorX << dy <<endl;
        monFluxPullsErrorX.close();

        ofstream monFluxPullsY(Form("%sPullsY%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsY << x <<endl;
        monFluxPullsY.close();

        ofstream monFluxPullsErrorY(Form("%sPullsErrorY%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsErrorY << dy <<endl;
        monFluxPullsErrorY.close();

    }

  }
  return ret;
}


int NbLignesFichier(string fichier)
{
        ifstream in(fichier.c_str()); //Ouverture en mode lecture de fichier

        string ligne; //Création d'une chaine de caractere
        int nbLignes = 0;

        while(std::getline(in, ligne)) nbLignes++;

        in.close(); //On ferme le fichier

        return nbLignes;
}




