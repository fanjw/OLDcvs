#include <TApplication.h>
#include <TLegend.h>

typedef TH1* pTH1;
TString eff_cmp_efficiency()
{
   gROOT->ProcessLine(".x setTDRStyle.C");
  

   TString const usage = "usage: root eff_cmp.C xxx.root [yyy.root ...]";

  TCanvas *c1 = new TCanvas("reconstruction","reconstruction");
  c1->cd();
  c1->SetGridx(1);
  c1->SetGridy(1);

   double ymin = 1;
   pTH1* rejh = new pTH1[ gApplication->Argc() ];
   double* eff90 = new double[ gApplication->Argc() ];
   int nHist = 0;
   for(int i=1;i<gApplication->Argc();++i)
   {
      TString n = gApplication->Argv(i);
      if( ! n.EndsWith(".root") )continue;

      TFile * file = TFile::Open( n );
      if( file==NULL || file->IsZombie() )continue;

      //TH1* hrej = dynamic_cast<TH1*>( file->Get("Method_MLP/MLP/MVA_MLP_rejBvsS") );
      TH1* hrej = dynamic_cast<TH1*>( file->Get("Method_BDT/BDT/MVA_BDT_rejBvsS") );
      //TH1* hrej = dynamic_cast<TH1*>( file->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS") );

      if( hrej==NULL )continue;         


      rejh[ nHist ] = hrej;
      eff90[ nHist ]= (hrej->GetBinContent( hrej->FindBin( 0.9 ) ) + hrej->GetBinContent( hrej->FindBin( 0.9 )-1) )/2.0;
      TString b = n.Replace(0,5,""); 
      TString fn = b.ReplaceAll(".root","");
      hrej->SetTitle( Form("%.4g: %s", eff90[nHist], fn.Data() ));
      hrej->SetName(fn.Data());
      //cout << fn.Data() <<endl;
      //hrej->SetTitle( Form("%.4g: %s", eff90[nHist], n.Data() ));
      //hrej->SetTitle( Form("%.4g (0.90 sig. eff.): %s", eff90[nHist], n.Data() ));
      ++ nHist;

      ymin = TMath::Min( ymin, hrej->GetBinContent( hrej->FindBin( 0.95 ) ) );
   }


   int* index = new int[ nHist ];
   TMath::Sort( nHist, eff90, index, true );
   TH1* hmin = rejh[ index[nHist-1] ];
   TString MinNom = hmin->GetName();
   TString RealMinNom1 = MinNom.Replace(0,8,"");
   int MinNomNumRest1 = RealMinNom1.Length();
   TString RealMinNom = RealMinNom1.Replace(MinNomNumRest1-3,3,"");
   cout<<RealMinNom<<endl;
   TString MinNom1 = hmin->GetName();
   TString MinNomNum = MinNom1.Replace(0,4,"");
   int MinNomNumRest = MinNomNum.Length();
   TString RealMinNomNum = MinNomNum.Replace(2,MinNomNumRest-2,"");
   cout << RealMinNomNum <<endl;
   cout << "\n\n\n\n\n" << endl; 
   cout << "MIN:";
   cout << RealMinNomNum<<" ";
   cout << RealMinNom.Data() << endl; 


   //TH2D *frame = new TH2D("frame_xx", "", 100, 0.7, 1, 100, ymin-1.0, 1 );
   TH2D *frame = new TH2D("frame_xx", "", 100, 0.5, 1, 100, 0.6, 1 );

   frame->GetXaxis()->SetTitle("Signal efficiency");  
   frame->GetYaxis()->SetTitle("Background rejection");  
   frame->GetXaxis()->CenterTitle(1);
   frame->GetYaxis()->CenterTitle(1);


   gStyle->SetOptStat(0);
   frame->Draw();
   TLegend* lg = new TLegend(0.18, 0.18, 0.65, 0.55);
   Color_t color = 1;  int LineStyle=1;
   cout << "\n\n\n\n\n";
   for( int i=0; i < nHist; ++i )
   {
      TH1* hrej = rejh[ index[i] ];
      hrej->SetLineColor( color );
      hrej->SetLineWidth(2);
      hrej->SetLineStyle( LineStyle );

      ++color;
      if (color == 5) color++;
      //if (color == 5 || color == 10 || color == 11) color++;
      //if (color == 18 ) color +=2;

      if(i!=0 && i%5==0) {
        LineStyle++;
        color = 1; 
      }

      lg->AddEntry( hrej, NULL, "L" );
      hrej->DrawCopy( "csame" );
      cout << hrej->GetTitle() << endl;
   }

   frame->Draw("sameaxis");
   lg->SetFillColor(0);
   lg->SetTextFont( lg->GetTextFont()/10*10 );
   lg->Draw("same");

   //c1->Print("MVA_Compare.ps");

   delete[] eff90;
   delete[] rejh;
   delete[] index;
   cout << "\n\n\n\n\n";
   return "done.";
}
