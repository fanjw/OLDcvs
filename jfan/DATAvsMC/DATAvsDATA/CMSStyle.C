//##########################################################################
// Initialization code executed at the start of a ROOT session.
//
// File: $Id: CMSStyle.C,v 1.1 2012/09/14 07:58:32 jfan Exp $
//##########################################################################

//#include <iostream>
//using namespace std;
 


void CMSstyle(){
  // use the 'plain' style for plots (white backgrounds, etc)
  //  cout << "...using style 'Plain'" << endl;
  gROOT->SetStyle("Plain");
  
  // Create the 'CMS' style for approved plots. Note that this style may need
  // some fine tuning in your macro depending on what you are plotting, e.g.
  //
  //  gStyle->SetMarkerSize(0.75);  // use smaller markers in a histogram with many bins
  //  gStyle->SetTitleOffset(0.65,"y");  // bring y axis label closer to narrow values
  
  TStyle *cmsStyle= new TStyle("CMS","CMS approved plots style");
	 
  // use plain black on white colors
  cmsStyle->SetFrameBorderMode(0);
  cmsStyle->SetCanvasBorderMode(0);
  cmsStyle->SetPadBorderMode(0);
  cmsStyle->SetPadColor(0);
  cmsStyle->SetCanvasColor(0);
  cmsStyle->SetTitleColor(1);
  cmsStyle->SetStatColor(0);
  cmsStyle->SetFrameFillColor(0);
  cmsStyle->SetLegendBorderSize(0);

  // set the paper & margin sizes
  cmsStyle->SetPaperSize(20,26);
  cmsStyle->SetPadTopMargin(0.05);
  cmsStyle->SetPadRightMargin(0.10);
  cmsStyle->SetPadBottomMargin(0.17);
  cmsStyle->SetPadLeftMargin(0.17);

  // use large Times-Roman fonts
  cmsStyle->SetTextFont(132);
  cmsStyle->SetTextSize(0.08);
  cmsStyle->SetLabelFont(132,"x");
  cmsStyle->SetLabelFont(132,"y");
  cmsStyle->SetLabelFont(132,"z");
  cmsStyle->SetTitleFont(132,"x");
  cmsStyle->SetTitleFont(132,"y");
  cmsStyle->SetTitleFont(132,"z");
  cmsStyle->SetLabelSize(0.05,"x");
  cmsStyle->SetTitleSize(0.06,"x");
  cmsStyle->SetLabelSize(0.05,"y");
  cmsStyle->SetTitleSize(0.06,"y");
  cmsStyle->SetLabelSize(0.05,"z");
  cmsStyle->SetTitleSize(0.06,"z");

  // use bold lines and markers
  cmsStyle->SetMarkerStyle(8);
  cmsStyle->SetHistLineWidth(1.85);
  cmsStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // do not display any of the standard histogram decorations
  cmsStyle->SetOptTitle(1);
  cmsStyle->SetOptStat(1);
  cmsStyle->SetOptFit(1);

  // put tick marks on top and RHS of plots
  cmsStyle->SetPadTickX(1);
  cmsStyle->SetPadTickY(1);

//  cout << endl << "    For approved plots use: gROOT->SetStyle(\"CMS\");"
//       << endl << "  To add a CMS label use: CMSLabel();"
//       << endl << endl;

  // restore the plain style
  gROOT->SetStyle("Plain");
  gROOT->SetStyle("CMS");
  return ;
}
