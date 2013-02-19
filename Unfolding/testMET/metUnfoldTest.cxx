//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "tdrstyle.C"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#include "RooUnfoldBinByBin.h"
#endif

//==============================================================================
// Example Unfolding
//==============================================================================

void metUnfoldTest()
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif

setTDRStyle();

double ttbarXsect =225.197;
int kval= 3;

//MET will need choice of variable at the top
 TString Variable ="_MET";
 int Nbins = 6;
 double width[6] = {25, 20, 25, 30, 50, 100};
 double xbins[7] = {0,25,45,70,100,150, 250}; 
 double xbinsTheo[7]= {3,28,48,73,103,153, 253};
 double xbinsMeas[7]= {-3,22,42,67,97,147, 247};
 TString folder = "unfolding_MET_analyser_muon_channel_patType1CorrectedPFMet/";
 TString Xtitle = "MET [GeV]";
 
//ST will need choice of variable at the top
//  TString Variable ="_ST";
//  int Nbins = 8;
//  double width[8] = {150,100,100,100,100,200,500,500};
//  double xbins[9] = {1,150,250,350,450,550,750,1250, 1750}; 
//  double xbinsTheo[9]= {1+10,150+10,250+10,350+10,450+10,550+10,750+10,1250+10, 1750+10};
//  double xbinsMeas[9]= {1-10,150-10,250-10,350-10,450-10,550-10,750-10,1250-10, 1750-10};
//  TString folder = "unfolding_ST_analyser_muon_channel_patType1CorrectedPFMet/";
//  TString Xtitle = "ST [GeV]";
 
 
  // Get histograms  
  cout << "==================================== READ =====================================" << endl;
  TString dir = "rootFiles/";
  TFile* unf_file = new TFile(dir + "unfolding_2012.root");
  TFile* meas_file = new TFile(dir + "diffResults"+Variable+".root");

  //folder for 2011 unfolding
  //TString folder = "unfoldingAnalyserMuonChannel/"
    
  TH2D* hResp = (TH2D*) unf_file->Get(folder+"response_AsymBins");

  TH1D *hRespMeas = (TH1D*) unf_file->Get(folder+"measured_AsymBins");
  //truth with extra for inefficiency 
  TH1D* hTrue = (TH1D*) unf_file->Get(folder+"truth_AsymBins");

  cout << "==================================== TRAIN ====================================" << endl;
  RooUnfoldResponse response (hRespMeas, hTrue, hResp, "response", "response");

  cout << "==================================== TEST =====================================" << endl;
  TH1D* hTrue_Assym= (TH1D*) unf_file->Get(folder+"truth_AsymBins");
  //TH1D* hMeas= (TH1D*) unf_file->Get(folder+"measured_AsymBins");    // this is used for closure test
  TH1D* hMeas= (TH1D*) meas_file->Get("central_dir/central_ttbar_fit");
  
  cout << "==================================== UNFOLD ===================================" << endl;
//RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
  RooUnfoldSvd     unfold (&response, hMeas, kval);    // OR
//RooUnfoldTUnfold unfold (&response, hMeas);
//RooUnfoldBinByBin unfold (&response, hMeas);

  TH1D* hReco= (TH1D*) unfold.Hreco();
   
   //differential histo
   TH1D *meas  = new TH1D("meas", "", Nbins, xbinsMeas); 
   TH1D *theo  = new TH1D("theor", "", Nbins, xbinsTheo); 
   TH1D *reco  = new TH1D("nominal", "", Nbins, xbins); 
   
   for(int i = 0; i < Nbins; i++){
   meas->SetBinContent(i+1,hMeas->GetBinContent(i+1));
   theo->SetBinContent(i+1,hTrue_Assym->GetBinContent(i+1));
   reco->SetBinContent(i+1,hReco->GetBinContent(i+1));
   meas->SetBinError(i+1,hMeas->GetBinError(i+1));
   theo->SetBinError(i+1,hTrue_Assym->GetBinError(i+1));
   reco->SetBinError(i+1,hReco->GetBinError(i+1));   
   }

   //pre fit used for    
   TH1D* tt_tot= (TH1D*) meas_file->Get("central_dir/central_ttbar_prefit");
   
   reco->Scale((hMeas->Integral()/(reco->Integral()*tt_tot->Integral()))*ttbarXsect);
   meas->Scale((1./tt_tot->Integral())*ttbarXsect);
   theo->Scale((1./hTrue_Assym->Integral())*ttbarXsect);

   for(int i = 0; i < Nbins; i++){
   cout << "meas: " << meas->GetBinContent(i+1) << endl;
   cout << "theor: " << theo->GetBinContent(i+1) << endl;
   cout << "reco: " << reco->GetBinContent(i+1) << endl;  
   }

  unfold.PrintTable (cout, hTrue_Assym);
     
  TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);

  reco->SetMinimum(0.0);
  reco->SetMaximum(theo->GetBinContent(hTrue_Assym->GetMaximumBin())*1.3);
  reco->Draw();
  
  reco->GetXaxis()->SetTitle(Xtitle); reco->GetXaxis()->SetTitleSize(0.05);
  reco->GetYaxis()->SetTitle("#partial #sigma [pb]");reco->GetYaxis()->SetTitleSize(0.05);
  
  meas->SetLineColor(kRed);
  meas->SetMarkerColor(kRed);
  meas->Draw("SAME");
  theo->SetLineColor(kGreen);
  theo->SetMarkerColor(kGreen);
  theo->Draw("SAME");
    	
	TLegend *tleg2;
	tleg2 = new TLegend(0.7,0.5,0.8,0.7);
	tleg2->SetTextSize(0.04);
	tleg2->SetBorderSize(0);
	tleg2->SetFillColor(10);
	tleg2->AddEntry(reco , "unfolded", "lpe");
 	tleg2->AddEntry(meas , "measured", "lpe");
	tleg2->AddEntry(theo , "truth", "lpe");
 	
	tleg2->Draw("same");	
	
  c1->SaveAs("plots/bin/Nevents"+Variable+".png");

   TCanvas *c2 = new TCanvas("Plot","Plot",900, 600);

   double measXsect = meas->Integral();
   double theoXsect = theo->Integral();
   double recoXsect = reco->Integral();
   
   for(int i = 0; i < Nbins; i++){
   meas->SetBinContent(i+1,meas->GetBinContent(i+1)/measXsect); 
   theo->SetBinContent(i+1,theo->GetBinContent(i+1)/theoXsect); 
   reco->SetBinContent(i+1,reco->GetBinContent(i+1)/recoXsect);
   meas->SetBinError(i+1,meas->GetBinError(i+1)/measXsect); 
   theo->SetBinError(i+1,theo->GetBinError(i+1)/theoXsect); 
   reco->SetBinError(i+1,reco->GetBinError(i+1)/recoXsect);
   }

  reco->SetMinimum(0.0);
  reco->SetMaximum(theo->GetBinContent(hTrue_Assym->GetMaximumBin())*1.3);
  reco->Draw();
  
  reco->GetXaxis()->SetTitle(Xtitle); reco->GetXaxis()->SetTitleSize(0.05);
  reco->GetYaxis()->SetTitle("#frac{1}{#sigma} #partial #sigma");reco->GetYaxis()->SetTitleSize(0.05);
  
  meas->SetLineColor(kRed);
  meas->SetMarkerColor(kRed);
  meas->Draw("SAME");
  theo->SetLineColor(kGreen);
  theo->SetMarkerColor(kGreen);
  theo->Draw("SAME");
 	
	tleg2->Draw("same");	
	
  c2->SaveAs("plots/bin/normalised"+Variable+".png");

   TCanvas *c3 = new TCanvas("Plot","Plot",900, 600);
   
   for(int i = 0; i < Nbins; i++){
   meas->SetBinContent(i+1,meas->GetBinContent(i+1)/width[i]); 
   theo->SetBinContent(i+1,theo->GetBinContent(i+1)/width[i]); 
   reco->SetBinContent(i+1,reco->GetBinContent(i+1)/width[i]);
   meas->SetBinError(i+1,meas->GetBinError(i+1)    /width[i]); 
   theo->SetBinError(i+1,theo->GetBinError(i+1)    /width[i]); 
   reco->SetBinError(i+1,reco->GetBinError(i+1)    /width[i]);
   }

  reco->SetMinimum(0.0);
  reco->SetMaximum(theo->GetBinContent(hTrue_Assym->GetMaximumBin())*1.3);
  reco->Draw();
  
  reco->GetXaxis()->SetTitle(Xtitle); reco->GetXaxis()->SetTitleSize(0.05);
  reco->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{#partial #sigma}{#partial MET} 				[GeV^{-1}]");reco->GetYaxis()->SetTitleSize(0.05);
  
  meas->SetLineColor(kRed);
  meas->SetMarkerColor(kRed);
  meas->Draw("SAME");
  theo->SetLineColor(kGreen);
  theo->SetMarkerColor(kGreen);
  theo->Draw("SAME");
 	
	tleg2->Draw("same");	
	
  c3->SaveAs("plots/bin/Nevents"+Variable+".png");
 
 for(int i = 0; i < Nbins; i++){
 cout <<  "bin: " << i+1 << " , val: " << reco->GetBinContent(i+1) << endl;
 }
 
  
}

#ifndef __CINT__
int main () { metUnfoldTest(); return 0; }  // Main program when run stand-alone
#endif
