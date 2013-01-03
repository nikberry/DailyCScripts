#include "TFile.h"
#include "TH1.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THStack.h"
#include <string.h>
#include <iostream>
#include <iomanip>
#include "tdrstyle.C"

void doPlotsMuon();
TH1D* getSample(TString sample, double weight);

double lumi = 5500;
//stuff to choose
bool logPlot = false; //true for log plot
int rebinFact = 5;

//isolation selection
//TString Isolation = "QCD No Iso/";
TString Isolation = "Ref selection/";
//TString Isolation = "QCD mu+jets PFRelIso/";
//TString Isolation = "QCD non iso mu+jets/";

// number of btags
TString Nbtags = "2btags";  //standard  "2btags" , qcd "0btag"

bool inclZ = false;
bool inclW = false;

//choose object
TString Obj = "Muon/";
//TString Obj = "MET/";

//muon variables
const int N = 6;
TString Variable;
TString Variables[N] = {"muon_AbsEta_", "muon_eta_", "muon_pfIsolation_04_", "muon_pT_", "muon_phi_", "muon_dB_"};
double MinXs[N] = {0,-2.6 ,0 , 0, -3.2, 0};
double MaxXs[N] = {2.6,2.6 ,0.15 , 400,3.2 , 0.3};
TString XTitles[N] = {"#left|#eta#right|_{#mu}", "#eta_{#mu}", "RelIso_{#mu}", "p_{T}(#mu) [GeV]", "#phi_{#mu}", "dB_{#mu}"};

//met variables
//TString Variable = "patType1CorrectedPFMet/MET_";


void doPlotsMuon(){
setTDRStyle();

//loop over variables
for(int i = 0; i<N; i++){
double MinX = MinXs[i];
double MaxX = MaxXs[i];
Variable = Variables[i];
TString Xtitle = XTitles[i];

//Data
TH1D* data = getSample("SingleMu", 1);

//MC
TH1D* tt = getSample("TTJet", lumi*225.2/6920475);

TH1D* wjets = getSample("WJetsToLNu", lumi*37509/57708550);
TH1D* w1jets = getSample("W1Jet", lumi*5400.0/23140779);
TH1D* w2jets = getSample("W2Jets", lumi*1750.0/34041404);
TH1D* w3jets = getSample("W3Jets", lumi*519.0/15536443);
TH1D* w4jets = getSample("W4Jets", lumi*214.0/13370904);

TH1D* zjets = getSample("DYJetsToLL", lumi*5745.25/30457954);
TH1D* z1jets = getSample("DY1JetsToLL", lumi*561.0/24042904);
TH1D* z2jets = getSample("DY2JetsToLL", lumi*181.0/21835749);
TH1D* z3jets = getSample("DY3JetsToLL", lumi*51.1/11010628);
TH1D* z4jets = getSample("DY4JetsToLL", lumi*23.04/6391785);

TH1D* qcd = getSample("QCD_Pt_20_MuEnrichedPt_15", lumi*34679.3/8500505);

TH1D* top_t = getSample("T_t-channel", lumi*56.4/3757707);
TH1D* top_tw = getSample("T_tW-channel", lumi*11.1/497395);
TH1D* top_s = getSample("T_s-channel", lumi*3.79/249516);
TH1D* tbar_t = getSample("Tbar_t-channel", lumi*30.7/1934817);
TH1D* tbar_tw = getSample("Tbar_tW-channel", lumi*11.1/493239);
TH1D* tbar_s = getSample("Tbar_s-channel", lumi*1.76/139948);

THStack *hs = new THStack("hs","test");
  hs->Add(qcd);
  
  if(inclZ == true){
  hs->Add(zjets);
  }else{
  hs->Add(z1jets);
  hs->Add(z2jets);
  hs->Add(z3jets);
  hs->Add(z4jets);  
  }
  
  if(inclW == true){
  hs->Add(wjets);
  }else{
  hs->Add(w1jets);
  hs->Add(w2jets);
  hs->Add(w3jets);
  hs->Add(w4jets);  
  }
  
  hs->Add(top_t);
  hs->Add(top_tw);
  hs->Add(top_s);
  hs->Add(tbar_t);
  hs->Add(tbar_tw);
  hs->Add(tbar_s);
  
  hs->Add(tt);

  //draw histos to files
  TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);
		
  hs->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.3);

  hs->Draw();
  data->Draw("E same");
  data->SetMarkerStyle(20);
  
  hs->GetXaxis()->SetLimits(MinX, MaxX);
  hs->GetXaxis()->SetTitle(Xtitle); hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitle("Number of Events");hs->GetYaxis()->SetTitleSize(0.05);
  
  
  	TLegend *tleg2;
	tleg2 = new TLegend(0.7,0.7,0.8,0.9);
	tleg2->SetTextSize(0.04);
	tleg2->SetBorderSize(0);
	tleg2->SetFillColor(10);
	tleg2->AddEntry(data , "2012 data", "lpe");
	tleg2->AddEntry(tt , "t#bar{t}", "lf");
	tleg2->AddEntry(top_t, "single top", "lf");
	tleg2->AddEntry(wjets , "w+jets", "lf");
	tleg2->AddEntry(zjets , "z+jets", "lf");
	tleg2->AddEntry(qcd , "QCD", "lf");
	
	//tleg2->AddEntry(singtEff, "single-t"      , "l");
	//tleg2->AddEntry(singtwEff, "single-tW"      , "l");
 	tleg2->Draw("same");	

  if(logPlot ==true){
  c1->SetLogy();
  }	
  
  TString plotName("plots/Control/Muon/");
  
  if(logPlot ==true){
    plotName += Variable+"Test_Log";
    plotName += Nbtags+".png";
    
  }else{
    plotName += Variable+"Test";  
    plotName += Nbtags+".png";
  }
 
 
  c1->SaveAs(plotName);
  delete c1;
  
  }
  	
}


TH1D* getSample(TString sample, double weight){
	TString dir = "rootFiles/";
	TFile* file = new TFile(dir + sample + "_10000pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");
	
	TH1D* plot = (TH1D*) file->Get("TTbarPlusMetAnalysis/MuPlusJets/"+Isolation+Obj+Variable+Nbtags);

        if(sample == "TTJet"){
	plot->SetFillColor(kRed+1);
        plot->SetLineColor(kRed+1);
	}else if(sample == "WJetsToLNu" || sample == "W1Jet" || sample == "W2Jets"|| sample == "W3Jets"|| sample == "W4Jets"){
	plot->SetLineColor(kGreen-3);	  
  	plot->SetFillColor(kGreen-3);
	}else if(sample == "DYJetsToLL" || sample == "DY1JetsToLL" || sample == "DY2JetsToLL" || sample == "DY3JetsToLL" || sample == "DY4JetsToLL"){
	plot->SetFillColor(kAzure-2);
	plot->SetLineColor(kAzure-2);
	}else if(sample == "QCD_Pt_20_MuEnrichedPt_15"){
	plot->SetFillColor(kYellow);
	plot->SetLineColor(kYellow);
	}else if(sample == "T_t-channel" || sample == "T_tW-channel" || sample == "T_s-channel" || sample == "Tbar_t-channel" || sample == "Tbar_tW-channel" || sample == "Tbar_s-channel"){
	plot->SetFillColor(kMagenta);
	plot->SetLineColor(kMagenta);
	}
	
	plot->Scale(weight);
	plot->Rebin(rebinFact);
	
	return plot;

}
