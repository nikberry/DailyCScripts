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

void doPlotsJets();
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

//choose object
TString Obj = "Jets/";

//muon variables
const int N = 8;
TString Variable;
TString Variables[N] = {"all_jet_pT_", "all_jet_phi_", "all_jet_eta_", "jet1_pT_", "jet2_pT_", "jet3_pT_", "jet4_pT_", "jet5_pT_"};
double MinXs[N] = {0, -3.4, -2.6,0,0,0,0,0};
double MaxXs[N] = {500, 3.4, 2.6,500,500,500,500,500};
TString XTitles[N] = {"p_{T} (Jets)", "#Phi (Jets)", "#eta (Jets)""p_{T} (Jet 1)","p_{T} (Jet 2)","p_{T} (Jet 3)","p_{T} (Jet 4)","p_{T} (Jet 5)"};


void doPlotsJets(){
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
TH1D* zjets = getSample("DYJetsToLL", lumi*5745.25/30457954);
TH1D* qcd = getSample("QCD_Pt-120to170_Mu", lumi*34679.3/8500505);


THStack *hs = new THStack("hs","test");
  hs->Add(qcd);
  hs->Add(zjets);
  hs->Add(wjets);
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
	tleg2->AddEntry(wjets , "w+jets", "lf");
	tleg2->AddEntry(zjets , "z+jets", "lf");
	//tleg2->AddEntry(singtEff, "single-t"      , "l");
	//tleg2->AddEntry(singtwEff, "single-tW"      , "l");
 	tleg2->Draw("same");	

  if(logPlot ==true){
  c1->SetLogy();
  }	
  
  TString plotName("plots/Control/"+Obj);
  
  if(logPlot == true){
    plotName += Variable+"_Log";
    plotName += Nbtags+".png";
    
  }else{
    plotName += Variable;  
    plotName += Nbtags+".png";
  }
 
 
  c1->SaveAs(plotName);
  delete c1;
  
  }
  	
}


TH1D* getSample(TString sample, double weight){
	TString dir = "rootFiles/";
	TFile* file = new TFile(dir + sample + "_5050pb_PFElectron_PFMuon_PF2PATJets_PFMET_FULL.root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");
	
	TH1D* plot = (TH1D*) file->Get("TTbarPlusMetAnalysis/MuPlusJets/"+Isolation+Obj+Variable+Nbtags);

        if(sample == "TTJet"){
	plot->SetFillColor(kRed+1);
        plot->SetLineColor(kRed+1);
	}else if(sample == "WJetsToLNu"){
	plot->SetLineColor(kGreen-3);	  
  	plot->SetFillColor(kGreen-3);
	}else if(sample == "DYJetsToLL"){
	plot->SetFillColor(kAzure-2);
	plot->SetLineColor(kAzure-2);
	}else if(sample == "QCD_Pt-120to170_Mu"){
	plot->SetFillColor(kYellow);
	plot->SetLineColor(kYellow);
	}
	
	plot->Scale(weight);
	plot->Rebin(rebinFact);
	
	return plot;

}

