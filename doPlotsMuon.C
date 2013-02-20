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
#include <sstream>
#include <iomanip>
#include "tdrstyle.C"
#include "getSamples.C"

void doPlotsMuon();
//TH1D* getSample(TString sample, double weight);
//TH1D* getQCD();


double lumi = 5800;
//stuff to choose
bool logPlot = false; //true for log plot
int rebinFact = 10;

//isolation selection
//TString Isolation = "QCD No Iso/";
TString Isolation = "Ref selection/";
//TString Isolation = "QCD mu+jets PFRelIso/";
//TString Isolation = "QCD non iso mu+jets/";

// number of btags
TString Nbtags = "2btags";  //standard  "2btags" , qcd "0btag"

bool inclZ = false;
bool inclW = false;
bool inclQ = false;
//choose object
TString Obj = "Muon/";
//TString Obj = "MET/";

//TString Systematic = "central";

//muon variables
const int N = 8;
TString Variable;
TString Variables[N] = {"muon_AbsEta_", "muon_eta_",  "muon_pT_", "muon_phi_", "muon_dB_", "muon_pfIsolation_03_", "muon_pfIsolation_04_", "muon_pfIsolation_04_"};
double MinXs[N] = {0, -2.6 , 0, -3.2, 0, 0, 0 , 0};
double MaxXs[N] = {2.6,2.6 , 400,3.2 , 0.3, 5, 5, 5};
TString XTitles[N] = {"#left|#eta#right|_{#mu}", "#eta_{#mu}", "p_{T}(#mu) [GeV]", "#phi_{#mu}", "dB_{#mu}", "PF Relative Isolation", "PF Relative Isolation", "PF Relative Isolation"};

void doPlotsMuon(){
setTDRStyle();

//loop over Systematics
const int S = 9;
TString Systematic;
TString Systematics[S] = {"BJet_down", "BJet_up", "central", "JES_down", "JES_up", "LightJet_down", "LightJet_up", "PU_down", "PU_up"};
for(int n = 0; n<S; n++){
Systematic = Systematics[n];

//loop over variables
for(int i = 0; i<N; i++){
double MinX = MinXs[i];
double MaxX = MaxXs[i];
Variable = Variables[i];
TString Xtitle = XTitles[i];

//Data
TH1D* data = getSample("SingleMu", 1, Obj, Variable, Isolation, rebinFact, "central"); //"central"

//MC
TH1D* tt = getSample("TTJet", lumi*225.2/6920475, Obj, Variable, Isolation, rebinFact, Systematic);

TH1D* wjets = getSample("W1Jet", lumi*37509/57708550, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* w1jets = getSample("W1Jet", lumi*5400.0/23140779, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* w2jets = getSample("W2Jets", lumi*1750.0/34041404, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* w3jets = getSample("W3Jets", lumi*519.0/15536443, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* w4jets = getSample("W4Jets", lumi*214.0/13370904, Obj, Variable, Isolation, rebinFact, Systematic);

TH1D* zjets = getSample("DY1JetsToLL", lumi*5745.25/30457954, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* z1jets = getSample("DY1JetsToLL", lumi*561.0/24042904, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* z2jets = getSample("DY2JetsToLL", lumi*181.0/21835749, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* z3jets = getSample("DY3JetsToLL", lumi*51.1/11010628, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* z4jets = getSample("DY4JetsToLL", lumi*23.04/6391785, Obj, Variable, Isolation, rebinFact, Systematic);

TH1D* qcd = getSample("QCD_Pt-15to20_MuEnrichedPt5",     lumi*34679.3/8500505, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd1 = getSample("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd2 = getSample("QCD_Pt-20to30_MuEnrichedPt5",   lumi*2.87e8 * 0.0065/8486893, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd3 = getSample("QCD_Pt-30to50_MuEnrichedPt5",   lumi*6.609e7 * 0.0122/8928999, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd4 = getSample("QCD_Pt-50to80_MuEnrichedPt5",   lumi*8082000.0 * 0.0218/7256011, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd5 = getSample("QCD_Pt-80to120_MuEnrichedPt5",  lumi*1024000.0 * 0.0395/9030624, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd6 = getSample("QCD_Pt-120to170_MuEnrichedPt5", lumi*157800.0 * 0.0473/8500505, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd7 = getSample("QCD_Pt-170to300_MuEnrichedPt5", lumi*34020.0 * 0.0676/7662483, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd8 = getSample("QCD_Pt-300to470_MuEnrichedPt5", lumi*1757.0 * 0.0864/7797481, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd9 = getSample("QCD_Pt-470to600_MuEnrichedPt5", lumi*115.2 * 0.1024/2995767, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd10 = getSample("QCD_Pt-800to1000_MuEnrichedPt5",lumi*3.57 * 0.1033/4047142, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* qcd11 = getSample("QCD_Pt-1000_MuEnrichedPt5",     lumi*0.774 * 0.1097/3807263, Obj, Variable, Isolation, rebinFact, Systematic);

  if(inclQ == true){
  qcd->Add(qcd);
  }else{
  qcd1->Add(qcd2);
  qcd1->Add(qcd3);
  qcd1->Add(qcd4);
  qcd1->Add(qcd5);
  qcd1->Add(qcd6);
  qcd1->Add(qcd7);
  qcd1->Add(qcd8);
  qcd1->Add(qcd9);
  qcd1->Add(qcd10);
  qcd1->Add(qcd11);
  }
  
TH1D* qcd_data = getQCD(Obj, Variable, rebinFact);
qcd_data->Scale(qcd1->Integral());

TH1D* top_t = getSample("T_t-channel", lumi*56.4/3757707, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* top_tw = getSample("T_tW-channel", lumi*11.1/497395, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* top_s = getSample("T_s-channel", lumi*3.79/249516, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* tbar_t = getSample("Tbar_t-channel", lumi*30.7/1934817, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* tbar_tw = getSample("Tbar_tW-channel", lumi*11.1/493239, Obj, Variable, Isolation, rebinFact, Systematic);
TH1D* tbar_s = getSample("Tbar_s-channel", lumi*1.76/139948, Obj, Variable, Isolation, rebinFact, Systematic);

THStack *hs = new THStack("hs","test");
  if(inclQ == true){
  hs->Add(qcd_data);
  }else{
  hs->Add(qcd_data);
  }
  
  
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

	TText* textPrelim = doPrelim(0.17,0.96);
	textPrelim->Draw();

  if(logPlot ==true){
  c1->SetLogy();
  }	
  
  TString plotName("plots/Control/Muon2/");
  
  if(logPlot ==true){
    plotName += Variable+"Log";
    plotName += Nbtags+".pdf";
    
  }else{
    plotName += "/"+Systematic+"/"+Variable;  
    plotName += Nbtags+".pdf";
  }
 
 
  c1->SaveAs(plotName);
  delete c1;
  
  }
  	
}

}
