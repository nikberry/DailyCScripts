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

void allCutFlow();
TH1D* getSample(TString sample, double weight);
TText* doPrelim(float x, float y);

double lumi = 5800;
//stuff to choose
bool logPlot = true; //true for log plot
int rebinFact = 1;

//isolation selection
//TString Isolation = "QCD No Iso/";
TString Isolation = "Ref selection/";
//TString Isolation = "QCD mu+jets PFRelIso/";
//TString Isolation = "QCD non iso mu+jets/";

// number of btags
TString Nbtags = "";  //standard  "2btags" , qcd "0btag"

bool inclZ = false;
bool inclW = false;
bool inclQ = false;
//choose object
TString Obj = "Muon/";
//TString Obj = "MET/";

//muon variables
const int N = 2;
TString Variable;
TString Variables[N] = {"TTbarMuPlusJetsRefSelection", "TTbarMuPlusJetsRefSelectionUnweighted"};
double MinXs[N] = {0,-2.6 };
double MaxXs[N] = {2.6,2.6 };
TString XTitles[N] = {"Cuts", "Cuts"};

//met variables
//TString Variable = "patType1CorrectedPFMet/MET_";


void allCutFlow(){
setTDRStyle();

//loop over variables
for(int i = 0; i<1; i++){
double MinX = MinXs[i];
double MaxX = MaxXs[i];
Variable = Variables[i];
TString Xtitle = XTitles[i];

//Data
TH1D* data = getSample("SingleMu", 1);

//MC
TH1D* tt = getSample("TTJet", lumi*225.2/6920475);

//TH1D* wjets = getSample("W1Jet", lumi*37509/57708550);
TH1D* wjets = getSample("W1Jet", lumi*5400.0/23140779);
TH1D* w2jets = getSample("W2Jets", lumi*1750.0/34041404);
TH1D* w3jets = getSample("W3Jets", lumi*519.0/15536443);
TH1D* w4jets = getSample("W4Jets", lumi*214.0/13370904);

//TH1D* zjets = getSample("DY1JetsToLL", lumi*5745.25/30457954);
TH1D* zjets = getSample("DY1JetsToLL", lumi*561.0/24042904);
TH1D* z2jets = getSample("DY2JetsToLL", lumi*181.0/21835749);
TH1D* z3jets = getSample("DY3JetsToLL", lumi*51.1/11010628);
TH1D* z4jets = getSample("DY4JetsToLL", lumi*23.04/6391785);

//TH1D* qcd = getSample("QCD_Pt-15to20_MuEnrichedPt5",     lumi*34679.3/8500505);
TH1D* qcd = getSample("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678);
TH1D* qcd2 = getSample("QCD_Pt-20to30_MuEnrichedPt5",   lumi*2.87e8 * 0.0065/8486893);
TH1D* qcd3 = getSample("QCD_Pt-30to50_MuEnrichedPt5",   lumi*6.609e7 * 0.0122/8928999);
TH1D* qcd4 = getSample("QCD_Pt-50to80_MuEnrichedPt5",   lumi*8082000.0 * 0.0218/7256011);
TH1D* qcd5 = getSample("QCD_Pt-80to120_MuEnrichedPt5",  lumi*1024000.0 * 0.0395/9030624);
TH1D* qcd6 = getSample("QCD_Pt-120to170_MuEnrichedPt5", lumi*157800.0 * 0.0473/8500505);
TH1D* qcd7 = getSample("QCD_Pt-170to300_MuEnrichedPt5", lumi*34020.0 * 0.0676/7662483);
TH1D* qcd8 = getSample("QCD_Pt-300to470_MuEnrichedPt5", lumi*1757.0 * 0.0864/7797481);
TH1D* qcd9 = getSample("QCD_Pt-470to600_MuEnrichedPt5", lumi*115.2 * 0.1024/2995767);
TH1D* qcd10 = getSample("QCD_Pt-800to1000_MuEnrichedPt5",lumi*3.57 * 0.1033/4047142);
TH1D* qcd11 = getSample("QCD_Pt-1000_MuEnrichedPt5",     lumi*0.774 * 0.1097/3807263);

TH1D* single_t = getSample("T_t-channel", lumi*56.4/3757707);
TH1D* top_t = getSample("T_t-channel", lumi*56.4/3757707);
TH1D* top_tw = getSample("T_tW-channel", lumi*11.1/497395);
TH1D* top_s = getSample("T_s-channel", lumi*3.79/249516);
TH1D* tbar_t = getSample("Tbar_t-channel", lumi*30.7/1934817);
TH1D* tbar_tw = getSample("Tbar_tW-channel", lumi*11.1/493239);
TH1D* tbar_s = getSample("Tbar_s-channel", lumi*1.76/139948);

THStack *hs = new THStack("hs","test");

  qcd->Add(qcd2);
  qcd->Add(qcd3);
  qcd->Add(qcd4);
  qcd->Add(qcd5);
  qcd->Add(qcd6);
  qcd->Add(qcd7);
  qcd->Add(qcd8);
  qcd->Add(qcd9);
  qcd->Add(qcd10);
  qcd->Add(qcd11);

  hs->Add(qcd);  

  zjets->Add(z2jets);
  zjets->Add(z3jets);
  zjets->Add(z4jets);  
    
  hs->Add(zjets);


  wjets->Add(w2jets);
  wjets->Add(w3jets);
  wjets->Add(w4jets);  

  hs->Add(wjets);
      
  hs->Add(top_t);
  hs->Add(top_tw);
  hs->Add(top_s);
  hs->Add(tbar_t);
  hs->Add(tbar_tw);
  hs->Add(tbar_s);
  
  single_t->Add(top_tw);
  single_t->Add(top_s);
  single_t->Add(tbar_t);
  single_t->Add(tbar_tw);
  single_t->Add(tbar_s);
  
  hs->Add(tt);

  //draw histos to files
  TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);

 TString step[11] = {"Skim" ,"Cleaning and HLT","one isolated #mu", "loose #mu veto", "loose e veto", "#geq 1 jets", "#geq 2 jets","#geq 3 jets", "#geq 4 jets", "#geq1 CSV b-tag", "#geq2 CSV b-tag" };
 		
  hs->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.3);
  
  if(logPlot ==true)
  hs->SetMinimum(100.);

  hs->Draw();
  for(int i =0; i<tt->GetNbinsX(); i++)
  hs->GetXaxis()->SetBinLabel(i+1, step[i]);
  
  data->Draw("E same");
  data->SetMarkerStyle(20);
  
//  hs->GetXaxis()->SetLimits(MinX, MaxX);
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
  
  TString plotName("plots/cutFlow/");
  
  if(logPlot ==true){
    plotName += Variable+"_Log";
    plotName += Nbtags+".png";
    
  }else{
    plotName += Variable+"";  
    plotName += Nbtags+".png";
  }
 
 
  c1->SaveAs(plotName);
  delete c1;
 
 //print out
 
 
 cout << " & ttbar & wjets & zjets & single-t & qcd & all MC & data " << endl; 
 
 for(int i = 0; i < tt->GetNbinsX(); i++){
 cout << step[i] << " & " << tt->GetBinContent(i+1) << " & " << wjets->GetBinContent(i+1) << " & " << zjets->GetBinContent(i+1) << " & " << single_t->GetBinContent(i+1) << " & " << qcd->GetBinContent(i+1) << " & " << tt->GetBinContent(i+1) +wjets->GetBinContent(i+1) +zjets->GetBinContent(i+1) +single_t->GetBinContent(i+1) +qcd->GetBinContent(i+1) << " & " << data->GetBinContent(i+1) << endl;
 
 }
 
  
  }
  	
}


TH1D* getSample(TString sample, double weight){
	TString dir = "rootFilesV2/central/";
	TFile* file = new TFile(dir + sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");
	
	TH1D* plot = (TH1D*) file->Get("EventCount/"+Variable);

        if(sample == "TTJet"){
	plot->SetFillColor(kRed+1);
        plot->SetLineColor(kRed+1);
	}else if(sample == "WJetsToLNu" || sample == "W1Jet" || sample == "W2Jets"|| sample == "W3Jets"|| sample == "W4Jets"){
	plot->SetLineColor(kGreen-3);	  
  	plot->SetFillColor(kGreen-3);
	}else if(sample == "DYJetsToLL" || sample == "DY1JetsToLL" || sample == "DY2JetsToLL" || sample == "DY3JetsToLL" || sample == "DY4JetsToLL"){
	plot->SetFillColor(kAzure-2);
	plot->SetLineColor(kAzure-2);
	}else if(sample == "QCD_Pt_20_MuEnrichedPt_15" || sample == "QCD_Pt-15to20_MuEnrichedPt5" || sample=="QCD_Pt-15to20_MuEnrichedPt5" || sample =="QCD_Pt-20to30_MuEnrichedPt5" || sample ==    "QCD_Pt-30to50_MuEnrichedPt5" || sample ==    "QCD_Pt-50to80_MuEnrichedPt5" || sample ==    "QCD_Pt-80to120_MuEnrichedPt5" || sample ==   "QCD_Pt-120to170_MuEnrichedPt5" || sample ==  "QCD_Pt-170to300_MuEnrichedPt5" || sample ==  "QCD_Pt-300to470_MuEnrichedPt5" || sample ==  "QCD_Pt-470to600_MuEnrichedPt5" || sample ==  "QCD_Pt-800to1000_MuEnrichedPt5" || sample =="QCD_Pt-1000_MuEnrichedPt5" 	 ){
	plot->SetFillColor(kYellow);
	plot->SetLineColor(kYellow);
	}else if(sample == "T_t-channel" || sample == "T_tW-channel" || sample == "T_s-channel" || sample == "Tbar_t-channel" || sample == "Tbar_tW-channel" || sample == "Tbar_s-channel"){
	plot->SetFillColor(kMagenta);
	plot->SetLineColor(kMagenta);
	}
	
	//plot->Scale(weight);
	plot->Rebin(rebinFact);
	
	return plot;

}

TText* doPrelim(float x, float y)
{
  std::ostringstream stream;
  stream  << "#mu, #geq 4 jets, #geq 2 b-tags               CMS Preliminary, L = 5.8 fb^{-1}";   

  TLatex* text = new TLatex(x, y, stream.str().c_str());
  //text->SetTextAlign(33);  //left
  //text->SetTextAlign(22);  //center
  //text->SetTextAlign(11);  //right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);  // for thesis

  return text;
}
