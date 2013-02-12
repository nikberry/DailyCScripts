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

void doPlotsQCD();
TH1D* getSample(TString sample, double weight);
TH1D* getCentral(TString sample, double weight);
TH1D* getQCD(double weight);
TText* doPrelim(float x, float y);

double lumi = 5800;
//stuff to choose
bool logPlot = false; //true for log plot
bool savePlots = false; //appart from the last plot
int rebinFact = 10;

//isolation selection
//TString Isolation = "QCD No Iso/";
//TString Isolation = "Ref selection/";
//TString Isolation = "QCD mu+jets PFRelIso/";
TString Isolation = "QCD non iso mu+jets/";

// number of btags
TString Nbtags = "0btag";  //standard  "2btags" , qcd "0btag"

bool inclZ = false;
bool inclW = false;
bool inclQ = false;

//choose object
TString Obj = "Muon/";

//muon variables
const int N = 6;
TString Variable;
TString Variables[N] = {"muon_AbsEta_", "muon_eta_", "muon_pfIsolation_04_", "muon_pT_", "muon_phi_", "muon_dB_"};
double MinXs[N] = {0,-2.6 ,0 , 0, -3.4, 0};
double MaxXs[N] = {2.6,2.6 ,0.5 , 400,3.2 , 0.3};
TString XTitles[N] = {"#left|#eta#right|_{#mu}", "#eta_{#mu}", "RelIso_{#mu}", "p_{T}(#mu)", "#phi_{#mu}", "dB_{#mu}"};


void doPlotsQCD(){
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

TH1D* wjets;
TH1D* w1jets = getSample("W1Jet", lumi*5400.0/23140779);
TH1D* w2jets = getSample("W2Jets", lumi*1750.0/34041404);
TH1D* w3jets = getSample("W3Jets", lumi*519.0/15536443);
TH1D* w4jets = getSample("W4Jets", lumi*214.0/13370904);

if(inclW ==true){
wjets = getSample("WJetsToLNu", lumi*37509/57708550);
}else{
wjets  = getSample("W1Jet", lumi*5400.0/23140779);
wjets->Add(w2jets);
wjets->Add(w3jets);
wjets->Add(w4jets);
}

TH1D* zjets;
TH1D* z1jets = getSample("DY1JetsToLL", lumi*561.0/24042904);
TH1D* z2jets = getSample("DY2JetsToLL", lumi*181.0/21835749);
TH1D* z3jets = getSample("DY3JetsToLL", lumi*51.1/11010628);
TH1D* z4jets = getSample("DY4JetsToLL", lumi*23.04/6391785);

if(inclZ ==true){
zjets = getSample("DYJetsToLL", lumi*5745.25/30457954);
}else{
zjets  = getSample("DY1JetsToLL", lumi*561.0/24042904);
zjets->Add(z2jets);
zjets->Add(z3jets);
zjets->Add(z4jets);
}

TH1D* qcd;

TH1D* qcd1 = getSample("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678);
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

if(inclQ ==true){
qcd = getSample("QCD_Pt_20_MuEnrichedPt_15", lumi*34679.3/8500505);
}else{
qcd  = getSample("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678);
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
}



TH1D* top_t = getSample("T_t-channel", lumi*56.4/3757707);
TH1D* top_tw = getSample("T_tW-channel", lumi*11.1/497395);
TH1D* top_s = getSample("T_s-channel", lumi*3.79/249516);
TH1D* tbar_t = getSample("Tbar_t-channel", lumi*30.7/1934817);
TH1D* tbar_tw = getSample("Tbar_tW-channel", lumi*11.1/493239);
TH1D* tbar_s = getSample("Tbar_s-channel", lumi*1.76/139948);

THStack *hs = new THStack("hs","test");
THStack *qcdstack = new THStack("hs","test");

  hs->Add(qcd);
  qcdstack->Add(qcd); 
  hs->Add(zjets);

  hs->Add(wjets);
      
  hs->Add(top_t);
  hs->Add(top_tw);
  hs->Add(top_s);
  hs->Add(tbar_t);
  hs->Add(tbar_tw);
  hs->Add(tbar_s);
  
  hs->Add(tt);

  //draw histos to files
  TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);
		
  hs->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.5);

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
  
  TString plotName("plots/Control/QCD/");
  
  if(logPlot ==true){
    plotName += Variable+"_Log";
    plotName += Nbtags+".pdf";
    
  }else{
    plotName += Variable;  
    plotName += Nbtags+".pdf";
  }
 
if(savePlots == true){ 
  c1->SaveAs(plotName);
  }
  delete c1;

 TH1D*dataclone = (TH1D*)data->Clone("dataclone");
 TH1D*mc_sub = (TH1D*)tt->Clone("mc_clone");
 
  mc_sub->Add(top_t); mc_sub->Add(top_tw); mc_sub->Add(top_s); mc_sub->Add(tbar_t); mc_sub->Add(tbar_tw); mc_sub->Add(tbar_s); mc_sub->Add(z1jets); mc_sub->Add(z2jets); mc_sub->Add(z3jets); mc_sub->Add(z4jets);  mc_sub->Add(w1jets); mc_sub->Add(w2jets); mc_sub->Add(w3jets); mc_sub->Add(w4jets); 
 
  //subtract samples
  data->Add(mc_sub, -1);

for(int i = 0; i<data->GetNbinsX(); i++){
data->SetBinError(i+1, sqrt(pow(dataclone->GetBinError(i+1),2)+pow(0.5*mc_sub->GetBinContent(i+1),2)));

}   
    
    //draw histos to files
  TCanvas *c2 = new TCanvas("Plot","Plot",900, 600);

  qcdstack->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.5);

  qcdstack->Draw();
  data->Draw("E same");
  data->SetMarkerStyle(20);		
    
     qcdstack->GetXaxis()->SetLimits(MinX, MaxX);

     qcdstack->GetXaxis()->SetTitle(Xtitle); qcdstack->GetXaxis()->SetTitleSize(0.05);
     qcdstack->GetYaxis()->SetTitle("Number of Events"); qcdstack->GetYaxis()->SetTitleSize(0.05); 
	
     TLegend *tleg3;
     tleg3 = new TLegend(0.7,0.7,0.8,0.9);
     tleg3->SetTextSize(0.04);
     tleg3->SetBorderSize(0);
     tleg3->SetFillColor(10);
     tleg3->AddEntry(data , "2012 data", "lpe");
     tleg3->AddEntry(qcd , "QCD", "lf");

     tleg2->Draw("same");    

     TText* textPrelim2 = doPrelim(0.17,0.96);
     textPrelim2->Draw();
      
     TString plotName2("plots/Control/QCD/");
  
  if(logPlot ==true){
    plotName2 += Variable+"Test_Log";
    plotName2 += Nbtags+".pdf";
    
  }else{
    plotName2 += Variable+"Subtract";  
    plotName2 += Nbtags+".pdf";
  }

if(savePlots == true){ 
  c2->SaveAs(plotName2);
  }
  delete c2;

 
 data->Scale(1./data->Integral());

TFile qcdfile("qcdest.root", "RECREATE", "comment");

  data->Write();
  qcdfile.Write();
  qcdfile.Close();
    
  //corrections
TH1D* central;

TH1D* central1 = getCentral("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678);
TH1D* central2 = getCentral("QCD_Pt-20to30_MuEnrichedPt5",   lumi*2.87e8 * 0.0065/8486893);
TH1D* central3 = getCentral("QCD_Pt-30to50_MuEnrichedPt5",   lumi*6.609e7 * 0.0122/8928999);
TH1D* central4 = getCentral("QCD_Pt-50to80_MuEnrichedPt5",   lumi*8082000.0 * 0.0218/7256011);
TH1D* central5 = getCentral("QCD_Pt-80to120_MuEnrichedPt5",  lumi*1024000.0 * 0.0395/9030624);
TH1D* central6 = getCentral("QCD_Pt-120to170_MuEnrichedPt5", lumi*157800.0 * 0.0473/8500505);
TH1D* central7 = getCentral("QCD_Pt-170to300_MuEnrichedPt5", lumi*34020.0 * 0.0676/7662483);
TH1D* central8 = getCentral("QCD_Pt-300to470_MuEnrichedPt5", lumi*1757.0 * 0.0864/7797481);
TH1D* central9 = getCentral("QCD_Pt-470to600_MuEnrichedPt5", lumi*115.2 * 0.1024/2995767);
TH1D* central10 = getCentral("QCD_Pt-800to1000_MuEnrichedPt5",lumi*3.57 * 0.1033/4047142);
TH1D* central11 = getCentral("QCD_Pt-1000_MuEnrichedPt5",     lumi*0.774 * 0.1097/3807263);

if(inclQ ==true){
central = getCentral("QCD_Pt_20_MuEnrichedPt_15", lumi*34679.3/8500505);
}else{
central  = getCentral("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678);
central->Add(central2);
central->Add(central3);
central->Add(central4);
central->Add(central5);
central->Add(central6);
central->Add(central7);
central->Add(central8);
central->Add(central9);
central->Add(central10);
central->Add(central11);
}

double qcdTot = central->Integral();
TH1D*qcd_mc = (TH1D*)central->Clone("qcd_mc_clone");

central->Scale(1./central->Integral());
qcd->Scale(1./qcd->Integral());

  TCanvas *c3 = new TCanvas("Plot","Plot",900, 600);
	
	central->Divide(qcd);
	central->Draw("E");
	central->SetMaximum(5);
	central->SetLineColor(kBlack);	
	central->SetMarkerStyle(20);
	central->SetAxisRange(MinX, MaxX);
  	central->GetYaxis()->SetTitle("C_{F}"); central->GetYaxis()->SetTitleSize(0.05);
        central->GetXaxis()->SetTitle(Xtitle); central->GetXaxis()->SetTitleSize(0.05);
if(savePlots ==true){  
  c3->SaveAs("plots/Control/QCD/QCD_Corrections.pdf");
  }
  delete c3;
  
  
  //Data
TH1D* data_ge4 = getCentral("SingleMu", 1);

//MC
TH1D* tt_ge4 = getCentral("TTJet", lumi*225.2/6920475);

TH1D* wjets_ge4;
TH1D* w1jets_ge4 = getCentral("W1Jet", lumi*5400.0/23140779);
TH1D* w2jets_ge4 = getCentral("W2Jets", lumi*1750.0/34041404);
TH1D* w3jets_ge4 = getCentral("W3Jets", lumi*519.0/15536443);
TH1D* w4jets_ge4 = getCentral("W4Jets", lumi*214.0/13370904);

if(inclW ==true){
wjets_ge4 = getCentral("WJetsToLNu", lumi*37509/57708550);
}else{
wjets_ge4  = getCentral("W1Jet", lumi*5400.0/23140779);
wjets_ge4->Add(w2jets_ge4);
wjets_ge4->Add(w3jets_ge4);
wjets_ge4->Add(w4jets_ge4);
}

TH1D* zjets_ge4;
TH1D* z1jets_ge4 = getCentral("DY1JetsToLL", lumi*561.0/24042904);
TH1D* z2jets_ge4 = getCentral("DY2JetsToLL", lumi*181.0/21835749);
TH1D* z3jets_ge4 = getCentral("DY3JetsToLL", lumi*51.1/11010628);
TH1D* z4jets_ge4 = getCentral("DY4JetsToLL", lumi*23.04/6391785);

if(inclZ ==true){
zjets_ge4 = getCentral("DYJetsToLL", lumi*5745.25/30457954);
}else{
zjets_ge4  = getCentral("DY1JetsToLL", lumi*561.0/24042904);
zjets_ge4->Add(z2jets_ge4);
zjets_ge4->Add(z3jets_ge4);
zjets_ge4->Add(z4jets_ge4);
}

TH1D* top_t_ge4 = getCentral("T_t-channel", lumi*56.4/3757707);
TH1D* top_tw_ge4 = getCentral("T_tW-channel", lumi*11.1/497395);
TH1D* top_s_ge4 = getCentral("T_s-channel", lumi*3.79/249516);
TH1D* tbar_t_ge4 = getCentral("Tbar_t-channel", lumi*30.7/1934817);
TH1D* tbar_tw_ge4 = getCentral("Tbar_tW-channel", lumi*11.1/493239);
TH1D* tbar_s_ge4 = getCentral("Tbar_s-channel", lumi*1.76/139948);

TH1D* qcd_ge4 = getQCD(lumi*34679.3/8500505);
//TH1D* qcd_mc = getCentral("QCD_Pt_20_MuEnrichedPt_15", lumi*34679.3/8500505);
qcd_ge4->Scale(qcdTot);

THStack *hsge4 = new THStack("hs","test");
  
  hsge4->Add(qcd_mc);
  //hsge4->Add(qcd_ge4);
  hsge4->Add(zjets_ge4);

  hsge4->Add(wjets_ge4);
      
  hsge4->Add(top_t_ge4);
  hsge4->Add(top_tw_ge4);
  hsge4->Add(top_s_ge4);
  hsge4->Add(tbar_t_ge4);
  hsge4->Add(tbar_tw_ge4);
  hsge4->Add(tbar_s_ge4);
  
  hsge4->Add(tt_ge4);

  //draw histos to files
  TCanvas *c4 = new TCanvas("Plot","Plot",900, 600);
		
  hsge4->SetMaximum(data_ge4->GetBinContent(data_ge4->GetMaximumBin())*1.5);

  hsge4->Draw();
  data_ge4->Draw("E same");
  
  cout << "tot data: " << data_ge4->Integral() << endl;
  
  data_ge4->SetMarkerStyle(20);
  
  hsge4->GetXaxis()->SetLimits(MinX, MaxX);
  hsge4->GetXaxis()->SetTitle(Xtitle); hsge4->GetXaxis()->SetTitleSize(0.05);
  hsge4->GetYaxis()->SetTitle("Number of Events");hsge4->GetYaxis()->SetTitleSize(0.05);
  
  
  	TLegend *tleg4;
	tleg4 = new TLegend(0.7,0.7,0.8,0.9);
	tleg4->SetTextSize(0.04);
	tleg4->SetBorderSize(0);
	tleg4->SetFillColor(10);
	tleg4->AddEntry(data , "2012 data", "lpe");
	tleg4->AddEntry(tt , "t#bar{t}", "lf");
	tleg4->AddEntry(top_t, "single top", "lf");
	tleg4->AddEntry(wjets , "w+jets", "lf");
	tleg4->AddEntry(zjets , "z+jets", "lf");
	tleg4->AddEntry(qcd , "QCD", "lf");
	
	//tleg4->AddEntry(singtEff, "single-t"      , "l");
	//tleg4->AddEntry(singtwEff, "single-tW"      , "l");
 	tleg4->Draw("same");	
	
	TText* textPrelim4 = doPrelim(0.17,0.96);
	textPrelim4->Draw();
	
  if(logPlot ==true){
  c1->SetLogy();
  }	
  
  TString plotName4("plots/Control/QCD/");
  
  if(logPlot ==true){
    plotName4 += Variable+"Test_Log";
    plotName4 += Nbtags+".pdf";
    
  }else{
    plotName4 += Variable+"QCD_MC";  
    plotName4 += ".pdf";
  }
	
  c4->SaveAs(plotName4);

  delete c4;


  }// loop over variables
  	
}


TH1D* getSample(TString sample, double weight){
	TString dir = "rootFilesV2/central/";
	TFile* file = new TFile(dir + sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");
	
	TH1D* plot = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/"+Isolation+Obj+Variable+"0btag");

        if(sample == "TTJet"){
	plot->SetFillColor(kRed+1);
        plot->SetLineColor(kRed+1);
	}else if(sample == "WJetsToLNu" || sample == "W1Jet" || sample == "W2Jets"|| sample == "W3Jets"|| sample == "W4Jets"){
	plot->SetLineColor(kGreen-3);	  
  	plot->SetFillColor(kGreen-3);
	}else if(sample == "DYJetsToLL" || sample == "DY1JetsToLL" || sample == "DY2JetsToLL" || sample == "DY3JetsToLL" || sample == "DY4JetsToLL"){
	plot->SetFillColor(kAzure-2);
	plot->SetLineColor(kAzure-2);
	}else if(sample == "QCD_Pt_20_MuEnrichedPt_15" || sample == "QCD_Pt-15to20_MuEnrichedPt5" || sample == "QCD_Pt-20to30_MuEnrichedPt5" || sample == "QCD_Pt-30to50_MuEnrichedPt5" || sample == "QCD_Pt-50to80_MuEnrichedPt5" || sample == "QCD_Pt-80to120_MuEnrichedPt5"|| sample == "QCD_Pt-120to170_MuEnrichedPt5" || sample == "QCD_Pt-170to300_MuEnrichedPt5" || sample == "QCD_Pt-300to470_MuEnrichedPt5" || sample == "QCD_Pt-470to600_MuEnrichedPt5" || sample == "QCD_Pt-800to1000_MuEnrichedPt5" || sample == "QCD_Pt-1000_MuEnrichedPt5"){
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

TH1D* getCentral(TString sample, double weight){
	TString dir = "rootFilesV2/central/";
	TFile* file = new TFile(dir + sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");
	
	TH1D* plot = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+Variable+"0btag");
	
	TH1D* plot1 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+Variable+"1btag");
	TH1D* plot2 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+Variable+"2btags");
        TH1D* plot3 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+Variable+"3btags");
	TH1D* plot4 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+Variable+"4orMoreBtags");
	
	plot->Add(plot1);
	plot->Add(plot2);
	plot->Add(plot3);
	plot->Add(plot4);
	
        if(sample == "TTJet"){
	plot->SetFillColor(kRed+1);
        plot->SetLineColor(kRed+1);
	}else if(sample == "WJetsToLNu" || sample == "W1Jet" || sample == "W2Jets"|| sample == "W3Jets"|| sample == "W4Jets"){
	plot->SetLineColor(kGreen-3);	  
  	plot->SetFillColor(kGreen-3);
	}else if(sample == "DYJetsToLL" || sample == "DY1JetsToLL" || sample == "DY2JetsToLL" || sample == "DY3JetsToLL" || sample == "DY4JetsToLL"){
	plot->SetFillColor(kAzure-2);
	plot->SetLineColor(kAzure-2);
	}else if(sample == "QCD_Pt_20_MuEnrichedPt_15" || sample == "QCD_Pt-15to20_MuEnrichedPt5" || sample == "QCD_Pt-20to30_MuEnrichedPt5" || sample == "QCD_Pt-30to50_MuEnrichedPt5" || sample == "QCD_Pt-50to80_MuEnrichedPt5" || sample == "QCD_Pt-80to120_MuEnrichedPt5"|| sample == "QCD_Pt-120to170_MuEnrichedPt5" || sample == "QCD_Pt-170to300_MuEnrichedPt5" || sample == "QCD_Pt-300to470_MuEnrichedPt5" || sample == "QCD_Pt-470to600_MuEnrichedPt5" || sample == "QCD_Pt-800to1000_MuEnrichedPt5" || sample == "QCD_Pt-1000_MuEnrichedPt5"){
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

TH1D* getQCD(double weight){
	TString dir = "rootFiles/";
	
	TFile* file = new TFile(dir +"qcdest.root");		
	TH1D* plot = (TH1D*) file->Get("muon_AbsEta_0btag");

// 	for(int i = 1; i <= plot->GetNbinsX(); i++){
// 	plot->SetBinError(i, 0.0);
// 	}

	plot->SetFillColor(kYellow);
	plot->SetLineColor(kYellow);
	plot->SetMarkerStyle(1);
		
	//plot->Scale(weight);	
		
	TH1D* copyplot = new TH1D("qcd plot", "qcd plot", 30, 0.0, 3.0);
	
	for(int i = 1; i <= plot->GetNbinsX(); i++){
	copyplot->SetBinContent(i, plot->GetBinContent(i));
	//copyplot->SetBinError(i, plot->GetBinError(i));
	}
	copyplot->SetFillColor(kYellow);
	copyplot->SetLineColor(kYellow);
	copyplot->SetMarkerStyle(1);

	
	return copyplot;

}

TText* doPrelim(float x, float y)
{
  std::ostringstream stream;
  stream  << "#mu, #geq 4 jets                      CMS Preliminary, L = 5.8 fb^{-1}";   

  TLatex* text = new TLatex(x, y, stream.str().c_str());
  //text->SetTextAlign(33);  //left
  //text->SetTextAlign(22);  //center
  //text->SetTextAlign(11);  //right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);  // for thesis

  return text;
}
