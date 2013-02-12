#include <TMinuit.h>
#include <TMath.h>
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
#include "TROOT.h"

double etaFit(TString bin, TString valReturn, TString dir);
TH1D* getSample(TString sample, double weight, int rebinFact, TString Obj, TString dir);
TText* doPrelim(float x, float y);
TH1D* getQCD(int rebinFact);
void fcn(int& npar, double* deriv, double& f, double par[], int flag);

using namespace std;

void runDiff();

double lumi = 5814;

//global histos for fit
TH1D* data; 
TH1D* top_fit;
TH1D* wjets_fit;
TH1D* zjets_fit;
TH1D* qcd_fit;
TH1D* bg_fit;

//these only need to be global is using constraints
double Nwjets, Nzjets, NQCD;

//currently set up to run first bin.. can change what to read in the read file bit
void runDiff(){
TString dir = "central";
etaFit("Binned_MET_Analysis/RecoMET_bin_0-25", "measured", dir);
}

double etaFit(TString bin, TString valReturn, TString dir){
setTDRStyle();

bool savePlots = true;


//these only need to be global is using constraints
double Nwjets, Nzjets, NQCD;

bool inclZ = false;
bool inclW = false;

if(dir == "Scale_up" || dir == "Scale_down" || dir == "Match_up" || dir == "Match_down"){
inclZ = true;
inclW = true;
}

//choose object
TString Obj = bin;
int rebinFact = 1;

double MinX = 0.;
double MaxX = 2.6;

TString Xtitle = "#left|#eta#right|_{#mu}";

//Data
data = getSample("SingleMu", 1, rebinFact, Obj, "central");

//MC
TH1D* tt = getSample("TTJet", lumi*225.2/6920475, rebinFact, Obj, dir);
TH1D* tt_tot = getSample("TTJet", lumi*225.2/6920475, rebinFact, "Muon", dir);

TH1D* wjets;
//TH1D* w1jets = getSample("W1Jet", lumi*5400.0/23140779, rebinFact, Obj, dir);
TH1D* w2jets = getSample("W2Jets", lumi*1750.0/34041404, rebinFact, Obj, dir);
TH1D* w3jets = getSample("W3Jets", lumi*519.0/15536443, rebinFact, Obj, dir);
TH1D* w4jets = getSample("W4Jets", lumi*214.0/13370904, rebinFact, Obj, dir);

TH1D* zjets;
//TH1D* z1jets = getSample("DY1JetsToLL", lumi*561.0/24042904, rebinFact, Obj, dir);
TH1D* z2jets = getSample("DY2JetsToLL", lumi*181.0/21835749, rebinFact, Obj, dir);
TH1D* z3jets = getSample("DY3JetsToLL", lumi*51.1/11010628, rebinFact, Obj, dir);
TH1D* z4jets = getSample("DY4JetsToLL", lumi*23.04/6391785, rebinFact, Obj, dir);

TH1D* qcd = getQCD(rebinFact);

TH1D* qcd_mc = getSample("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678, rebinFact, Obj, dir);
TH1D* qcd2 = getSample("QCD_Pt-20to30_MuEnrichedPt5",   lumi*2.87e8 * 0.0065/8486893, rebinFact, Obj, dir);
TH1D* qcd3 = getSample("QCD_Pt-30to50_MuEnrichedPt5",   lumi*6.609e7 * 0.0122/8928999, rebinFact, Obj, dir);
TH1D* qcd4 = getSample("QCD_Pt-50to80_MuEnrichedPt5",   lumi*8082000.0 * 0.0218/7256011, rebinFact, Obj, dir);
TH1D* qcd5 = getSample("QCD_Pt-80to120_MuEnrichedPt5",  lumi*1024000.0 * 0.0395/9030624, rebinFact, Obj, dir);
TH1D* qcd6 = getSample("QCD_Pt-120to170_MuEnrichedPt5", lumi*157800.0 * 0.0473/8500505, rebinFact, Obj, dir);
TH1D* qcd7 = getSample("QCD_Pt-170to300_MuEnrichedPt5", lumi*34020.0 * 0.0676/7662483, rebinFact, Obj, dir);
TH1D* qcd8 = getSample("QCD_Pt-300to470_MuEnrichedPt5", lumi*1757.0 * 0.0864/7797481, rebinFact, Obj, dir);
TH1D* qcd9 = getSample("QCD_Pt-470to600_MuEnrichedPt5", lumi*115.2 * 0.1024/2995767, rebinFact, Obj, dir);
TH1D* qcd10 = getSample("QCD_Pt-800to1000_MuEnrichedPt5",lumi*3.57 * 0.1033/4047142, rebinFact, Obj, dir);
TH1D* qcd11 = getSample("QCD_Pt-1000_MuEnrichedPt5",     lumi*0.774 * 0.1097/3807263, rebinFact, Obj, dir);

qcd_mc->Add(qcd2);
qcd_mc->Add(qcd3);
qcd_mc->Add(qcd4);
qcd_mc->Add(qcd5);
qcd_mc->Add(qcd6);
qcd_mc->Add(qcd7);
qcd_mc->Add(qcd8);
qcd_mc->Add(qcd9);
qcd_mc->Add(qcd10);
qcd_mc->Add(qcd11);

qcd->Scale(qcd_mc->Integral());
cout << "NQCD: " << qcd_mc->Integral() << endl;

TH1D* top_t = getSample("T_t-channel", lumi*56.4/3757707, rebinFact, Obj, dir);
TH1D* top_tw = getSample("T_tW-channel", lumi*11.1/497395, rebinFact, Obj, dir);
TH1D* top_s = getSample("T_s-channel", lumi*3.79/249516, rebinFact, Obj, dir);
TH1D* tbar_t = getSample("Tbar_t-channel", lumi*30.7/1934817, rebinFact, Obj, dir);
TH1D* tbar_tw = getSample("Tbar_tW-channel", lumi*11.1/493239, rebinFact, Obj, dir);
TH1D* tbar_s = getSample("Tbar_s-channel", lumi*1.76/139948, rebinFact, Obj, dir);

//make combined top and single top template
TH1D* top = (TH1D*)tt->Clone("top");
top->Add(top_t); top->Add(top_tw);top->Add(top_s); top->Add(tbar_t); top->Add(tbar_tw);top->Add(tbar_s);

//sum single top into one
TH1D* single_top = (TH1D*)top_t->Clone("single top");
single_top->Add(top_tw);single_top->Add(top_s); single_top->Add(tbar_t); single_top->Add(tbar_tw);single_top->Add(tbar_s);

  
THStack *hs = new THStack("hs","test");

  hs->Add(qcd);
    
  if(inclZ == true){
  zjets = getSample("DYJetsToLL", lumi*5745.25/30457954, rebinFact, Obj, dir);

  }else{
  zjets  = getSample("DY1JetsToLL", lumi*561.0/24042904, rebinFact, Obj, dir);
  zjets->Add(z2jets);
  zjets->Add(z3jets);
  zjets->Add(z4jets);  
  }
  
  if(inclW == true){
  wjets = getSample("WJetsToLNu", lumi*37509/57708550, rebinFact, Obj, dir);
  }else{
  wjets = getSample("W1Jet", lumi*5400.0/23140779, rebinFact, Obj, dir); 
  wjets->Add(w2jets);
  wjets->Add(w3jets);
  wjets->Add(w4jets);
  }
  
  hs->Add(zjets);
  hs->Add(wjets);
      
  hs->Add(top_t);
  hs->Add(top_tw);
  hs->Add(top_s);
  hs->Add(tbar_t);
  hs->Add(tbar_tw);
  hs->Add(tbar_s);
  
  hs->Add(tt);

//combined histo for pseudo?
TH1D* allMC = (TH1D*)top->Clone("allMC");
allMC->Add(wjets); allMC->Add(zjets); allMC->Add(qcd);

if(savePlots ==true){

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

  
  TString plotName("plots/Control/Muon/");
  
    plotName += "absEta";  
    plotName += "_ge2btags.pdf";
  
//  c1->SaveAs(plotName);
  delete c1;
}

//clone and scale
top_fit = (TH1D*)top->Clone("top fit");
wjets_fit = (TH1D*)wjets->Clone("wjets fit");
zjets_fit = (TH1D*)zjets->Clone("zjets fit");
qcd_fit = (TH1D*)qcd->Clone("qcd fit");

bg_fit = (TH1D*)wjets_fit->Clone("bg fit");
bg_fit->Add(zjets_fit);
bg_fit->Add(qcd_fit);

top_fit->Scale(1./ top_fit->Integral());
wjets_fit->Scale(1./ wjets_fit->Integral()); 
zjets_fit->Scale(1./ zjets_fit->Integral()); 
qcd_fit->Scale(1./ qcd_fit->Integral());
bg_fit->Scale(1./ bg_fit->Integral());
  
  
  if(savePlots == true){
  //draw histos to files
  TCanvas *c2 = new TCanvas("Plot","Plot",900, 600);
  
  top_fit->SetFillColor(kWhite); wjets_fit->SetFillColor(kWhite); zjets_fit->SetFillColor(kWhite); qcd_fit->SetFillColor(kWhite); bg_fit->SetFillColor(kWhite); bg_fit->SetLineColor(kBlack);
  top_fit->Draw();
  bg_fit->Draw("same");
  wjets_fit->Draw("same");
  zjets_fit->Draw("same");
  qcd_fit->Draw("same");  
  
  top_fit->SetAxisRange(MinX, MaxX);
  top_fit->GetXaxis()->SetTitle(Xtitle); top_fit->GetXaxis()->SetTitleSize(0.05);
  top_fit->GetYaxis()->SetTitle("Normalised Events");top_fit->GetYaxis()->SetTitleSize(0.05);
    
  	TLegend *tleg3;
	tleg3 = new TLegend(0.65,0.7,0.8,0.9);
	tleg3->SetTextSize(0.04);
	tleg3->SetBorderSize(0);
	tleg3->SetFillColor(10);
	
	tleg3->AddEntry(top_fit , "signal", "l");
	tleg3->AddEntry(bg_fit , "background", "l");
	tleg3->AddEntry(wjets_fit , "w+jets", "l");
	tleg3->AddEntry(zjets_fit , "z+jets", "l");
	tleg3->AddEntry(qcd_fit , "QCD", "l");
 	tleg3->Draw("same");	
	
	TText* textPrelim2 = doPrelim(0.17,0.96);
	textPrelim2->Draw();
  c2->SaveAs("plots/Fits/"+bin+"_Template.pdf");
  delete c2;
 }
 
 
 
int Ntotal = data->Integral();
double Nsignal = top->Integral();
Nwjets = wjets->Integral();
Nzjets = zjets->Integral();
NQCD = qcd->Integral();

  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);

  minuit.SetPrintLevel(-1);
  minuit.SetErrorDef(1.);
  
  
  int ierflg = 0;
  double Nbg=  wjets->Integral()+zjets->Integral()+qcd->Integral();
  string parName[npar] = {"ttbar+single-top", "background"}; //background parameters
  double par[npar] = {top->Integral(), Nbg};               //using the MC estimation as the start values 1fb
  
  cout << "total data events: " << Ntotal << endl;

  for(int i=0; i<npar; i++){

    //minuit.mnparm(i, parName[i], par[i],10., -1.e6, 1.e6, ierflg);
    minuit.mnparm(i, parName[i], par[i], 10., 0, Ntotal, ierflg);

  }
  
   //the following is copied from Fabian's fitting code to improve minimum, but you can comment it, it won't affect the fitting results.
  // 1 standard
  // 2 try to improve minimum (slower)
  double arglist[10];
  arglist[0]=2;
  minuit.mnexcm("SET STR",arglist,1,ierflg);
  minuit.Migrad();
  
  double outpar[npar], err[npar];
    
  for (int i=0; i<npar; i++){
    minuit.GetParameter(i,outpar[i],err[i]);
  }

  for (int i=0; i<top_fit->GetNbinsX() ; i++){
  
 	top_fit->SetBinContent(i+1, top_fit->GetBinContent(i+1)*outpar[0]);
	bg_fit->SetBinContent(i+1, wjets_fit->GetBinContent(i+1)*outpar[1]);
	//zjets_fit->SetBinContent(i+1, zjets_fit->GetBinContent(i+1)*outpar[2]);
	//qcd_fit->SetBinContent(i+1, qcd_fit->GetBinContent(i+1)*outpar[3]);
  
  }
  
  //print out the results for all templates
//   cout <<" \n Total number of events after the fit" << endl;
//   cout<<"   & ttbar+single top & w+jets & z+jets & qcd "<<endl;
//   cout <<  " & " << Nsignal <<  " & " << Nwjets << " & " <<  Nzjets << " & " <<  NQCD  <<endl;
//   cout<< " & "<<outpar[0]<<"+-"<<err[0]<<" & "<<outpar[1]<<"+-"<<err[1]<<" & "<<outpar[2]<<"+-"<<err[2]<<" & "<<outpar[3]<<"+-"<<err[3]<<endl;

  cout <<" \n Total number of events after the fit" << endl;
  cout<<"   & ttbar+single top & w+jets & z+jets & qcd "<<endl;
  cout <<  " & " << Nsignal <<  " & " << Nwjets + Nzjets + NQCD  <<endl;
  cout<< " & "<<outpar[0] << "+-" <<err[0] << " & " <<outpar[1]<<"+-"<<err[1] <<endl; 

  if(savePlots == true){
    TCanvas *c3 = new TCanvas("Plot","Plot",900, 600);
  
  THStack* sum_fit = new THStack("sum fit","stacked histograms"); //used for stack plot
  qcd_fit->SetFillColor(kYellow); zjets_fit->SetFillColor(kBlue);   wjets_fit->SetFillColor(kGreen);  top_fit->SetFillColor(kRed);
  //sum_fit->Add(qcd_fit); sum_fit->Add(zjets_fit);  sum_fit->Add(wjets_fit);  
  top_fit->SetLineColor(kBlack);
  
  bg_fit->SetFillColor(kGreen);
  sum_fit->Add(bg_fit);sum_fit->Add(top_fit);
 
  sum_fit->Draw();
  data->Draw("E same");
  
  sum_fit->GetXaxis()->SetLimits(MinX, MaxX);
  sum_fit->GetXaxis()->SetTitle(Xtitle); sum_fit->GetXaxis()->SetTitleSize(0.05);
  sum_fit->GetYaxis()->SetTitle("Number of Events");sum_fit->GetYaxis()->SetTitleSize(0.05);
    
  	TLegend *tleg4;
	tleg4 = new TLegend(0.65,0.7,0.8,0.9);
	tleg4->SetTextSize(0.04);
	tleg4->SetBorderSize(0);
	tleg4->SetFillColor(10);
	
	tleg4->AddEntry(top_fit , "signal", "ef");
	tleg4->AddEntry(bg_fit , "background", "ef");
 	tleg4->Draw("same");	
	
	TText* textPrelim3 = doPrelim(0.17,0.96);
	textPrelim3->Draw();
    c3->SaveAs("plots/Fits/"+bin+"_Fit.pdf");
    delete c3;

 }

//cout << "cross section is:  " <<  ((outpar[0]-single_top->Integral())/ tt_tot->Integral())*225.2 << endl;

if(valReturn == "measured"){
return outpar[0];
}else if(valReturn == "measuredErr"){
return err[0];
}else if(valReturn == "bgscale"){
return (outpar[1]/Nbg);

}

return 0;

}


TH1D* getSample(TString sample, double weight, int rebinFact, TString Obj, TString dir){
		
	TFile* file = new TFile();
	
	if(dir == "central" || (dir == "Scale_up_tt" && sample != "TTJet") || (dir == "Scale_down_tt" && sample != "TTJet") || ((dir == "Scale_up" && sample != "WJetsToLNu") && (dir == "Scale_up"  && sample != "DYJetsToLL"))  || ((dir == "Scale_down" && sample != "WJetsToLNu") && (dir == "Scale_down"  && sample != "DYJetsToLL")) || (dir == "Match_up_tt" && sample != "TTJet") || (dir == "Match_down_tt" && sample != "TTJet") || ((dir == "Match_up" && sample != "WJetsToLNu") && (dir == "Match_up"  && sample != "DYJetsToLL"))  || ((dir == "Match_down" && sample != "WJetsToLNu") && (dir == "Match_down"  && sample != "DYJetsToLL")) || dir == "UnclusteredEnUp" || dir == "UnclusteredEnDown" || dir == "JetEnUp" || dir == "JetEnDown" || dir == "JetResUp" || dir == "JetResDown" || dir == "TauEnUp" || dir == "TauEnDown" || dir == "MuonEnUp" || dir == "MuonEnDown" || dir == "ElectronEnUp" || dir == "ElectronEnDown"){
	file = new TFile("rootFilesV2/central/"+ sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	}else if(dir == "JES_up")
	file = new TFile("rootFilesV2/"+ dir +"/"+ sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_plusJES.root");
	else if(dir == "JES_down")
	file = new TFile("rootFilesV2/"+ dir +"/"+ sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_minusJES.root");
	else if(dir == "BJet_up")
	file = new TFile("rootFilesV2/"+ dir +"/"+ sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_plusBjet.root");
	else if(dir == "BJet_down")
	file = new TFile("rootFilesV2/"+ dir +"/"+ sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_minusBJet.root");
	else if(dir == "PU_up")
	file = new TFile("rootFilesV2/"+ dir +"/"+ sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_PU_72765mb.root");
	else if(dir == "PU_down")
	file = new TFile("rootFilesV2/"+ dir +"/"+ sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_PU_65835mb.root");
	else if(dir == "Scale_up_tt" && sample == "TTJet")
	file = new TFile("rootFilesV2/"+ dir +"/TTJets-scaleup_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Scale_down_tt" && sample == "TTJet")
	file = new TFile("rootFilesV2/"+ dir +"/TTJets-scaledown_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Scale_up" && sample == "WJetsToLNu")
	file = new TFile("rootFilesV2/"+ dir +"/WJets-scaleup_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Scale_up" && sample == "DYJetsToLL")
	file = new TFile("rootFilesV2/"+ dir +"/ZJets-scaleup_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Scale_down" && sample == "WJetsToLNu")
	file = new TFile("rootFilesV2/"+ dir +"/WJets-scaledown_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Scale_down" && sample == "DYJetsToLL")
	file = new TFile("rootFilesV2/"+ dir +"/ZJets-scaledown_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Match_up_tt" && sample == "TTJet")
	file = new TFile("rootFilesV2/"+ dir +"/TTJets-matchingup_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Match_down_tt" && sample == "TTJet")
	file = new TFile("rootFilesV2/"+ dir +"/TTJets-matchingdown_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Match_up" && sample == "WJetsToLNu")
	file = new TFile("rootFilesV2/"+ dir +"/WJets-matchingup_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Match_up" && sample == "DYJetsToLL")
	file = new TFile("rootFilesV2/"+ dir +"/ZJets-matchingup_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Match_down" && sample == "WJetsToLNu")
	file = new TFile("rootFilesV2/"+ dir +"/WJets-matchingdown_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(dir == "Match_down" && sample == "DYJetsToLL")
	file = new TFile("rootFilesV2/"+ dir +"/ZJets-matchingdown_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	else if(sample == "TTJet_MCNLO")
	file = new TFile("rootFilesV2/central/TTJet_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_MCatNLO.root");
	else if(sample == "TTJet_POWHEG")
	file = new TFile("rootFilesV2/central/TTJet_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET_POWHEG.root");
	
	TH1D* plot; 
	TH1D* plot2;
	TH1D* plot3;
	
	
	
	if(Obj=="Muon"){
	plot = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+"/muon_AbsEta_2btags");
	plot2 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+"/muon_AbsEta_3btags");
	plot3 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+"/muon_AbsEta_4orMoreBtags");
	}else{
	plot = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+"/muon_absolute_eta_2btags");
	plot2 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+"/muon_absolute_eta_3btags");
	plot3 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/"+Obj+"/muon_absolute_eta_4orMoreBtags");
	}


	plot->Add(plot2);
	plot->Add(plot3);

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
	}else if(sample == "SingleMu")
	plot->SetLineColor(kBlack);
	
	//plot->Scale(weight);
	plot->Rebin(rebinFact);
	
	plot->SetDirectory(gROOT);
	
	file->Close();
	return plot;

}

TH1D* getQCD(int rebinFact){
	TString dir = "rootFiles/";
	
	TFile* file = new TFile(dir +"qcdest.root");		
	TH1D* plot = (TH1D*) file->Get("muon_AbsEta_0btag");

// 	for(int i = 1; i <= plot->GetNbinsX(); i++){
// 	plot->SetBinError(i, 0.0);
// 	}

	plot->SetFillColor(kYellow);
	plot->SetLineColor(kYellow);
	plot->SetMarkerStyle(1);
		
	TH1D* copyplot = new TH1D("qcd plot", "qcd plot", 30, 0.0, 3.0);
	
	for(int i = 1; i <= plot->GetNbinsX(); i++){
	copyplot->SetBinContent(i, plot->GetBinContent(i));
	//copyplot->SetBinError(i, plot->GetBinError(i));
	}
	copyplot->SetFillColor(kYellow);
	copyplot->SetLineColor(kYellow);
	copyplot->SetMarkerStyle(1);
	copyplot->Scale(1./copyplot->Integral());	
	copyplot->Rebin(rebinFact);
	
	//file->Close("R");
	return copyplot;
	//file->Close();
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

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

   double lnL = 0.0;


  for (int i=0; i< top_fit->GetNbinsX(); i++){

    //data_i is the observed number of events in each bin
    int data_i = data->GetBinContent(i+1);
    //xi is the expected number of events in each bin
    double xi = par[0]*top_fit->GetBinContent(i+1) + par[1]*bg_fit->GetBinContent(i+1);

   
    if(data_i !=0 && xi != 0){
      lnL += log(TMath::Poisson(data_i, xi));
    }
    
  }

  //W+jets, Z+jets constraints
  f = -2.0 * lnL;
  
// 
//   double nwjets = Nwjets;
// //  double nwjets_err = nwjets*0.3;
//   //double nwjets_err = nwjets*0.02;
//    
//   double nzjets = Nzjets;
// //  double nzjets_err = nzjets*0.1;
// 
// double nqcd = 0;
// if(NQCD>0){
//   nqcd = NQCD;
//   }
//   else{
//   nqcd = 0.00000001;
//   }
//   double nqcd_err = nqcd*1.;

  //ratio constraints
   //f += ( (par[2]/par[1] - nzjets/nwjets) / (0.05 *nzjets/nwjets) )  * ( (par[2]/par[1] - nzjets/nwjets) / (0.05*nzjets/nwjets) ); 


   //f += ((par[3]-nqcd)*(par[3]-nqcd))/nqcd_err/nqcd_err;



}                         
