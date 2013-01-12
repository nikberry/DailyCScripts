#include "TFile.h"
#include "TH1.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include <string.h>
#include <iostream>
#include <iomanip>
#include "tdrstyle.C"


void cutFlow();
TH1D* getSample(TString sample, bool muon);
void printCutflow(string step, int bin, TH1D* ttbar);
void weightHisto(TH1D* sample, double lumi, double xsect, double NEvts);

void cutFlow(){
setTDRStyle();
	bool muon = true;
	bool scale = true;

double lumi = 5800;

TH1D* tt = getSample("TTJet", muon);
TH1D* ttEff = new TH1D("cut eff","cut eff",10,0,10);
TH1D* wjets = getSample("WJetsToLNu", muon);
TH1D* wjetsEff = new TH1D("cut eff","cut eff",10,0,10);
TH1D* zjets = getSample("DYJetsToLL", muon);
TH1D* zjetsEff = new TH1D("cut eff","cut eff",10,0,10);
TH1D* qcd = getSample("QCD_Pt-120to170_Mu", muon);
TH1D* qcdEff = new TH1D("cut eff","cut eff",10,0,10);
TH1D* singt = getSample("T_t-channel", muon);
TH1D* singtEff = new TH1D("cut eff","cut eff",10,0,10);
TH1D* singtw = getSample("T_tW-channel", muon);
TH1D* singtwEff = new TH1D("cut eff","cut eff",10,0,10);

for(int i =1; i<10; i++){
ttEff->SetBinContent(i+1, tt->GetBinContent(i+1)/tt->GetBinContent(i));
wjetsEff->SetBinContent(i+1, wjets->GetBinContent(i+1)/wjets->GetBinContent(i));
zjetsEff->SetBinContent(i+1, zjets->GetBinContent(i+1)/zjets->GetBinContent(i));
qcdEff->SetBinContent(i+1, qcd->GetBinContent(i+1)/qcd->GetBinContent(i));
singtEff->SetBinContent(i+1, singt->GetBinContent(i+1)/singt->GetBinContent(i));
singtwEff->SetBinContent(i+1, singtw->GetBinContent(i+1)/singtw->GetBinContent(i));

}	

//weight by xsect and lumi etc..
if(scale == true){

weightHisto(tt, lumi, 225.2, 6920475);
weightHisto(wjets, lumi,37509 , 57708550);
weightHisto(zjets, lumi, 5745.25, 30457954);
//weightHisto(qcd, lumi, 84679.3, 8500505);
weightHisto(singt, lumi, 87.1, 3757707);
weightHisto(singtw, lumi, 22.2, 497395);
}

 //draw histos to files
  TCanvas *c1 = new TCanvas("cutflow","cutflow",900, 600);
        
	c1->SetLogy(1);
	wjets->GetYaxis()->SetTitle("Events");
	wjets->GetXaxis()->SetTitle("Selection Step");
	
	wjets->SetMinimum(50);
	wjets->SetLineColor(kGreen);
	wjets->Draw();
//	qcd->SetLineColor(kYellow);
//	qcd->Draw("same");
  	tt->SetLineColor(kRed);
	tt->Draw("same");
	zjets->SetLineColor(kBlue);
	zjets->Draw("same");
	singt->SetLineColor(kMagenta);
	singt->Draw("same");
	singtw->SetLineColor(kPink+2);
	singtw->Draw("same");
	
	TLegend *tleg;
	  tleg = new TLegend(0.2,0.15,0.5,0.45);
	  tleg->SetTextSize(0.04);
	  tleg->SetBorderSize(0);
	  tleg->SetFillColor(10);
	  tleg->AddEntry(tt , "t#bar{t}"      , "l");
	  tleg->AddEntry(wjets , "w+jets"      , "l");
	  tleg->AddEntry(zjets , "z+jets"      , "l");
	  tleg->AddEntry(singt, "single-t"      , "l");
	  tleg->AddEntry(singtw, "single-tW"      , "l");
 	tleg->Draw("same");	
	
	c1->SaveAs("plots/cutFlow/cutFlow.png");

  TCanvas *c2 = new TCanvas("cutflow eff","cutflow eff",900, 600);
	
	ttEff->GetYaxis()->SetTitle("Cut Efficiency");
	ttEff->GetXaxis()->SetTitle("Selection Step");
	ttEff->SetLineColor(kRed);
	ttEff->Draw();
	wjetsEff->SetLineColor(kGreen);
	wjetsEff->Draw("same");
	//qcdEff->SetLineColor(kYellow);
	//qcdEff->Draw("same");
  	ttEff->SetLineColor(kRed);
	zjetsEff->SetLineColor(kBlue);
	zjetsEff->Draw("same");
	singtEff->SetLineColor(kMagenta);
	singtEff->Draw("same");
	singtwEff->SetLineColor(kPink+2);
	singtwEff->Draw("same");
	
	TLegend *tleg2;
	  tleg2 = new TLegend(0.4,0.15,0.7,0.45);
	  tleg2->SetTextSize(0.04);
	  tleg2->SetBorderSize(0);
	  tleg2->SetFillColor(10);
	  tleg2->AddEntry(ttEff , "t#bar{t}"      , "l");
	  tleg2->AddEntry(wjetsEff , "w+jets"      , "l");
	  tleg2->AddEntry(zjetsEff , "z+jets"      , "l");
	  tleg2->AddEntry(singtEff, "single-t"      , "l");
	  tleg2->AddEntry(singtwEff, "single-tW"      , "l");
 	tleg2->Draw("same");	
	
	c2->SaveAs("plots/cutFlow/cutEff.png");

//write out cutflow
if(muon == true){

printCutflow("skim &  ", 1, tt);
printCutflow("trigger and clean &  ", 2, tt);
printCutflow("1==mu &  ", 3, tt);
printCutflow("mu veto &  ", 4, tt);
printCutflow("e veto &  ", 5, tt);
printCutflow("$\\geq$3 jets &  ", 8, tt);
printCutflow("$\\geq$4 jets &  ", 9, tt);
printCutflow("$\\geq$1 btags &  ", 10, tt);
printCutflow("$\\geq$2 btags &  ", 11, tt);

}else{

printCutflow("skim &  ", 1, tt);
printCutflow("trigger and clean &  ", 2, tt);
printCutflow("1==ele &  ", 3, tt);
printCutflow("mu veto &  ", 4, tt);
printCutflow("e veto &  ", 5, tt);
printCutflow("conv. rej &  ", 6, tt);
printCutflow("$\\geq$3 jets &  ", 9, tt);
printCutflow("$\\geq$4 jets &  ", 10, tt);
printCutflow("$\\geq$1 btags &  ", 11, tt);
printCutflow("$\\geq$2 btags &  ", 12, tt);

}
}

TH1D* getSample(TString sample, bool muon){
	TString dir = "rootFiles/";
	TFile* tt_file = new TFile(dir + sample + "_5050pb_PFElectron_PFMuon_PF2PATJets_PFMET_FULL.root");
	TDirectoryFile* tt_folder = (TDirectoryFile*) tt_file->Get("EventCount");
	TH1D* tt_cutflow;
	if(muon == true){
	tt_cutflow = (TH1D*) tt_folder->Get("TTbarMuPlusJetsRefSelection");
	}else{
  	  tt_cutflow = (TH1D*) tt_folder->Get("TTbarEplusJetsRefSelection");
	}
    return tt_cutflow;
}

void printCutflow(string step, int bin, TH1D* ttbar){
std::cout  << setprecision(6)<< step  << ttbar->GetBinContent(bin) << "  $\\pm$ " << setprecision(3) << ttbar->GetBinError(bin) << "  \\\\ " << std::endl;
}

void weightHisto(TH1D* sample, double lumi, double xsect, double NEvts){
sample->Scale(lumi*xsect/NEvts);
}
