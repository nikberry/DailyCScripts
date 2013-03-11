#include "TFile.h"
#include "TH1.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TText.h"
#include <string.h>
#include <iostream>
#include <iomanip>
#include "tdrstyle.C"


void cutFlowTwoSamp();
TH1D* getSample(TString sample, bool muon);
void printCutflow(string step, int bin, TH1D* ttbar);
void weightHisto(TH1D* sample, double lumi, double xsect, double NEvts);

void cutFlowTwoSamp(){

	bool muon = true;
	bool scale = false;

double lumi = 5800;

setTDRStyle();
//TH1D* tt = getSample("TTJet", muon);
TH1D* tt = getSample("TTJet", muon);

TH1D* WJets = getSample("WJetsToLNu", muon);
WJets->SetLineColor(kRed);

TH1D* ttEff = new TH1D("cut eff","cut eff",10,0,10);
ttEff->SetLineColor(kBlue);
//ttEff->GetYAxis->SetTitle("Events");
ttEff->SetTitle("Cut Efficiency");

TH1D* WJetsEff = new TH1D("WJets cut eff","cut eff",10,0,10);
WJetsEff->SetLineColor(kRed);
//WJets->GetYAxis->SetTitle("Events");
WJets->SetTitle("Cut Efficiency");

for(int i =1; i<10; i++){
ttEff->SetBinContent(i+1, tt->GetBinContent(i+1)/tt->GetBinContent(i));
}	

for(int i =1; i<10; i++){
WJetsEff->SetBinContent(i+1, WJets->GetBinContent(i+1)/WJets->GetBinContent(i));
}	

//weight by xsect and lumi etc..
if(scale == true){
weightHisto(tt, lumi, 157.5, 6920475);
}

 //draw histos to files
/*  TCanvas *c1 = new TCanvas("cutflow","cutflow",600, 500);    
	
	tt->Draw();
	WJets->Draw("same");
	c1->SaveAs("plots/cutFlow/cutFlow.png");
*/
  TCanvas *c2 = new TCanvas("cutflow eff","cutflow eff",600, 500);
  	
/*	TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(ttEff , "t#bar{t} Efficiency", "lf");
        leg->AddEntry(WJetsEff , "W+Jets Efficiency", "lf"); 
 */ 	
 	TText* text = new TText(0.3,0.95,"p_T > 70");
	text->SetTextFont(1);
	
	ttEff->Draw();
	WJetsEff->Draw("same");
	text->Draw("same");
	//leg->Draw("same");
	c2->SaveAs("plots/cutFlow/cutEff.png");

if(muon == true){

printCutflow("skim &  ", 1, tt);
printCutflow("trigger and clean &  ", 2, tt);
printCutflow("1==mu &  ", 3, tt);
printCutflow("mu veto &  ", 4, tt);
printCutflow("e veto &  ", 5, tt);
printCutflow("$\\geq$1 jets &  ", 6, tt);
printCutflow("$\\geq$2 jets &  ", 7, tt);
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
	TString dir = "rootFilesV2/"; //"rootFilesV2/central/"
	//TFile* tt_file = new TFile(dir + sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	//TFile* tt_file = new TFile(dir + sample + "_5050pb_PFElectron_PFMuon_PF2PATJets_PFMET_TEST.root");
	TFile* tt_file = new TFile(dir + sample + "_5050pb_PFElectron_PFMuon_PF2PATJets_PFMET_pt70.root");
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
