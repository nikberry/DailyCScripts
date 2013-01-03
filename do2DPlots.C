#include <iomanip>
#include "tdrstyle.C"
#include "TTree.h"
#include "TFile.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLine.h"
#include "TChain.h"
#include "TLatex.h"
#include <iostream>
#include <sstream>
#include "TLegend.h"
using namespace std;

void do2DPlots();
void do2DPlots(bool muon, TString variable, TString xtitle, TString ytitle);
void getBinning(bool muon, TString variable, TString xtitle, TString ytitle);
void writeStab(TH2D * tt_2d, string write, bool choice);
void do2DPlots(){

bool muon = false;

TString variable[10] = {"GenHT_vs_RecoHT", "GenHT_lep_vs_RecoHT_lep", "GenHT_lep_met_vs_RecoHT_lep_met", "Genleptonic_W_pt_vs_Recoleptonic_W_pt", "GenM3_vs_RecoM3", "GenNu_vs_RecoMET", "GenParton_vs_RecoHT", "GenJetHT_vs_GenParton", "GenLepPlusMETPt_vs_RecoLepPlusMetPt", "GendPhiLepMet_vs_RecodPhiLepMetPt"};
TString xtitle[10] = {"Gen HT (GeV)", "Gen HT+lep (GeV)", "Gen HT+lep+met (GeV)", "Gen leptonic W pt (GeV)", "Gen M3 (GeV)", "Reco Nu (GeV)", "Gen HT GeVReco Parton HT (GeV)", "Gen Jet HT (GeV)", "Gen #mu pt + MET  (GeV)", "Gen #Delta#Phi(#mu,MET)"};
TString ytitle[10] = {"Reco HT (GeV)", "Reco HT+lep (GeV)", "Reco HT+lep+met (GeV)", "Reco leptonic W pt (GeV)", "Reco M3 (GeV)", "Reco MET (GeV)", "Reco HT (GeV)", "Gen Parton HT (GeV)", "Reco #mu pt + MET  (GeV)", "Reco #Delta#Phi(#mu,MET)"};



for(int i =2; i<3; i++){
//do2DPlots(muon, variable[i], xtitle[i], ytitle[i]);
getBinning(muon, variable[i], xtitle[i], ytitle[i]);
}

}

void do2DPlots(bool muon, TString variable, TString xtitle, TString ytitle){

	TString leptonFolder;
	if(muon == true){
		leptonFolder = "MuPlusJets/";
	}else{
		leptonFolder = "EPlusJets/";
		}	
		
  	setTDRStyle();
  	gStyle->SetPalette(1);

	TString dir = "rootFiles/";
	TFile* tt_file = new TFile(dir + "TTJet_5050pb_PFElectron_PFMuon_PF2PATJets_PFMET_FULL.root");  



TString Nbtags[5] = {"_0btag","_1btag", "_2btags","_3btags",  "_4orMoreBtags"};

for(int i = 2; i < 3; i++){

	TH2D* tt_2d = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+Nbtags[i]);

	tt_2d->Rebin2D(10,10);
  	tt_2d->GetYaxis()->SetTitle(ytitle);
  	tt_2d->GetXaxis()->SetTitle(xtitle); 
	tt_2d->GetYaxis()->SetTitleOffset(1.8);
	tt_2d->GetXaxis()->SetTitleOffset(1.5);	

        TCanvas *c= new TCanvas("c","c",10,10,800,600);
	tt_2d->Draw("COLZ");
	
  	TString plotName("plots/"+leptonFolder);
        plotName += variable;
        plotName += Nbtags[i]+".png";
 
  c->SaveAs(plotName);
  delete c;
	
}


}


void getBinning(bool muon, TString variable, TString xtitle, TString ytitle){

	TString leptonFolder;
	if(muon == true){
		leptonFolder = "MuPlusJets/";
	}else{
		leptonFolder = "EPlusJets/";
		}

	TString dir = "rootFiles/";
	TFile* tt_file = new TFile(dir + "TTJet_5050pb_PFElectron_PFMuon_PF2PATJets_PFMET_FULL.root");

	TH2D* tt_2d = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+"_2btags");

//writeStab(tt_2d, "bin", false);
	
writeStab(tt_2d, "bin", true);
writeStab(tt_2d, "events", true);
writeStab(tt_2d, "eventsWeight", true);
writeStab(tt_2d, "purity", true);
writeStab(tt_2d, "stability", true);

}


void writeStab(TH2D * tt_2d, string write, bool choice){

int i = 0;
int binMin[60];
binMin[0] = 0;
//	int binCho[9] = {55, 70, 80, 95, 115, 135, 160, 200, 499}; //hT
//	int binCho[6] = {0, 60, 90, 120, 160, 240}; //mwt
//	int binCho[5] = {0, 100, 140, 190, 299}; //m3
	int binCho[11] = {75, 90, 105, 120, 135, 155, 175, 200, 225, 275, 499};

for(int bin = 0; bin<tt_2d->GetNbinsX(); bin++){
	double purity[60];
	double stability[60];

	stability[i] = tt_2d->Integral(binMin[i], bin+1,binMin[i],bin+1)/tt_2d->Integral(binMin[i], bin+1,0,tt_2d->GetNbinsX()+1);
	purity[i] = tt_2d->Integral(binMin[i], bin+1,binMin[i],bin+1)/tt_2d->Integral(0,tt_2d->GetNbinsX()+1 ,binMin[i],bin+1);

double weight  = 20000*210.5/6920475;


	if(choice == false){

	if(purity[i]>=0.6 && stability[i]>=0.6 && tt_2d->Integral(0,tt_2d->GetNbinsX()+1 ,binMin[i],bin+1)*weight >1500){
	cout << "purity: " << purity[i] << ", stability: " << stability[i] << " bin: " << bin << " evts: " << tt_2d->Integral(0,tt_2d->GetNbinsX()+1 ,binMin[i],bin+1)*weight << endl;
	i++;
	binMin[i] = bin;
	}
		
		}else{
		if(bin == binCho[i]){
		
		if(write == "bin"){
		cout << bin*4 << " & " ;
	        }else if(write == "events"){
		cout<<setprecision(2) << tt_2d->Integral(0,tt_2d->GetNbinsX()+1 ,binMin[i],bin+1) << " & " ;
		}else if(write == "eventsWeight"){
		cout<<setprecision(2) << tt_2d->Integral(0,tt_2d->GetNbinsX()+1 ,binMin[i],bin+1)*weight << " & " ;
		}else if(write == "purity"){
		cout<<setprecision(2) << purity[i] << " & " ;
		}else if(write == "stability"){
		cout<<setprecision(2) << stability[i] << " & " ;
		}
		
		i++;
		binMin[i] = bin;
		}
		
		}

}

cout << " " << endl;

}





