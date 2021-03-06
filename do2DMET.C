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
#include <iomanip>
#include <sstream>
#include "TLegend.h"
using namespace std;

void do2DMET();
void do2DPlots(bool muon, TString variable, TString ytitle);
void getBinning(bool muon, TString variable);


void do2DMET(){

bool muon = true;

TString variable[41] = {"patMETsPFlow" , "GenMET" , "patType1CorrectedPFMet" , "patType1p2CorrectedPFMet" , "patType1p2CorrectedPFMetElectronEnUp" , "patType1p2CorrectedPFMetElectronEnDown" , "patType1p2CorrectedPFMetMuonEnUp" , "patType1p2CorrectedPFMetMuonEnDown" , "patType1p2CorrectedPFMetTauEnUp" , "patType1p2CorrectedPFMetTauEnDown" , "patType1p2CorrectedPFMetJetResUp" , "patType1p2CorrectedPFMetJetResDown" , "patType1p2CorrectedPFMetJetEnUp" , "patType1p2CorrectedPFMetJetEnDown" , "patType1p2CorrectedPFMetUnclusteredEnUp" , "patType1p2CorrectedPFMetUnclusteredEnDown" , "patPFMetElectronEnUp" , "patPFMetElectronEnDown" , "patPFMetMuonEnUp" , "patPFMetMuonEnDown" , "patPFMetTauEnUp" , "patPFMetTauEnDown" , "patPFMetJetResUp" , "patPFMetJetResDown" , "patPFMetJetEnUp" , "patPFMetJetEnDown" , "patPFMetUnclusteredEnUp" , "patPFMetUnclusteredEnDown" , "patType1CorrectedPFMetElectronEnUp" , "patType1CorrectedPFMetElectronEnDown" , "patType1CorrectedPFMetMuonEnUp" , "patType1CorrectedPFMetMuonEnDown" , "patType1CorrectedPFMetTauEnUp" , "patType1CorrectedPFMetTauEnDown" , "patType1CorrectedPFMetJetResUp" , "patType1CorrectedPFMetJetResDown" , "patType1CorrectedPFMetJetEnUp" , "patType1CorrectedPFMetJetEnDown" , "patType1CorrectedPFMetUnclusteredEnUp" , "patType1CorrectedPFMetUnclusteredEnDown" , "recoMetPFlow"};

TString ytitle[41] = {"patMETsPFlow" , " GenMET" , " Type1 Corrected PFMet (GeV)" , " patType1p2CorrectedPFMet" , " patType1p2CorrectedPFMetElectronEnUp" , "patType1p2CorrectedPFMetElectronEnDown" , "patType1p2CorrectedPFMetMuonEnUp" , "patType1p2CorrectedPFMetMuonEnDown" , "patType1p2CorrectedPFMetTauEnUp" , "patType1p2CorrectedPFMetTauEnDown" , "patType1p2CorrectedPFMetJetResUp" , "patType1p2CorrectedPFMetJetResDown" , "patType1p2CorrectedPFMetJetEnUp" , "patType1p2CorrectedPFMetJetEnDown" , "patType1p2CorrectedPFMetUnclusteredEnUp" , "patType1p2CorrectedPFMetUnclusteredEnDown" , "patPFMetElectronEnUp" , "patPFMetElectronEnDown" , "patPFMetMuonEnUp" , "patPFMetMuonEnDown" , "patPFMetTauEnUp" , "patPFMetTauEnDown" , "patPFMetJetResUp" , "patPFMetJetResDown" , "patPFMetJetEnUp" , "patPFMetJetEnDown" , "patPFMetUnclusteredEnUp" , "patPFMetUnclusteredEnDown" , "patType1CorrectedPFMetElectronEnUp" , "patType1CorrectedPFMetElectronEnDown" , "patType1CorrectedPFMetMuonEnUp" , "patType1CorrectedPFMetMuonEnDown" , "patType1CorrectedPFMetTauEnUp" , "patType1CorrectedPFMetTauEnDown" , "patType1CorrectedPFMetJetResUp" , "patType1CorrectedPFMetJetResDown" , "patType1CorrectedPFMetJetEnUp" , "patType1CorrectedPFMetJetEnDown" , "patType1CorrectedPFMetUnclusteredEnUp" , "patType1CorrectedPFMetUnclusteredEnDown" , "recoMetPFlow"};



for(int i =2; i<3; i++){
do2DPlots(muon, variable[i], ytitle[i]);
getBinning(muon, variable[i]);
}

}

void do2DPlots(bool muon, TString variable, TString ytitle){

	TString leptonFolder;
	if(muon == true){
		leptonFolder = "MuonMET/";
	}else{
		leptonFolder = "ElectronMET/";
		}	
		
  	setTDRStyle();
  	gStyle->SetPalette(1);

	TString dir = "rootFilesV2/central/";
	TFile* tt_file = new TFile(dir + "TTJet_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");

cout << "Getting histo: " << "Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET"<<endl;
	TH2D* tt_2d = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET_2btags");
	TH2D* tt_2d_3b = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET_3btags");
	TH2D* tt_2d_4b = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET_4orMoreBtags");

	tt_2d->Add(tt_2d_3b);
	tt_2d->Add(tt_2d_4b);

	tt_2d->Rebin2D(10,10);
  	tt_2d->GetYaxis()->SetTitle(ytitle);
  	tt_2d->GetXaxis()->SetTitle("Gen MET (GeV)");
	tt_2d->GetYaxis()->SetTitleOffset(1.8);
	tt_2d->GetXaxis()->SetTitleOffset(1.5);

        TCanvas *c= new TCanvas("c","c",10,10,800,600);
	tt_2d->Draw("COLZ");
			
	int bin[5] = {25, 45, 70, 100, 150};//MET
	for(int i = 0; i < 5; i++){
	TLine *line = new TLine(bin[i],0,bin[i],300);
	TLine *liney = new TLine(0,bin[i],300,bin[i]);
	//TLine *line = new TLine(bin[i],0,bin[i],500);
	//TLine *liney = new TLine(0,bin[i],500,bin[i]);
	line->SetLineWidth(2);
	liney->SetLineWidth(2);
	liney->Draw();
	line->Draw();
	}
	
  	TString plotName("plots/"+leptonFolder);
        plotName += variable;
        plotName += "_2btags.pdf";
 
  c->SaveAs(plotName);
  delete c;


}


void getBinning(bool muon, TString variable){

	TString leptonFolder;
	if(muon == true){
		leptonFolder = "MuonMET/";
	}else{
		leptonFolder = "ElectronMET/";
		}

	TString dir = "rootFilesV2/central/";

TFile* tt_file = new TFile(dir + "TTJet_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
TH2D* tt_2d = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET_2btags");
TH2D* tt_2d_3b = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET_3btags");
TH2D* tt_2d_4b = (TH2D*) tt_file->Get("Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET_4orMoreBtags");

tt_2d->Add(tt_2d_3b);
tt_2d->Add(tt_2d_4b);

cout << "Getting histo: " << "Binning/"+leptonFolder+variable+"/GenMET_vs_RecoMET_2btags" <<endl;


int bin2 = 25;
int bin3 = 45;
int bin4 = 70;
int bin5 = 100;
int bin6 = 150;
int binEnd =301;

//double weight  = 20000*210.5/6920475;
//double weight = 5050*157.5/6920475;

  double Evts1 = tt_2d->Integral(0,301,0,bin2-1);
  double Evts2 = tt_2d->Integral(0,301,bin2,bin3-1);
  double Evts3 = tt_2d->Integral(0,301,bin3,bin4-1);
  double Evts4 = tt_2d->Integral(0,301,bin4,bin5-1);
  double Evts5 = tt_2d->Integral(0,301,bin5,bin6-1);
  double Evts6 = tt_2d->Integral(0,301,bin6,301);

  double p1 = tt_2d->Integral(0,bin2-1,0,bin2-1)/tt_2d->Integral(0,301,0,bin2-1);
  double p2 = tt_2d->Integral(bin2,bin3-1,bin2,bin3-1)/tt_2d->Integral(0,301,bin2,bin3-1);
  double p3 = tt_2d->Integral(bin3,bin4-1,bin3,bin4-1)/tt_2d->Integral(0,301,bin3,bin4-1);
  double p4 = tt_2d->Integral(bin4,bin5-1,bin4,bin5-1)/tt_2d->Integral(0,301,bin4,bin5-1);
  double p5 = tt_2d->Integral(bin5,bin6-1,bin5,bin6-1)/tt_2d->Integral(0,301,bin5,bin6-1);
  double p6 = tt_2d->Integral(bin6,301,bin6,301)/tt_2d->Integral(0,301,bin6,301);
  
  double s1 = tt_2d->Integral(0,bin2-1,0,bin2-1)/tt_2d->Integral(0,bin2-1,0,301);
  double s2 = tt_2d->Integral(bin2,bin3-1,bin2,bin3-1)/tt_2d->Integral(bin2,bin3-1,0,301);
  double s3 = tt_2d->Integral(bin3,bin4-1,bin3,bin4-1)/tt_2d->Integral(bin3,bin4-1,0,301);
  double s4 = tt_2d->Integral(bin4,bin5-1,bin4,bin5-1)/tt_2d->Integral(bin4,bin5-1,0,301);
  double s5 = tt_2d->Integral(bin5,bin6-1,bin5,bin6-1)/tt_2d->Integral(bin5,bin6-1,0,301);
  double s6 = tt_2d->Integral(bin6,301,bin6,301)/tt_2d->Integral(bin6,301,0,301);
  
  double dp1 = p1*sqrt(pow(tt_2d->Integral(0,bin2-1,0,bin2-1),-1)+pow(tt_2d->Integral(0,301,0,bin2-1),-1));
  double dp2 = p2*sqrt(pow(tt_2d->Integral(bin2,bin3-1,bin2,bin3-1),-1)+pow(tt_2d->Integral(0,301,bin2,bin3-1),-1));
  double dp3 = p3*sqrt(pow(tt_2d->Integral(bin3,bin4-1,bin3,bin4-1),-1)+pow(tt_2d->Integral(0,301,bin3,bin4-1),-1));
  double dp4 = p4*sqrt(pow(tt_2d->Integral(bin4,bin5-1,bin4,bin5-1),-1)+pow(tt_2d->Integral(0,301,bin4,bin5-1),-1));
  double dp5 = p5*sqrt(pow(tt_2d->Integral(bin5,301,bin5,301),-1)+pow(tt_2d->Integral(0,301,bin5,301),-1));

  double ds1 = s1*sqrt(pow(tt_2d->Integral(0,bin2-1,0,bin2-1),-1)+pow(tt_2d->Integral(0,bin2-1,0,301),-1));
  double ds2 = s2*sqrt(pow(tt_2d->Integral(bin2,bin3-1,bin2,bin3-1),-1)+pow(tt_2d->Integral(bin2,bin3-1,0,301),-1));
  double ds3 = s3*sqrt(pow(tt_2d->Integral(bin3,bin4-1,bin3,bin4-1),-1)+pow(tt_2d->Integral(bin3,bin4-1,0,301),-1));
  double ds4 = s4*sqrt(pow(tt_2d->Integral(bin4,bin5-1,bin4,bin5-1),-1)+pow(tt_2d->Integral(bin4,bin5-1,0,301),-1));
  double ds5 = s5*sqrt(pow(tt_2d->Integral(bin5,301,bin5,301),-1)+pow(tt_2d->Integral(bin5,301,0,301),-1));

//cout << " bin1 stab: " << s1 << " bin2 pure: " << p1 << endl;
//cout << " bin2 stab: " << s2 << " bin2 pure: " << p2 << endl;
//cout << " bin3 stab: " << s3 << " bin3 pure: " << p3 << endl;
//cout << " bin4 stab: " << s4 << " bin4 pure: " << p4 << endl;
//cout << " bin5 stab: " << s5 << " bin5 pure: " << p5 << endl;

cout << setprecision(4) << "MC Evts & " << Evts1 << " & " << Evts2 << " & " << Evts3 << " & " << Evts4 << " & " << Evts5 << " & " << Evts6<< endl;
cout << setprecision(2) << "purity & " << p1 << " & " << p2 << " & " << p3 << " & " << p4 << " & " << p5<< " & " << p6 << endl;
cout << setprecision(2) << "stability & " << s1 << " & " << s2 << " & " << s3 << " & " << s4 << " & " << s5 << " & " << s6<< endl;

}
