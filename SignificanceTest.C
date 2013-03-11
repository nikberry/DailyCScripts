#include "TFile.h"
#include "TH1.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TText.h"
#include "TMath.h"
#include <string.h>
#include <iostream>
#include <iomanip>
#include "tdrstyle.C"


void SignificanceTest();
TH1D* getSample(TString sample, bool muon, TString pt);
void printCutflow(string step, int bin, TH1D* ttbar);
void weightHisto(TH1D* sample, double lumi, double xsect, double NEvts);

void SignificanceTest(){

	bool muon = true;

setTDRStyle();

TH1D* sigTest = new TH1D("sigTest", "Jet pt Significance Test", 6, 0, 6);

const int N = 6;
TString pts[N] = {"20.root", "30.root", "40.root", "50.root", "60.root", "70.root"};

for(int i = 0; i<N; i++){
TString pt = pts[i];

TH1D* tt = getSample("TTJet", muon, pt);

TH1D* WJets = getSample("WJetsToLNu", muon, pt);

sigTest->SetBinContent(i+1, tt->GetBinContent(9) / sqrt(tt->GetBinContent(9) + WJets->GetBinContent(9)) );
}


}
*/
 //draw histos to files
  TCanvas *c1 = new TCanvas("sigTest","Cutflow Significance Test",600, 500);    
	
	sigTest->Draw();
	c1->SaveAs("plots/cutFlow/sigTest.png");

 /* TCanvas *c2 = new TCanvas("cutflow eff","cutflow eff",600, 500);
  	
	TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(ttEff , "t#bar{t} Efficiency", "lf");
        leg->AddEntry(WJetsEff , "W+Jets Efficiency", "lf"); 
  	
 	TText* text = new TText(0.3,0.95,"p_T > 70");
	text->SetTextFont(1);
	
	ttEff->Draw();
	WJetsEff->Draw("same");
	text->Draw("same");
	//leg->Draw("same");
	c2->SaveAs("plots/cutFlow/cutEff.png");
*/

}

	
TH1D* getSample(TString sample, bool muon, TString pt){
	TString dir = "rootFilesV2/"; //"rootFilesV2/central/"
	TFile* tt_file = new TFile(dir + sample + "_5050pb_PFElectron_PFMuon_PF2PATJets_PFMET_pt" + pt);
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

