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
#include "TROOT.h"

TH1D* getSample(TString sample, double weight, TString Obj, TString Variable, TString Isolation, int RebinFact, TString Systematic);
TH1D* getQCD(TString Obj, TString Variable, int RebinFact);
TText* doPrelim(float x, float y);


TH1D* getSample(TString sample, double weight, TString Obj, TString Variable, TString Isolation, int rebinFact, TString Systematic){
	TString dir = "rootFilesV2/"+ Systematic +"/";
	
	TString syst = "";
	
	if(Systematic == "BJet_down")
		syst = "_minusBJet";
	else if(Systematic == "BJet_up")
		syst = "_plusBjet";
	else if(Systematic == "JES_down")
		syst = "_minusJES";
	else if(Systematic == "JES_up")
		syst = "_plusJES";
	else if(Systematic == "LightJet_up")
		syst = "_plusLightJet";	
	else if(Systematic == "LightJet_down")
		syst = "_minusLightJet";							
	else if(Systematic == "PU_down")
		syst = "_PU_65835mb";
	else if(Systematic == "PU_up")
		syst = "_PU_72765mb";
	else	
		syst = "";
		
	TFile* file = new TFile(dir + sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET"+syst+".root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");

	cout << dir + sample + "_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET"+syst+".root" << endl;

	TH1D* plot = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/"+Isolation+Obj+Variable+"2btags");
	TH1D* plot2 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/"+Isolation+Obj+Variable+"3btags");
	TH1D* plot3 = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/"+Isolation+Obj+Variable+"4orMoreBtags");

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
	}else if(sample == "QCD_Pt_20_MuEnrichedPt_15" || sample == "QCD_Pt-15to20_MuEnrichedPt5" || sample=="QCD_Pt-15to20_MuEnrichedPt5" || sample =="QCD_Pt-20to30_MuEnrichedPt5" || sample ==    "QCD_Pt-30to50_MuEnrichedPt5" || sample ==    "QCD_Pt-50to80_MuEnrichedPt5" || sample ==    "QCD_Pt-80to120_MuEnrichedPt5" || sample ==   "QCD_Pt-120to170_MuEnrichedPt5" || sample ==  "QCD_Pt-170to300_MuEnrichedPt5" || sample ==  "QCD_Pt-300to470_MuEnrichedPt5" || sample ==  "QCD_Pt-470to600_MuEnrichedPt5" || sample ==  "QCD_Pt-800to1000_MuEnrichedPt5" || sample =="QCD_Pt-1000_MuEnrichedPt5" 	 ){
	plot->SetFillColor(kYellow);
	plot->SetLineColor(kYellow);
	}else if(sample == "T_t-channel" || sample == "T_tW-channel" || sample == "T_s-channel" || sample == "Tbar_t-channel" || sample == "Tbar_tW-channel" || sample == "Tbar_s-channel"){
	plot->SetFillColor(kMagenta);
	plot->SetLineColor(kMagenta);
	}else if(sample == "SingleMu"){
	plot->SetLineColor(kBlack);	  
	}

	//plot->Scale(weight);
	plot->Rebin(rebinFact);
	
	plot->SetDirectory(gROOT);
	file->Close();
	
	return plot;

}

//always data
TH1D* getQCD(TString Obj, TString Variable, int rebinFact){
	TString dir = "rootFilesV2/central/";
	TFile* file = new TFile(dir + "SingleMu_5814pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");
	
	TH1D* plot = (TH1D*) file->Get("TTbar_plus_X_analysis/MuPlusJets/QCD non iso mu+jets/"+Obj+Variable+"0btag");

	plot->SetFillColor(kYellow);
	plot->SetLineColor(kYellow);
    	
	plot->Scale(1/plot->Integral());
	plot->Rebin(rebinFact);
	
	plot->SetDirectory(gROOT);
	file->Close();
	
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
