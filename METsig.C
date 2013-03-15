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

void METsig();

TText* doPrelim(float x, float y);

double lumi = 5800;

bool logplot = false;

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

TString dir = "rootFilesV2/central/";

//choose object
//TString Obj = "Muon/";
TString Obj = "MET/patType1CorrectedPFMet/";
//TString Obj = "MET/GenMET/";

//met variables
TString Variable = "MET_";

TString Systematic = "central";

//Re-bin Factor. Divides the number of bins rebin factor.
int rebinFact = 1;

void METsig(){

setTDRStyle();

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

TH1D* allMET = getSample("TTJet" , lumi*225.2/6920475, Obj, Variable, Isolation, rebinFact, Systematic);

allMET->Add(wjets);
allMET->Add(w1jets);
allMET->Add(w2jets);
allMET->Add(w3jets);
allMET->Add(w4jets);
allMET->Add(zjets);
allMET->Add(z1jets);
allMET->Add(z2jets);
allMET->Add(z3jets);
allMET->Add(z4jets);
allMET->Add(qcd);
allMET->Add(qcd1);
allMET->Add(qcd2);
allMET->Add(qcd3);
allMET->Add(qcd4);
allMET->Add(qcd5);
allMET->Add(qcd6);
allMET->Add(qcd7);
allMET->Add(qcd8);
allMET->Add(qcd9);
allMET->Add(qcd10);
allMET->Add(qcd11);

TH1D* sigTest = new TH1D("sigTest", "MET Significance Test", tt->GetNbinsX(), 0, tt->GetNbinsX());
sigTest->GetYaxis()->SetTitle("Significance");
sigTest->SetFillColor(kRed+1);
sigTest->SetLineColor(kRed+1);

for(int i = 0; i<tt->GetNbinsX(); i++){
sigTest->SetBinContent(i+1, tt->Integral(i+1, tt->GetNbinsX()) / sqrt(allMET->Integral(i+1, allMET->GetNbinsX())) );
//cout << allMET->Integral(i+1, allMET->GetNbinsX()) << endl;
}

TCanvas* sigtest = new TCanvas("sigtest", "MET Significance Test", 600, 500);

sigTest->Draw();
sigtest->SaveAs("plots/cutFlow/sigTestMET.png");


}
