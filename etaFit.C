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

void etaFit();
TH1D* getSample(TString sample, double weight);
TText* doPrelim(float x, float y);

void fcn(int& npar, double* deriv, double& f, double par[], int flag);

double lumi = 5800;
//stuff to choose
bool logPlot = false; //true for log plot
int rebinFact = 10;

//isolation selection
//TString Isolation = "QCD No Iso/";
TString Isolation = "Ref selection/";
//TString Isolation = "QCD mu+jets PFRelIso/";
//TString Isolation = "QCD non iso mu+jets/";

// number of btags
TString Nbtags = "2btags";  //standard  "2btags" , qcd "0btag"

bool inclZ = false;
bool inclW = false;
bool inclQ = true;
//choose object
TString Obj = "Muon/";
//TString Obj = "MET/";

//muon variables
const int N = 1;
TString Variable;
TString Variables[N] = {"muon_AbsEta_"};
double MinXs[N] = {0};
double MaxXs[N] = {2.6};
TString XTitles[N] = {"#left|#eta#right|_{#mu}"};

//global histos for fit
TH1D* data; 
TH1D* top_fit;
TH1D* wjets_fit;
TH1D* zjets_fit;
TH1D* qcd_fit;

double Nwjets, Nzjets, NQCD;

void etaFit(){
setTDRStyle();

//loop over variables
for(int i = 0; i<N; i++){
double MinX = MinXs[i];
double MaxX = MaxXs[i];
Variable = Variables[i];
TString Xtitle = XTitles[i];

//Data
data = getSample("SingleMu", 1);

//MC
TH1D* tt = getSample("TTJet", lumi*225.2/6920475);

TH1D* wjets = getSample("WJetsToLNu", lumi*37509/57708550);
TH1D* w1jets = getSample("W1Jet", lumi*5400.0/23140779);
TH1D* w2jets = getSample("W2Jets", lumi*1750.0/34041404);
TH1D* w3jets = getSample("W3Jets", lumi*519.0/15536443);
TH1D* w4jets = getSample("W4Jets", lumi*214.0/13370904);

TH1D* zjets = getSample("DYJetsToLL", lumi*5745.25/30457954);
TH1D* z1jets = getSample("DY1JetsToLL", lumi*561.0/24042904);
TH1D* z2jets = getSample("DY2JetsToLL", lumi*181.0/21835749);
TH1D* z3jets = getSample("DY3JetsToLL", lumi*51.1/11010628);
TH1D* z4jets = getSample("DY4JetsToLL", lumi*23.04/6391785);

TH1D* qcd = getSample("QCD_Pt_20_MuEnrichedPt_15",     lumi*34679.3/8500505);
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

TH1D* top_t = getSample("T_t-channel", lumi*56.4/3757707);
TH1D* top_tw = getSample("T_tW-channel", lumi*11.1/497395);
TH1D* top_s = getSample("T_s-channel", lumi*3.79/249516);
TH1D* tbar_t = getSample("Tbar_t-channel", lumi*30.7/1934817);
TH1D* tbar_tw = getSample("Tbar_tW-channel", lumi*11.1/493239);
TH1D* tbar_s = getSample("Tbar_s-channel", lumi*1.76/139948);

TH1D* top = (TH1D*)tt->Clone("top");
top->Add(top_t); top->Add(top_tw);top->Add(top_s); top->Add(tbar_t); top->Add(tbar_tw);top->Add(tbar_s);
TH1D* allMC = (TH1D*)top->Clone("allMC");
allMC->Add(wjets); allMC->Add(zjets); allMC->Add(qcd);

top_fit = (TH1D*)top->Clone("top fit");
wjets_fit = (TH1D*)wjets->Clone("wjets fit");
zjets_fit = (TH1D*)zjets->Clone("zjets fit");
qcd_fit = (TH1D*)qcd->Clone("qcd fit");
top_fit->Scale(1./ top_fit->Integral());
wjets_fit->Scale(1./ wjets_fit->Integral()); 
zjets_fit->Scale(1./ zjets_fit->Integral()); 
qcd_fit->Scale(1./ qcd_fit->Integral());
  
  //draw histos to files
  TCanvas *c2 = new TCanvas("Plot","Plot",900, 600);
  
  top_fit->SetFillColor(kWhite); wjets_fit->SetFillColor(kWhite); zjets_fit->SetFillColor(kWhite); qcd_fit->SetFillColor(kWhite);
  top_fit->Draw();
  wjets_fit->Draw("same");
  zjets_fit->Draw("same");
  qcd_fit->Draw("same");
  c2->SaveAs("plots/Fits/Templates.png");
  delete c2;
  
  
THStack *hs = new THStack("hs","test");
  if(inclQ == true){
  hs->Add(qcd);
  }else{
  hs->Add(qcd1);
  hs->Add(qcd2);
  hs->Add(qcd3);
  hs->Add(qcd4);
  hs->Add(qcd5);
  hs->Add(qcd6);
  hs->Add(qcd7);
  hs->Add(qcd8);
  hs->Add(qcd9);
  hs->Add(qcd10);
  hs->Add(qcd11);
  }
  
  
  if(inclZ == true){
  hs->Add(zjets);
  }else{
  hs->Add(z1jets);
  hs->Add(z2jets);
  hs->Add(z3jets);
  hs->Add(z4jets);  
  }
  
  if(inclW == true){
  hs->Add(wjets);
  }else{
  hs->Add(w1jets);
  hs->Add(w2jets);
  hs->Add(w3jets);
  hs->Add(w4jets);  
  }
      
  hs->Add(top_t);
  hs->Add(top_tw);
  hs->Add(top_s);
  hs->Add(tbar_t);
  hs->Add(tbar_tw);
  hs->Add(tbar_s);
  
  hs->Add(tt);

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
	
  if(logPlot ==true){
  c1->SetLogy();
  }	
  
  TString plotName("plots/Control/Muon/");
  
  if(logPlot ==true){
    plotName += Variable+"Test_Log";
    plotName += Nbtags+".png";
    
  }else{
    plotName += Variable+"Test";  
    plotName += Nbtags+".png";
  }
 
 
  c1->SaveAs(plotName);
  delete c1;
  
int Ntotal = data->Integral();
double Nsignal = top->Integral();
Nwjets = wjets->Integral();
Nzjets = zjets->Integral();
NQCD = qcd->Integral();

  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 4;              // the number of parameters
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);

  minuit.SetPrintLevel(1);
  minuit.SetErrorDef(1.);
  
  
  int ierflg = 0;
  string parName[npar] = {"ttbar+single-top", "wjets", "zjets", "qcd"}; //background parameters
  double par[npar] = {top->Integral(), wjets->Integral(), zjets->Integral(), qcd->Integral()};               //using the MC estimation as the start values 1fb
  
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
	wjets_fit->SetBinContent(i+1, wjets_fit->GetBinContent(i+1)*outpar[1]);
	zjets_fit->SetBinContent(i+1, zjets_fit->GetBinContent(i+1)*outpar[2]);
	qcd_fit->SetBinContent(i+1, qcd_fit->GetBinContent(i+1)*outpar[3]);
  
  }
  //print out the results
  cout <<" \n Total number of events after the fit" << endl;
  cout<<"   & ttbar+single top & w+jets & z+jets & qcd "<<endl;
  cout <<  " & " << Nsignal <<  " & " << Nwjets << " & " <<  Nzjets << " & " <<  NQCD  <<endl;
  cout<< " & "<<outpar[0]<<"+-"<<err[0]<<" & "<<outpar[1]<<"+-"<<err[1]<<" & "<<outpar[2]<<"+-"<<err[2]<<" & "<<outpar[3]<<"+-"<<err[3]<<endl;
  
  
    TCanvas *c3 = new TCanvas("Plot","Plot",900, 600);
  
  THStack* sum_fit = new THStack("sum fit","stacked histograms"); //used for stack plot
  qcd_fit->SetFillColor(kYellow); zjets_fit->SetFillColor(kBlue);   wjets_fit->SetFillColor(kGreen);  top_fit->SetFillColor(kRed);
  sum_fit->Add(qcd_fit); sum_fit->Add(zjets_fit);  sum_fit->Add(wjets_fit);  sum_fit->Add(top_fit);
  

  sum_fit->Draw();
  data->Draw("E same");
  
   c3->SaveAs("plots/Fits/Fit.png");
    delete c3;
  
  
  }//end of loop over variables
  	
}


TH1D* getSample(TString sample, double weight){
	TString dir = "rootFiles/";
	TFile* file = new TFile(dir + sample + "_10000pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	//TDirectoryFile* folder = (TDirectoryFile*) file->Get("TTbarPlusMetAnalysis/QCD No Iso/Muon/");
	
	TH1D* plot = (TH1D*) file->Get("TTbarPlusMetAnalysis/MuPlusJets/"+Isolation+Obj+Variable+"2btags");
	TH1D* plot2 = (TH1D*) file->Get("TTbarPlusMetAnalysis/MuPlusJets/"+Isolation+Obj+Variable+"3btags");
	TH1D* plot3 = (TH1D*) file->Get("TTbarPlusMetAnalysis/MuPlusJets/"+Isolation+Obj+Variable+"4orMoreBtags");

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
	}
	
	plot->Scale(weight);
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

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

   double lnL = 0.0;


  for (int i=0; i< top_fit->GetNbinsX(); i++){

    //data_i is the observed number of events in each bin
    int data_i = data->GetBinContent(i+1);
    //xi is the expected number of events in each bin
    double xi = par[0]*top_fit->GetBinContent(i+1) + par[1]*wjets_fit->GetBinContent(i+1) + par[2]*zjets_fit->GetBinContent(i+1) + par[3]*qcd_fit->GetBinContent(i+1);

   
    if(data_i !=0 && xi != 0){
      lnL += log(TMath::Poisson(data_i, xi));
    }
    
  }

  //W+jets, Z+jets constraints
  f = -2.0 * lnL;
  

  double nwjets = Nwjets;
//  double nwjets_err = nwjets*0.3;
  //double nwjets_err = nwjets*0.02;
   
  double nzjets = Nzjets;
//  double nzjets_err = nzjets*0.1;

double nqcd = 0;
if(NQCD>0){
  nqcd = NQCD;
  }
  else{
  nqcd = 0.00000001;
  }
  double nqcd_err = nqcd*1.;

  //ratio constraints
   f += ( (par[2]/par[1] - nzjets/nwjets) / (0.05 *nzjets/nwjets) )  * ( (par[2]/par[1] - nzjets/nwjets) / (0.05*nzjets/nwjets) ); 


   f += ((par[3]-nqcd)*(par[3]-nqcd))/nqcd_err/nqcd_err;



}                         

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

