#include "etaFit.C"
#include "TString.h"
#include "TLatex.h"
#include "TDirectory.h"
#include <vector>

void diffXsect();

using namespace std;

void diffXsect(){

//choice of systematic to look at
int choice = 0;

//MET will need choice of variable at the top
TString Variable ="_MET";
int Nbins = 6;
TString bins[6] = {"_bin_0-25", "_bin_25-45", "_bin_45-70", "_bin_70-100", "_bin_100-150", "_bin_150-inf"};
double width[6] = {25, 20, 25, 30, 50, 100};
double xbins[7] = {1,25,45,70,100,150, 250}; 
TString varBin = "Binned_MET_Analysis/patType1CorrectedPFMet";
TString Xtitle = "MET [GeV]";

//HT
// TString Variable ="_HT";
// int Nbins = 8;
// TString bins[8] = {"0-50", "50-150", "150-250", "250-350", "350-450", "450-650", "650-1100", "1100-inf"};
// double width[8] = {50,100,100,100,100,200,450,400};
// double xbins[9] = {1,50,150,250,350,450,650,1100, 1500}; 
// TString varBin = "Binned_HT_Analysis/HT_bin_";
// TString Xtitle = "HT [GeV]";

//ST
// TString Variable ="_ST";
// int Nbins = 8;
// TString bins[8] = {"0-150", "150-250", "250-350", "350-450", "450-550", "550-750", "750-1250", "1250-inf"};
// double width[8] = {150,100,100,100,100,200,500,500};
// double xbins[9] = {1,150,250,350,450,550,750,1250, 1750}; 
// TString varBin = "Binned_ST_Analysis/ST_with_RecoMET_bin_";
// TString Xtitle = "ST [GeV]";
//MT
// TString Variable ="_MT";
// int Nbins = 5;
// TString bins[5] = {"0-40", "40-65", "65-85", "85-150", "150-inf"};
// double width[5] = {40,25,20,65,50};
// double xbins[6] = {1,40,65,85,150,200}; 
// TString varBin = "Binned_MT_Analysis/MT_with_RecoMET_bin_";
// TString Xtitle = "M(W)_{T} [GeV]";

int Nsys = 1;
//int Nsys = 15;

//from fit
double NfitVal[Nbins][Nsys];
double NfitErr[Nbins][Nsys];
double NttbarVal[Nbins][Nsys];
double NttbarErr[Nbins][Nsys];
double sigmaVal[Nbins][Nsys];
double sigmaErr[Nbins][Nsys];

//before fit
double madgraphVals[Nbins][Nsys];
double mcatnloVals[Nbins][Nsys];
double powhegVals[Nbins][Nsys];

double ttbarPre[Nbins][Nsys];
double singletPre[Nbins][Nsys];
double wjetsPre[Nbins][Nsys];
double zjetsPre[Nbins][Nsys];
double qcdPre[Nbins][Nsys];
double BGscale[Nbins][Nsys];
double totXsect[Nsys]; 

//sample
//TString dir = "central";
//TString dir = "JES_up";
TString dirs[27] = {"central","JES_up","JES_down", "BJet_up", "BJet_down", "PU_up", "PU_down", "Scale_up_tt", "Scale_down_tt", "Scale_up", "Scale_down", "Match_up_tt", "Match_down_tt", "Match_up", "Match_down", "UnclusteredEnUp", "UnclusteredEnDown", "JetEnUp", "JetEnDown", "JetResUp", "JetResDown", "TauEnUp", "TauEnDown", "MuonEnUp", "MuonEnDown", "ElectronEnUp", "ElectronEnDown"};

//loop over systematics
for(int sys = 0; sys < Nsys; sys++){
TString dir  = dirs[sys];

totXsect[sys]= 0;

int rebinFact = 1;
//for ttbar total
TH1D* tt_tot = getSample("TTJet", lumi*225.197/6920475, rebinFact, "Muon", dir);
TH1D* mcnlo_tot = getSample("TTJet_MCNLO", lumi*225.197/6920475, rebinFact, "Muon", dir);
TH1D* powheg_tot = getSample("TTJet_POWHEG", lumi*225.197/6920475, rebinFact, "Muon", dir);
//loop over bins of distribution
for(int i = 0; i < Nbins; i++){

TString bin = varBin;

if(dir == "UnclusteredEnUp" || dir == "UnclusteredEnDown" || dir == "JetEnUp" || dir == "JetEnDown" || dir == "JetResUp" || dir == "JetResDown" || dir == "TauEnUp" || dir == "TauEnDown" || dir == "MuonEnUp" || dir == "MuonEnDown" || dir == "ElectronEnUp" || dir == "ElectronEnDown")
bin += dir;

bin += bins[i];


bool inclZ = false;
bool inclW = false;

if(dir == "Scale_up" || dir == "Scale_down" || dir == "Match_up" || dir == "Match_down"){
inclZ = true;
inclW = true;
}

cout <<  bin << endl;
TH1D* tt = getSample("TTJet", lumi*225.197/6920475, rebinFact, bin, dir);
TH1D* tt_mcnlo = getSample("TTJet_MCNLO", lumi*225.197/6920475, rebinFact, bin, dir);
TH1D* tt_powheg = getSample("TTJet_POWHEG", lumi*225.197/6920475, rebinFact, bin, dir);

TH1D* top_t = getSample("T_t-channel", lumi*56.4/3757707, rebinFact, bin, dir);
TH1D* top_tw = getSample("T_tW-channel", lumi*11.1/497395, rebinFact, bin, dir);
TH1D* top_s = getSample("T_s-channel", lumi*3.79/249516, rebinFact, bin, dir);
TH1D* tbar_t = getSample("Tbar_t-channel", lumi*30.7/1934817, rebinFact, bin, dir);
TH1D* tbar_tw = getSample("Tbar_tW-channel", lumi*11.1/493239, rebinFact, bin, dir);
TH1D* tbar_s = getSample("Tbar_s-channel", lumi*1.76/139948, rebinFact, bin, dir);


TH1D* wjets;
//TH1D* w1jets = getSample("W1Jet", lumi*5400.0/23140779, rebinFact, bin, dir);
TH1D* w2jets = getSample("W2Jets", lumi*1750.0/34041404, rebinFact, bin, dir);
TH1D* w3jets = getSample("W3Jets", lumi*519.0/15536443, rebinFact, bin, dir);
TH1D* w4jets = getSample("W4Jets", lumi*214.0/13370904, rebinFact, bin, dir);

TH1D* zjets;
//TH1D* z1jets = getSample("DY1JetsToLL", lumi*561.0/24042904, rebinFact, bin, dir);
TH1D* z2jets = getSample("DY2JetsToLL", lumi*181.0/21835749, rebinFact, bin, dir);
TH1D* z3jets = getSample("DY3JetsToLL", lumi*51.1/11010628, rebinFact, bin, dir);
TH1D* z4jets = getSample("DY4JetsToLL", lumi*23.04/6391785, rebinFact, bin, dir);

TH1D* qcd_mc = getSample("QCD_Pt-15to20_MuEnrichedPt5",   lumi*7.022e8 * 0.0039/1722678, rebinFact, bin, dir);
TH1D* qcd2 = getSample("QCD_Pt-20to30_MuEnrichedPt5",   lumi*2.87e8 * 0.0065/8486893, rebinFact, bin, dir);
TH1D* qcd3 = getSample("QCD_Pt-30to50_MuEnrichedPt5",   lumi*6.609e7 * 0.0122/8928999, rebinFact, bin, dir);
TH1D* qcd4 = getSample("QCD_Pt-50to80_MuEnrichedPt5",   lumi*8082000.0 * 0.0218/7256011, rebinFact, bin, dir);
TH1D* qcd5 = getSample("QCD_Pt-80to120_MuEnrichedPt5",  lumi*1024000.0 * 0.0395/9030624, rebinFact, bin, dir);
TH1D* qcd6 = getSample("QCD_Pt-120to170_MuEnrichedPt5", lumi*157800.0 * 0.0473/8500505, rebinFact, bin, dir);
TH1D* qcd7 = getSample("QCD_Pt-170to300_MuEnrichedPt5", lumi*34020.0 * 0.0676/7662483, rebinFact, bin, dir);
TH1D* qcd8 = getSample("QCD_Pt-300to470_MuEnrichedPt5", lumi*1757.0 * 0.0864/7797481, rebinFact, bin, dir);
TH1D* qcd9 = getSample("QCD_Pt-470to600_MuEnrichedPt5", lumi*115.2 * 0.1024/2995767, rebinFact, bin, dir);
TH1D* qcd10 = getSample("QCD_Pt-800to1000_MuEnrichedPt5",lumi*3.57 * 0.1033/4047142, rebinFact, bin, dir);
TH1D* qcd11 = getSample("QCD_Pt-1000_MuEnrichedPt5",     lumi*0.774 * 0.1097/3807263, rebinFact, bin, dir);

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

  if(inclZ == true){
  zjets = getSample("DYJetsToLL", lumi*5745.25/30457954, rebinFact, bin, dir);

  }else{
  zjets  = getSample("DY1JetsToLL", lumi*561.0/24042904, rebinFact, bin, dir);
  zjets->Add(z2jets);
  zjets->Add(z3jets);
  zjets->Add(z4jets);  
  }
  
  if(inclW == true){
  wjets = getSample("WJetsToLNu", lumi*37509/57708550, rebinFact, bin, dir);
  }else{
  wjets = getSample("W1Jet", lumi*5400.0/23140779, rebinFact, bin, dir); 
  wjets->Add(w2jets);
  wjets->Add(w3jets);
  wjets->Add(w4jets);
  }

//sum single top into one
TH1D* single_top = (TH1D*)top_t->Clone("single top");
single_top->Add(top_tw);single_top->Add(top_s); single_top->Add(tbar_t); single_top->Add(tbar_tw);single_top->Add(tbar_s);

NfitVal[i][sys] = etaFit(bin, "measured",dir);
NfitErr[i][sys] = etaFit(bin, "measuredErr",dir);
BGscale[i][sys] = etaFit(bin,"bgscale",dir);

NttbarVal[i][sys] = NfitVal[i][sys]-single_top->Integral();
NttbarErr[i][sys] = NfitErr[i][sys] ;

sigmaVal[i][sys] = ((NfitVal[i][sys]-single_top->Integral())/tt_tot->Integral())*225.197;
sigmaErr[i][sys] = (((NfitVal[i][sys]+NfitErr[i][sys]-single_top->Integral())/tt_tot->Integral())*225.197)-sigmaVal[i][sys];

//cout << "error: " << sigmaErr[i][sys] << endl;

madgraphVals[i][sys] = (tt->Integral()/tt_tot->Integral())*225.197;
mcatnloVals[i][sys] = (tt_mcnlo->Integral()/mcnlo_tot->Integral())*225.197;
powhegVals[i][sys] = (tt_powheg->Integral()/powheg_tot->Integral())*225.197;

ttbarPre[i][sys] = tt->Integral();
singletPre[i][sys] = single_top->Integral();
wjetsPre[i][sys] = wjets->Integral();
zjetsPre[i][sys] = zjets->Integral();
qcdPre[i][sys] = qcd_mc->Integral();

totXsect[sys] += sigmaVal[i][sys];
}

}

       //will need to change MET for other variables
       TFile resultsfile("outFiles/diffResults"+Variable+".root", "UPDATE", "comment");
	
	for(int sys = 0; sys < Nsys; sys++){
	TString dir  = dirs[sys];
	
	resultsfile.mkdir(dir+"_dir");
	resultsfile.cd(dir+"_dir");
	
	//measured histo will have to put name of systematic in here to write into file
	TH1D *muon_part  = new TH1D(dir+"_signal_fit", "", Nbins, xbins);  //muon
	TH1D *ttbar_fit  = new TH1D(dir+"_ttbar_fit", "", Nbins, xbins);  //muon
	TH1D *ttbar_prefit  = new TH1D(dir+"_ttbar_prefit", "", Nbins, xbins);  //muon
	TH1D *singlet_prefit  = new TH1D(dir+"_singlet_prefit", "", Nbins, xbins);  //muon
	TH1D *wjets_prefit  = new TH1D(dir+"_wjets_prefit", "", Nbins, xbins);  //muon
	TH1D *zjets_prefit  = new TH1D(dir+"_zjets_prefit", "", Nbins, xbins);  //muon
	TH1D *qcd_prefit  = new TH1D(dir+"_qcd_prefit", "", Nbins, xbins);  //muon
	TH1D *wjets_fit  = new TH1D(dir+"_wjets_fit", "", Nbins, xbins);  //mu
	TH1D *zjets_fit  = new TH1D(dir+"_zjets_fit", "", Nbins, xbins);  //mu
	TH1D *qcd_fit  = new TH1D(dir+"_qcd_fit", "", Nbins, xbins);  //muon
	
//	muon_part->SetDirectory(directory);
	//loop over bins of distribution
	for(int i = 0; i < Nbins; i++){
	muon_part->SetBinContent(i+1,NfitVal[i][sys]);	
	muon_part->SetBinError(i+1, NfitErr[i][sys]);   
	ttbar_fit->SetBinContent(i+1,NttbarVal[i][sys]);  
	ttbar_fit->SetBinError(i+1, NttbarErr[i][sys]); 
	ttbar_prefit->SetBinContent(i+1,ttbarPre[i][sys]);  
	singlet_prefit->SetBinContent(i+1,singletPre[i][sys]); 
	wjets_prefit->SetBinContent(i+1, wjetsPre[i][sys]);
	zjets_prefit->SetBinContent(i+1, zjetsPre[i][sys]);
	qcd_prefit->SetBinContent(i+1,  qcdPre[i][sys]); 
	wjets_fit->SetBinContent(i+1, wjetsPre[i][sys]*BGscale[i][sys]);
	zjets_fit->SetBinContent(i+1, zjetsPre[i][sys]*BGscale[i][sys]);
	qcd_fit->SetBinContent(i+1,  qcdPre[i][sys]*BGscale[i][sys]); 
	     
       	}
		
	resultsfile.cd();
	}
	     
       //resultsfile.Write();
       resultsfile.Close();


cout << "N fit: " << endl;
for(int i = 0; i < Nbins; i++){
cout << bins[i] << ": " << " = " << NfitVal[i][choice] << " +- " << NfitErr[i][choice] << endl;
}

cout << "partial xsects: " << endl;
for(int i = 0; i < Nbins; i++){
cout << bins[i] << ": " << " = " << sigmaVal[i][choice] << " +- " << sigmaErr[i][choice] << endl;
//cout << bins[i] << ": " << " = " << madgraphVals[i] << endl;
}

cout << "normalised xsects: " << endl;
for(int i = 0; i < Nbins; i++){
cout << bins[i] << ": " << " = " << sigmaVal[i][choice]/totXsect[choice] << endl;
}

cout << "normalised differential: " << endl;
for(int i = 0; i < Nbins; i++){
cout << bins[i] << ": " << " = " << sigmaVal[i][choice]/(totXsect[choice]*width[i]) << endl;
}
cout << "cross section is:  " <<  totXsect[choice] << endl;  

   //measured histo will have to put name of systematic in here to write into file
   TH1D *muon_part  = new TH1D("central", "", Nbins, xbins);  //muon
   
   //different generators
   TH1D *madgraph  = new TH1D("madgraph", "", Nbins, xbins); 
   TH1D *mcatnlo  = new TH1D("mcatnlo", "", Nbins, xbins);  
   TH1D *powheg  = new TH1D("powheg", "", Nbins, xbins);  
	
	madgraph->SetLineColor(kRed);
	mcatnlo->SetLineColor(kBlue);
	powheg->SetLineColor(kGreen+2);
	for(int i = 0; i < Nbins; i++){
	muon_part->SetBinContent(i+1,sigmaVal[i][choice]);	
	muon_part->SetBinError(i+1, sigmaErr[i][choice]);
	
	madgraph->SetBinContent(i+1,madgraphVals[i][choice]);
	madgraph->SetBinError(i+1,0.0);
	mcatnlo->SetBinContent(i+1,mcatnloVals[i][choice]);
	mcatnlo->SetBinError(i+1,0.0);
	powheg->SetBinContent(i+1,powhegVals[i][choice]);
	powheg->SetBinError(i+1,0.0);
	}
 
 
	//do the plots  
 	TCanvas *c= new TCanvas("c","c",10,10,800,600);
  
  	muon_part->SetMarkerStyle(20);	 
  	muon_part->Draw("E");
	muon_part->SetMaximum(muon_part->GetBinContent(muon_part->GetMaximumBin())*1.5);
	muon_part->SetMinimum(0.);
	
	
	mcatnlo->Draw("same");
	powheg->Draw("same");
	madgraph->Draw("same");
	TLegend *tleg;
	tleg = new TLegend(0.6,0.75,0.85,0.9);
	tleg->SetTextSize(0.03);
	tleg->SetBorderSize(0);
	tleg->SetFillColor(10);

	tleg->AddEntry(muon_part  , "data '12'"      , "lep"); 
	tleg->AddEntry(madgraph  , "MadGraph"      , "l");
	tleg->AddEntry(mcatnlo  , "MC@NLO"      , "l");
	tleg->AddEntry(powheg  , "POWHEG"      , "l");
	tleg->Draw();
  	//titles
	muon_part->GetYaxis()->SetTitle("#partial #sigma [pb]"); muon_part->GetYaxis()->SetTitleSize(0.05);
	muon_part->GetXaxis()->SetTitle(Xtitle); muon_part->GetXaxis()->SetTitleSize(0.05);
   
   
       TText* textPrelim = doPrelim(0.16,0.96); 
       textPrelim->Draw();
       
       c->SaveAs("plots/Measurments/partialXsect"+Variable+".pdf");
	
	//normailise
	for(int i = 0; i < Nbins; i++){

	muon_part->SetBinContent(i+1,sigmaVal[i][choice]/totXsect[choice]);	
	muon_part->SetBinError(i+1,sigmaErr[i][choice]/totXsect[choice]);
//	muon_norm_diff->SetBinContent(i+1,sigmaVal[i]/(totXsect*width));
	
	madgraph->SetBinContent(i+1,madgraphVals[i][choice]/225.197);
	madgraph->SetBinError(i+1,0.0/225.197);
	mcatnlo->SetBinContent(i+1,mcatnloVals[i][choice]/225.197);
	mcatnlo->SetBinError(i+1,0.0/225.197);
	powheg->SetBinContent(i+1,powhegVals[i][choice]/225.197);
	powheg->SetBinError(i+1,0.0/225.197);
	
	}	
	
       
       TCanvas *c2= new TCanvas("c2","c2",10,10,800,600);
  
  	muon_part->SetMarkerStyle(20);
	 
  	muon_part->Draw("E");
	muon_part->SetMaximum(muon_part->GetBinContent(muon_part->GetMaximumBin())*1.5);
	
	
	mcatnlo->Draw("same");
	powheg->Draw("same");
	madgraph->Draw("same");
	
	tleg->Draw();
	
  	//titles
	muon_part->GetYaxis()->SetTitle("#frac{1}{#sigma} #partial #sigma [pb GeV^{-1}]"); 
	muon_part->GetYaxis()->SetTitleSize(0.05);
	//muon_part->GetXaxis()->SetTitle(Xtitle); 
	//muon_part->GetXaxis()->SetTitleSize(0.05);
   
       textPrelim->Draw();
       
       c2->SaveAs("plots/Measurments/partialXsectNorm"+Variable+".pdf");  
       
       	//normailise and differential
	for(int i = 0; i < Nbins; i++){
	double width = muon_part->GetBinWidth(i+1);
	
	muon_part->SetBinContent(i+1,sigmaVal[i][choice]/(totXsect[choice]*width));	
	muon_part->SetBinError(i+1,sigmaErr[i][choice]/(totXsect[choice]*width));
	
	madgraph->SetBinContent(i+1,madgraphVals[i][choice]/(225.197*width));
	madgraph->SetBinError(i+1,0.0/(225.197*width));
	mcatnlo->SetBinContent(i+1,mcatnloVals[i][choice]/(225.197*width));
	mcatnlo->SetBinError(i+1,0.0/(225.197*width));
	powheg->SetBinContent(i+1,powhegVals[i][choice]/(225.197*width));
	powheg->SetBinError(i+1,0.0/(225.197*width));
	}
       
       
        TCanvas *c3= new TCanvas("c3","c3",10,10,800,600);
  
  	muon_part->SetMarkerStyle(20);	 
  	muon_part->Draw("E");
	muon_part->SetMaximum(muon_part->GetBinContent(muon_part->GetMaximumBin())*1.5);
	
	
	mcatnlo->Draw("same");
	powheg->Draw("same");
	madgraph->Draw("same");
	
	tleg->Draw();
	
  	//titles
	muon_part->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{#partial #sigma}{#partial MET} 				[GeV^{-1}]"); 
	muon_part->GetYaxis()->SetTitleSize(0.05);
	//muon_part->GetXaxis()->SetTitle("MET [GeV]"); 
	//muon_part->GetXaxis()->SetTitleSize(0.05);
   
       textPrelim->Draw();
       
       c3->SaveAs("plots/Measurments/partialXsectNormDiff"+Variable+".pdf"); 


       
}



