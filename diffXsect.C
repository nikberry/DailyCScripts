#include "etaFit.C"
#include "TString.h"
#include "TLatex.h"
#include <vector>

void diffXsect();

void diffXsect(){

//MET will need choice of variable at the top
TString Variable ="_MET";
int Nbins = 6;
TString bins[6] = {"0-25", "25-45", "45-70", "70-100", "100-150", "150-inf"};
double width[6] = {25, 20, 25, 30, 50, 100};
double xbins[7] = {1,25,45,70,100,150, 250}; 
TString varBin = "Binned_MET_Analysis/RecoMET_bin_";

//HT
//TString Variable ="_HT";
// int Nbins = 8;
// TString bins[8] = {"0-50", "50-150", "150-250", "250-350", "350-450", "450-650", "650-1100", "1100-inf"};
// double width[8] = {50,100,100,100,100,200,450,400};
// double xbins[9] = {1,50,150,250,350,450,650,1100, 1500}; 
// TString varBin = "Binned_HT_Analysis/HT_bin_";

//ST
//TString Variable ="_ST";
// int Nbins = 8;
// TString bins[8] = {"0-150", "150-250", "250-350", "350-450", "450-550", "550-750", "750-1250", "1250-inf"};
// double width[8] = {150,100,100,100,100,200,500,500};
// double xbins[9] = {1,150,250,350,450,550,750,1250, 1750}; 
// TString varBin = "Binned_ST_Analysis/ST_with_RecoMET_bin_";

//MT
//TString Variable ="_MT";
// int Nbins = 5;
// TString bins[5] = {"0-40", "40-65", "65-85", "85-150", "150-inf"};
// double width[5] = {40,25,20,65,50};
// double xbins[6] = {1,40,65,85,150,200}; 
// TString varBin = "Binned_MT_Analysis/MT_with_RecoMET_bin_";

double sigmaVal[Nbins];
double sigmaErr[Nbins];
double madgraphVals[Nbins];

double totXsect = 0;
//sample
//TString dir = "central";
//TString dir = "JES_up";
TString dirs[3] = {"central","JES_up","JES_down"};

//loop over systematics
for(int sys = 0; sys < 3; sys++){
TString dir  = dirs[sys];

//loop over bins of distribution
for(int i = 0; i < Nbins; i++){

TString bin = varBin;
bin += bins[i];

cout <<  bin << endl;

sigmaVal[i] = etaFit(bin, "measured",dir);
sigmaErr[i] = etaFit(bin, "measuredErr",dir);

madgraphVals[i] = etaFit(bin, "madgraph",dir);

totXsect += sigmaVal[i];
}

cout << "partial xsects: " << endl;
for(int i = 0; i < Nbins; i++){
cout << bins[i] << ": " << " = " << sigmaVal[i] << endl;
}

cout << "normalised xsects: " << endl;
for(int i = 0; i < Nbins; i++){
cout << bins[i] << ": " << " = " << sigmaVal[i]/totXsect << endl;
}

cout << "normalised differential: " << endl;
for(int i = 0; i < Nbins; i++){
cout << bins[i] << ": " << " = " << sigmaVal[i]/(totXsect*width[i]) << endl;
}
cout << "cross section is:  " <<  totXsect << endl;  

   //measured histo will have to put name of systematic in here to write into file
   TH1D *muon_part  = new TH1D(dir, "", Nbins, xbins);  //muon
   
   //different generators
   TH1D *madgraph  = new TH1D("madgraph", "", Nbins, xbins); 
	madgraph->SetLineColor(kRed);
	
	for(int i = 0; i < Nbins; i++){
	muon_part->SetBinContent(i+1,sigmaVal[i]);	
	muon_part->SetBinError(i+1, sigmaErr[i]);
	
	madgraph->SetBinContent(i+1,madgraphVals[i]);
	madgraph->SetBinError(i+1,0.0);
	}
 
 
	//do the plots  
 	TCanvas *c= new TCanvas("c","c",10,10,800,600);
  
  	muon_part->SetMarkerStyle(20);	 
  	muon_part->Draw("E");
	muon_part->SetMaximum(muon_part->GetBinContent(muon_part->GetMaximumBin())*1.5);
	muon_part->SetMinimum(0.);
	
	madgraph->Draw("same");
	
	TLegend *tleg;
	tleg = new TLegend(0.65,0.75,0.9,0.9);
	tleg->SetTextSize(0.03);
	tleg->SetBorderSize(0);
	tleg->SetFillColor(10);

	tleg->AddEntry(muon_part  , "data '12'"      , "lep"); 
	tleg->AddEntry(madgraph  , "MadGraph"      , "l");
	tleg->Draw();
  	//titles
	muon_part->GetYaxis()->SetTitle("#partial #sigma [pb]"); muon_part->GetYaxis()->SetTitleSize(0.05);
	muon_part->GetXaxis()->SetTitle("MET [GeV]"); muon_part->GetXaxis()->SetTitleSize(0.05);
   
   
       TText* textPrelim = doPrelim(0.16,0.96); 
       textPrelim->Draw();
       
       c->SaveAs("plots/Measurments/partialXsect"+Variable+".png");
	
	//normailise
	for(int i = 0; i < Nbins; i++){

	muon_part->SetBinContent(i+1,sigmaVal[i]/totXsect);	
	muon_part->SetBinError(i+1,sigmaErr[i]/totXsect);
//	muon_norm_diff->SetBinContent(i+1,sigmaVal[i]/(totXsect*width));
	
	madgraph->SetBinContent(i+1,madgraphVals[i]/225.2);
	madgraph->SetBinError(i+1,0.0/225.2);
	}	
	
       
       TCanvas *c2= new TCanvas("c2","c2",10,10,800,600);
  
  	muon_part->SetMarkerStyle(20);
	 
  	muon_part->Draw("E");
	muon_part->SetMaximum(muon_part->GetBinContent(muon_part->GetMaximumBin())*1.5);
	
	madgraph->Draw("same");
	
	tleg->Draw();
	
  	//titles
	muon_part->GetYaxis()->SetTitle("#frac{1}{#sigma} #partial #sigma [pb GeV^{-1}]"); 
	muon_part->GetYaxis()->SetTitleSize(0.05);
	muon_part->GetXaxis()->SetTitle("MET [GeV]"); 
	muon_part->GetXaxis()->SetTitleSize(0.05);
   
       textPrelim->Draw();
       
       c2->SaveAs("plots/Measurments/partialXsectNorm"+Variable+".png");  
       
       	//normailise and differential
	for(int i = 0; i < Nbins; i++){
	double width = muon_part->GetBinWidth(i+1);
	
	muon_part->SetBinContent(i+1,sigmaVal[i]/(totXsect*width));	
	muon_part->SetBinError(i+1,sigmaErr[i]/(totXsect*width));
	
	madgraph->SetBinContent(i+1,madgraphVals[i]/(225.2*width));
	madgraph->SetBinError(i+1,0.0/(225.2*width));
	}
       
       
        TCanvas *c3= new TCanvas("c3","c3",10,10,800,600);
  
  	muon_part->SetMarkerStyle(20);	 
  	muon_part->Draw("E");
	muon_part->SetMaximum(muon_part->GetBinContent(muon_part->GetMaximumBin())*1.5);
	
	madgraph->Draw("same");
	
	tleg->Draw();
	
  	//titles
	muon_part->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{#partial #sigma}{#partial MET} 				[GeV^{-1}]"); 
	muon_part->GetYaxis()->SetTitleSize(0.05);
	muon_part->GetXaxis()->SetTitle("MET [GeV]"); 
	muon_part->GetXaxis()->SetTitleSize(0.05);
   
       textPrelim->Draw();
       
       c3->SaveAs("plots/Measurments/partialXsectNormDiff"+Variable+".png"); 
       
       //will need to change MET for other variables
       TFile resultsfile("outFiles/diffResults"+Variable+".root", "RECREATE", "comment");

       muon_part->Write();
       resultsfile.Write();
       resultsfile.Close();

}//end loop over systematic
       
}



