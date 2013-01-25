#include "etaFit.C"
#include "TString.h"
#include "TLatex.h"

void diffXsect();

void diffXsect(){

int Nbins = 5;

//MET will need choice of variable at the top
double sigmaVal[Nbins];
double sigmaErr[Nbins];
double madgraphVals[Nbins];

TString bins[5] = {"0-25", "25-45", "45-70", "70-100", "100-inf"};
double width[5] = {25, 20, 25, 30, 50};
double xbins[6]= {1,25,45,70,100,150}; 

//TString Variables[N] = {"DeltaPhi_lepton_MET_", "Transverse_Mass_", "METsignificance_", "MET_", "MET_phi_"};
double totXsect = 0;


for(int i = 0; i < Nbins; i++){

TString bin = "BinnedMETAnalysis/Muon_RecoMET_bin_";

bin += bins[i];
sigmaVal[i] = etaFit(bin, "measured");
sigmaErr[i] = etaFit(bin, "measuredErr");

madgraphVals[i] = etaFit(bin, "madgraph");

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

   //measured histo
   TH1D *muon_part  = new TH1D("muon part", "", 5, xbins);  //muon
   
   //different generators
   TH1D *madgraph  = new TH1D("madgraph", "", 5, xbins); 
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
       
       c->SaveAs("plots/Measurments/partialXsect.png");
	
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
       
       c2->SaveAs("plots/Measurments/partialXsectNorm.png");  
       
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
       
       c3->SaveAs("plots/Measurments/partialXsectNormDiff.png"); 
       
       
       
}



