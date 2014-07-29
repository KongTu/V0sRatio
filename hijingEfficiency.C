#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEfficiency(){


	TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HIJING_TH3D_July28_2014.root");

	TH3D* ksHist;
	TH3D* laHist;
	TH3D* genksHist;
	TH3D* genlaHist;

	ksHist = (TH3D*)file->Get("ana/InvMass_ks_underlying");
	laHist = (TH3D*)file->Get("ana/InvMass_la_underlying");
	genksHist = (TH3D*)file->Get("ana/genKS_underlying");
	genlaHist = (TH3D*)file->Get("ana/genLA_underlying");

	TH1D* ks_mass[6][20];
	TH1D* la_mass[6][20];
	TH1D* genks_mass[6][20];
	TH1D* genla_mass[6][20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};

    stringstream ksHistName;
    stringstream laHistName;

    for ( int eta = 0; eta < 6; eta++){

	    for (int pt = 0; pt < 20; pt++){

	        ksHistName.str("");
	        laHistName.str("");

	        ksHistName << "ks1_";
	        ksHistName << pt;

	        laHistName << "la1_";
	        laHistName << pt;

	        ks_mass[eta][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );
	        la_mass[eta][pt] = laHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );

	        ksHistName.str("");
	        laHistName.str("");

	        ksHistName << "genks1_";
	        ksHistName << pt;

	        laHistName << "genla1_";
	        laHistName << pt;

	        genks_mass[eta][pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );
	        genla_mass[eta][pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );
	   	
	   	}
   }

   	int genksYield[6][20];
   	int genlaYield[6][20];

   	for ( eta = 0; eta < 6; eta++){
	   	
	   	for (int i = 0; i < 20; i++){

	   		genksYield[eta][i] = genks_mass[eta][i]->GetEntries();
	   		genlaYield[eta][i] = genla_mass[eta][i]->GetEntries();

	   	}
   	}	

   	double ksYield[6][20];
   	double laYield[6][20];

   	for ( eta = 0; eta < 6; eta++){

	   	for ( i = 0; i < 20; i++){

	   		ksYield[eta][i] = ks_YieldCal( ks_mass[eta][i] );
	   		laYield[eta][i] = la_YieldCal( la_mass[eta][i] );
	   	}
   	}

   	double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
   	double etabins[] = {-2.4,-1.6,-0.8,0.0,0.8,1.6,2.4};

   	TH2D* ks_eff = new TH2D("ks_eff","ks_eff",6,etabins,20,ptbins);
   	TH2D* la_eff = new TH2D("la_eff","la_eff",6,etabins,20,ptbins);

   	TH2D* ratio_eff = new TH2D("ratio_eff","ratio_eff",6,etabins,20,ptbins);

   	for ( eta = 0; eta < 6; eta++){

	   	for ( i = 0; i < 20; i++){

	   		double temp_ks = ksYield[eta][i]/genksYield[eta][i];
	   		double temp_la = laYield[eta][i]/genlaYield[eta][i];

	   		ks_eff->SetBinContent(eta+1,i+1, temp_ks );
	   		la_eff->SetBinContent(eta+1,i+1, temp_la );

	   		ratio_eff->SetBinContent(eta+1,i+1, temp_ks/temp_la );

	   	}
   	}


   	TFile f1("eff_2Dtable.root","new");
   	ks_eff->Write();
   	la_eff->Write();
   	ratio_eff->Write();


   	TCanvas* c1 = new TCanvas();
   	ks_eff->Draw("lego2");

   	TCanvas* c2 = new TCanvas();
   	la_eff->Draw("lego2");

   	TCanvas* c3 = new TCanvas();
   	ratio_eff->Draw("lego2");







}