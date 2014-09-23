#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEfficiency(){


	TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HIJING_TH3D_Sep22_5M_v1_2014.root");
	//TFile* file = new TFile("~/Desktop/HIJING_TH3D_Sep16_10M_v1_ptDist_2014.root");
	
	TH3D* ksHist;
	TH3D* laHist;
	TH3D* genksHist;
	TH3D* genlaHist;

	TH3D* xiHist;

	ksHist = (TH3D*)file->Get("ana/InvMass_ks_underlying");
	laHist = (TH3D*)file->Get("ana/InvMass_la_underlying");
	genksHist = (TH3D*)file->Get("ana/genKS_underlying");
	genlaHist = (TH3D*)file->Get("ana/genLA_underlying");
	xiHist = (TH3D*)file->Get("ana/XiDaughter");

	TH1D* ks_mass[6][20];
	TH1D* la_mass[6][20];
	TH1D* genks_mass[6][20];
	TH1D* genla_mass[6][20];
	TH1D* xiHist_mass[6][20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for ( int eta = 0; eta < 6; eta++){

	    for (int pt = 3; pt < 15; pt++){

	        ksHistName.str("");
	        laHistName.str("");
	        xiHistName.str("");

	        ksHistName << "ks_";
	        ksHistName << eta+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "la_";
	        laHistName << eta+1;
	        laHistName << "_";
	        laHistName << pt;

	        xiHistName << "xiDau_";
	        xiHistName << eta+1;
	        xiHistName << "_";
	        xiHistName << pt;

	        ks_mass[eta][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        la_mass[eta][pt] = laHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        xiHist_mass[eta][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

	        ksHistName.str("");
	        laHistName.str("");

	        ksHistName << "genks_";
	        ksHistName << eta+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "genla_";
          	laHistName << eta+1;
	        laHistName << "_";
	        laHistName << pt;

	        genks_mass[eta][pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        genla_mass[eta][pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	   	
	   	}
   }

	double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
	double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
	double etabins[] = {-2.4,-1.6,-0.8,0.0,0.8,1.6,2.4};

   	int genksYield[6][20];
   	int genlaYield[6][20];

   	for ( eta = 0; eta < 6; eta++){
	   	
	   	for (int i = 3; i < 15; i++){

	   		genksYield[eta][i] = genks_mass[eta][i]->GetEntries();
	   		genlaYield[eta][i] = genla_mass[eta][i]->GetEntries();

	   	}
   	}	

   	double ksYield[6][20];
   	double laYield[6][20];

   	for ( eta = 0; eta < 6; eta++){

	   	for ( i = 3; i < 15; i++){

	   		ksYield[eta][i] = ks_YieldCal( ks_mass[eta][i] );
	   			la_mass[eta][i]->Add( xiHist_mass[eta][i],-1 );
	   		laYield[eta][i] = la_YieldCal( la_mass[eta][i] );
	   	}
   	}


/**
 * Efficiency in different eta bins:
 */
/*
 	TH1D* ks_eff_eta[6];
 	TH1D* la_eff_eta[6];
 	
 	stringstream ks_effName;
 	stringstream la_effName;

 	for (eta = 0; eta < 6; eta++){

 		ks_effName.str("");
 		ks_effName << "ks_eff_Eta_";
 		ks_effName << eta + 1;
 		ks_eff_eta[eta] = new TH1D(ks_effName.str().c_str(), ks_effName.str().c_str(),15,ptbins_1);

 		la_effName.str("");
 		la_effName << "la_eff_Eta_";
 		la_effName << eta + 1;
 		la_eff_eta[eta] = new TH1D(la_effName.str().c_str(), la_effName.str().c_str(),15,ptbins_1);

 		for (pt = 3; pt < 15; pt++){
 			
 			ks_eff_eta[eta]->SetBinContent(pt+1, ksYield[eta][pt]/genksYield[eta][pt]);
 				ks_eff_eta[eta]->SetBinError(pt+1, errorCal( ksYield[eta][pt], genksYield[eta][pt] ) );

 			la_eff_eta[eta]->SetBinContent(pt+1, laYield[eta][pt]/genlaYield[eta][pt]);
 				la_eff_eta[eta]->SetBinError(pt+1, errorCal( laYield[eta][pt], genlaYield[eta][pt] ) );
 		}
 	}
*/
/**
 * 2D efficiency Table:
 */


   	TH2D* ks_eff = new TH2D("ks_eff","ks_eff",6,etabins,15,ptbins_1);
   	TH2D* la_eff = new TH2D("la_eff","la_eff",6,etabins,15,ptbins_1);

   	//TH2D* ratio_eff = new TH2D("ratio_eff","ratio_eff",6,etabins,15,ptbins_1);

   	for ( eta = 0; eta < 6; eta++){

	   	for ( i = 3; i < 15; i++){

	   		double temp_ks = ksYield[eta][i]/genksYield[eta][i];
	   		double temp_la = laYield[eta][i]/genlaYield[eta][i];

	   		ks_eff->SetBinContent(eta+1,i+1, temp_ks );
	   			ks_eff->SetBinError(eta+1,i+1, errorCal(ksYield[eta][i], genksYield[eta][i]) );
	   		la_eff->SetBinContent(eta+1,i+1, temp_la );
	   			la_eff->SetBinError(eta+1,i+1, errorCal(laYield[eta][i], genlaYield[eta][i]) );

	   		//ratio_eff->SetBinContent(eta+1,i+1, temp_la/temp_ks );

	   	}
   	}


   	TFile f1("HIJING_withXiRemoval_effKongNew2DTable_5M_Sep22_v1_12pTbins.root","new");
   	ks_eff->Write();
   	la_eff->Write();


   	/*for ( eta = 0; eta < 6; eta++){

   		ks_eff_eta[eta]->Write();
   		la_eff_eta[eta]->Write();
   	}*/
   

}