#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hiRegitEfficiency(){


	TFile* file = new TFile("~/Desktop/RegitPYTHIA80_Oct3_inJet30Cone_2014.root");
	
	TH3D* ksHist;
	TH3D* laHist;
	TH3D* genksHist;
	TH3D* genlaHist;

	TH3D* xiHist;

	ksHist = (TH3D*)file->Get("v0analyzerHI/InvMass_ks_underlying");
	laHist = (TH3D*)file->Get("v0analyzerHI/InvMass_la_underlying");
	genksHist = (TH3D*)file->Get("v0analyzerHI/genKS_underlying");
	genlaHist = (TH3D*)file->Get("v0analyzerHI/genLA_underlying");
	xiHist = (TH3D*)file->Get("v0analyzerHI/XiDaughter");

	TH1D* ks_mass[20];
	TH1D* la_mass[20];
	TH1D* genks_mass[20];
	TH1D* genla_mass[20];
	TH1D* xiHist_mass[20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for (int pt = 3; pt < 15; pt++){

        ksHistName.str("");
        laHistName.str("");
        xiHistName.str("");

        ksHistName << "ks_";
        ksHistName << pt;

        laHistName << "la_";
        laHistName << pt;

        xiHistName << "xiDau_";
        xiHistName << pt;

        ks_mass[pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
        la_mass[pt] = laHist->ProjectionZ( laHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
        //xiHist_mass[pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

        ksHistName.str("");
        laHistName.str("");

        ksHistName << "genks_";
        ksHistName << pt;

        laHistName << "genla_";
        laHistName << pt;

        genks_mass[pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
        genla_mass[pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
   	
   	}
 

	double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
	double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
	double etabins[] = {-2.4,-1.6,-0.8,0.0,0.8,1.6,2.4};

   	int genksYield[20];
   	int genlaYield[20];

   
   	
   	for (int i = 3; i < 15; i++){

   		genksYield[i] = genks_mass[i]->GetEntries();
   		genlaYield[i] = genla_mass[i]->GetEntries();

   	}


   	double ksYield[20];
   	double laYield[20];



   	for ( i = 3; i < 15; i++){

   		ksYield[i] = ks_YieldCal( ks_mass[i] );
   			//la_mass[i]->Add( xiHist_mass[i],-1 );
   		laYield[i] = la_YieldCal( la_mass[i] );
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


   	TH1D* ks_eff = new TH1D("ks_eff","ks_eff",15,ptbins_1);
   	TH1D* la_eff = new TH1D("la_eff","la_eff",15,ptbins_1);

   	for ( i = 3; i < 15; i++){

   		double temp_ks = ksYield[i]/genksYield[i];
   		double temp_la = laYield[i]/genlaYield[i];

   		ks_eff->SetBinContent(i+1, temp_ks );
   			ks_eff->SetBinError(i+1, errorCal(ksYield[i], genksYield[i]) );
   		la_eff->SetBinContent(i+1, temp_la );
   			la_eff->SetBinError(i+1, errorCal(laYield[i], genlaYield[i]) );

   	}
   	

   	TFile f1("hiRegitTracking_HI_12pTbins.root","new");
   	ks_eff->Write();
   	la_eff->Write();


   	/*for ( eta = 0; eta < 6; eta++){

   		ks_eff_eta[eta]->Write();
   		la_eff_eta[eta]->Write();
   	}*/
   

}