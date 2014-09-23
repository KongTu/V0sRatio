#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEfficiencyRapidity(){


	TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_rpyDependent_3Dhisto/HIJING_TH3D_Sep22_10M_rpyDependent_v1_2014.root");
	//TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/EPOS_TH3D_Sep10_5M_v1_vtx_2014.root");
	
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

	TH1D* ks_mass[5][20];
	TH1D* la_mass[5][20];
	TH1D* genks_mass[5][20];
	TH1D* genla_mass[5][20];
	TH1D* xiHist_mass[5][20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};

    double rpybins[6] = {33,91,145,199,250,314};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for ( int rpy = 0; rpy < 5; rpy++){

	    for (int pt = 3; pt < 15; pt++){

	        ksHistName.str("");
	        laHistName.str("");
	        xiHistName.str("");

	        ksHistName << "ks_";
	        ksHistName << rpy+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "la_";
	        laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt;

	        xiHistName << "xiDau_";
	        xiHistName << rpy+1;
	        xiHistName << "_";
	        xiHistName << pt;

	        ks_mass[rpy][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        la_mass[rpy][pt] = laHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        xiHist_mass[rpy][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

	        ksHistName.str("");
	        laHistName.str("");

	        ksHistName << "genks_";
	        ksHistName << rpy+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "genla_";
          	laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt;

	        genks_mass[rpy][pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        genla_mass[rpy][pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	   	
	   	}
   }

	double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
	double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
	//double etabins[] = {-2.4,-1.6,-0.8,0.0,0.8,1.6,2.4};

   	int genksYield[5][20];
   	int genlaYield[5][20];

   	for ( rpy = 0; rpy < 5; rpy++){
	   	
	   	for (int i = 3; i < 15; i++){

	   		genksYield[rpy][i] = genks_mass[rpy][i]->GetEntries();
	   		genlaYield[rpy][i] = genla_mass[rpy][i]->GetEntries();

	  

	   	}
   	}	

   	double ksYield[5][20];
   	double laYield[5][20];

   	for ( rpy = 0; rpy < 5; rpy++){

	   	for ( i = 3; i < 15; i++){

	   		ksYield[rpy][i] = ks_YieldCal( ks_mass[rpy][i] );
	   			//la_mass[rpy][i]->Add( xiHist_mass[rpy][i],-1 );
	   		laYield[rpy][i] = la_YieldCal( la_mass[rpy][i] );
	   	}
   	}

   	for ( rpy = 0; rpy < 5; rpy++){

   		for ( i = 3; i < 15; i++){

			cout << "ks___rpy = " << rpy << " and " << "pt = " << i << ",yield: " << ksYield[rpy][i] << endl;
			cout << "la___rpy = " << rpy << " and " << "pt = " << i <<  ",yield: " << laYield[rpy][i] << endl;

			cout << "genks___rpy = " << rpy << " and " << "pt = " << i << ",genksyield: " << genksYield[rpy][i] << endl;
			cout << "genla___rpy = " << rpy << " and " << "pt = " << i <<  ",genlayield: " << genlaYield[rpy][i] << endl;
   			
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


   	TH2D* ks_eff = new TH2D("ks_eff","ks_eff",5,rpybins,15,ptbins_1);
   	TH2D* la_eff = new TH2D("la_eff","la_eff",5,rpybins,15,ptbins_1);

   	//TH2D* ratio_eff = new TH2D("ratio_eff","ratio_eff",6,etabins,15,ptbins_1);

   	for ( rpy = 0; rpy < 5; rpy++){

	   	for ( i = 3; i < 15; i++){

	   		double temp_ks = ksYield[rpy][i]/genksYield[rpy][i];
	   		double temp_la = laYield[rpy][i]/genlaYield[rpy][i];

	   		ks_eff->SetBinContent(rpy+1,i+1, temp_ks );
	   			ks_eff->SetBinError(rpy+1,i+1, errorCal(ksYield[rpy][i], genksYield[rpy][i]) );
	   		la_eff->SetBinContent(rpy+1,i+1, temp_la );
	   			la_eff->SetBinError(rpy+1,i+1, errorCal(laYield[rpy][i], genlaYield[rpy][i]) );

	   		//ratio_eff->SetBinContent(eta+1,i+1, temp_la/temp_ks );

	   	}
   	}



   	TFile f1("effKongNew2DTable_withoutXiRemoval_10M_Sep22_rapidity_v1_12pTbins.root","new");
   	ks_eff->Write();
   	la_eff->Write();
 

   	/*for ( eta = 0; eta < 6; eta++){

   		ks_eff_eta[eta]->Write();
   		la_eff_eta[eta]->Write();
   	}*/

   	//ratio_eff->Write();
   	//
   

}