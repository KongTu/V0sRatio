#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void printMassDist(){

	//TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/EPOS_TH3D_Sep10_5M_v1_vtx_2014.root");
	TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HIJING_TH3D_Sep9_15M_v1_vtx_2014.root");

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
	TH1D* xiHist_mass[6][20];
	
	TH1D* xiFracDenom[20];
	TH1D* xiFracNumer[20];
	
	TH1D* genks_mass[20];
	TH1D* genla_mass[20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};
    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};
    double binwidth_1[15] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.6,0.6,1.0,1.8,3.0};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    stringstream tempName1;
    stringstream tempName2;


    for ( int eta = 0; eta < 6; eta++){

	    for (int pt = 3; pt < 15; pt++){

	        ksHistName.str("");
	        laHistName.str("");
	        xiHistName.str("");
	        tempName1.str("");
	        tempName2.str("");

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

	        tempName1 << "temp1_";
	        tempName1 << eta+1;
	        tempName1 << "_";
	        tempName1 << pt;

	        tempName2 << "temp2_";
	        tempName2 << eta+1;
	        tempName2 << "_";
	        tempName2 << pt;

	        ks_mass[eta][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        la_mass[eta][pt] = laHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        xiHist_mass[eta][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

	        xiFracDenom[pt] = laHist->ProjectionZ( tempName1.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        xiFracNumer[pt] = xiHist->ProjectionZ( tempName2.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

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

	        genks_mass[pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        genla_mass[pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	   	
	   	}
   }


	double ksYield[6][20];
	double laYield[6][20];

	TCanvas* c1 = new TCanvas();

	c1->Print("laMass_hijing3bins.pdf[");

	for (eta = 0; eta < 6; eta++){

		for (pt = 3; pt < 6; pt++){

			la_mass[eta][pt]->Add( xiHist_mass[eta][pt], -1);
    		laYield[eta][pt] = la_YieldCal( la_mass[eta][pt] );

    		c1->Print("laMass_hijing3bins.pdf");



		}
	}

	c1->Print("laMass_hijing3bins.pdf]");



}