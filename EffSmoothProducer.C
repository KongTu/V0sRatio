#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include <stdio.h>
#include <iomanip>
#include <vector>
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;

void EffSmoothProducer(){

	gStyle->SetErrorX(0);

	double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_ptbincenter[28] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
  	double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Nov3_rapidity_v2_28ks_pTbins.root");

	TH1D* ks_eff_rpy[5];
	TH1D* la_eff_rpy[5];

	stringstream ksName;
	stringstream laName;

	for (int rpy = 0; rpy < 5; rpy++){

		ksName.str("");
		laName.str("");

		ksName << "ks_eff_rpy_";
		ksName << rpy+1;

		laName << "la_eff_rpy_";
		laName << rpy+1;

		ks_eff_rpy[rpy] = (TH1D*) file->Get( ksName.str().c_str() );
		la_eff_rpy[rpy] = (TH1D*) file->Get( laName.str().c_str() );

	}


	TH1D* ks_eff_rpy_new[5];
	TH1D* la_eff_rpy_new[5];

	for(rpy = 0; rpy < 5; rpy++){

		ksName.str("");
		laName.str("");

		ksName << "ks_eff_rpy_new_";
		ksName << rpy+1;

		laName << "la_eff_rpy_new_";
		laName << rpy+1;

		ks_eff_rpy_new[rpy] = new TH1D(ksName.str().c_str(),ksName.str().c_str(),28,ks_ptbins);
		la_eff_rpy_new[rpy] = new TH1D(laName.str().c_str(),laName.str().c_str(),20,la_ptbins);

		/*
		smoothing the efficiency and only store the values after 2GeV:
		 */
		
		ks_eff_rpy[rpy]->Smooth(10);
		la_eff_rpy[rpy]->Smooth(10);

			for(pt = 0; pt < 28; pt++){

				ks_eff_rpy_new[rpy]->SetBinContent(pt+1, ks_eff_rpy[rpy]->GetBinContent(pt+1) );
					ks_eff_rpy_new[rpy]->SetBinError(pt+1, ks_eff_rpy[rpy]->GetBinError(pt+1) );

			}

			for(pt = 0; pt < 20; pt++){

				la_eff_rpy_new[rpy]->SetBinContent(pt+1, la_eff_rpy[rpy]->GetBinContent(pt+1) );
					la_eff_rpy_new[rpy]->SetBinError(pt+1, la_eff_rpy[rpy]->GetBinError(pt+1) );
			}

		
	}

	TFile f1("EffSmoothProducer_v2.root","new");	

	for(rpy = 0; rpy < 5; rpy++){

		ks_eff_rpy_new[rpy]->Write();
		la_eff_rpy_new[rpy]->Write();
	}
		




}