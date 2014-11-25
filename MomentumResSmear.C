#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void MomentumResSmear(){

	TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/newEPOSsample/EPOS_TH3D_Nov21_vtx_2014.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/ROOT_file/newEPOSsample/EPOS_TH3D_Nov21_vtx_pTsmear_2014.root");

	TH3D* genKS = (TH3D*)file1->Get("ana/genKS_underlying");
	TH3D* genLA = (TH3D*)file1->Get("ana/genLA_underlying");

	TH3D* genKS_smear = (TH3D*)file2->Get("ana/genKS_underlying");
	TH3D* genLA_smear = (TH3D*)file2->Get("ana/genLA_underlying");

	//cout << "number of entries: " << genKS->GetEntries() << endl;
	//cout << "number of entries: " << genKS_smear->GetEntries() << endl;


	TH1D* genks_mass[5][28];
	TH1D* genla_mass[5][20];

	TH1D* genks_smear_mass[5][28];
	TH1D* genla_smear_mass[5][28];



	double ks_pTbinsBound[29] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_binwidth[28] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};
    double ks_ptbincenter[28] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

    double la_pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
    double la_binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};

    double rpybins[6] = {6,16,26,35,44,55};
    double rpybinwidth[5] = {1.0,0.97,0.9,0.93,1.0};

	stringstream ksName;
	stringstream laName;

	for(int rpy = 0; rpy < 5; rpy++){

		for(int pt = 0; pt < 28; pt++){

			ksName.str("");
			ksName << "genks_mass";
			ksName << "_";
			ksName << rpy+1;
			ksName << "_";
			ksName << pt+1;

			genks_mass[rpy][pt] = genKS->ProjectionZ(ksName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1]);

			ksName.str("");
			ksName << "genks_smear_mass";
			ksName << "_";
			ksName << rpy+1;
			ksName << "_";
			ksName << pt+1;

			genks_smear_mass[rpy][pt] = genKS_smear->ProjectionZ(ksName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1]);
		}

		for(pt = 0; pt < 20; pt++){

			laName.str("");
			laName << "genla_mass";
			laName << "_";
			laName << rpy+1;
			laName << "_";
			laName << pt+1;

			genla_mass[rpy][pt] = genLA->ProjectionZ(laName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1]);

			laName.str("");
			laName << "genla_smear_mass";
			laName << "_";
			laName << rpy+1;
			laName << "_";
			laName << pt+1;

			genla_smear_mass[rpy][pt] = genLA_smear->ProjectionZ(laName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1]);


		}
	}

   	int genksYield[5][28];
   	int genlaYield[5][20];

  	int genksSmearYield[5][28];
  	int genlaSmearYield[5][20];

   	for(rpy = 0; rpy < 5; rpy++){

   		for( pt = 0; pt < 28; pt++){

   			genksYield[rpy][pt] = genks_mass[rpy][pt]->GetEntries();
   			genksSmearYield[rpy][pt] = genks_smear_mass[rpy][pt]->GetEntries();
   		}
	   	
	   	for (pt = 0; pt < 20; pt++){

	   		genlaYield[rpy][pt] = genla_mass[rpy][pt]->GetEntries();
	   		genlaSmearYield[rpy][pt] = genla_smear_mass[rpy][pt]->GetEntries();

	   	}
	}
	   		
	

/*
GEN spectra distribution;
 */

	TH1D* ksSpectra_gen[5];
	TH1D* laSpectra_gen[5];

	TH1D* ksSpectra_gen_smear[5];
	TH1D* laSpectra_gen_smear[5];

	for(rpy = 0; rpy < 5; rpy++){

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_epos_gen_";
		ksName << rpy+1;

		laName << "laSpectra_epos_gen_";
		laName << rpy+1;

		ksSpectra_gen[rpy] = new TH1D(ksName.str().c_str(),ksName.str().c_str(),28,ks_ptbins);
		laSpectra_gen[rpy] = new TH1D(laName.str().c_str(),laName.str().c_str(),20,la_ptbins);

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_epos_gen_smear_";
		ksName << rpy+1;

		laName << "laSpectra_epos_gen_smear_";
		laName << rpy+1;

		ksSpectra_gen_smear[rpy] = new TH1D(ksName.str().c_str(),ksName.str().c_str(),28,ks_ptbins);
		laSpectra_gen_smear[rpy] = new TH1D(laName.str().c_str(),laName.str().c_str(),20,la_ptbins);

		for(pt = 0; pt < 28; pt++){

			ksSpectra_gen[rpy]->SetBinContent(pt+1, genksYield[rpy][pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));
				ksSpectra_gen[rpy]->SetBinError(pt+1, sqrt( genksYield[rpy][pt] )/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));

			ksSpectra_gen_smear[rpy]->SetBinContent(pt+1, genksSmearYield[rpy][pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));
				ksSpectra_gen_smear[rpy]->SetBinError(pt+1, sqrt( genksSmearYield[rpy][pt] )/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));

		}

		for(pt = 0; pt < 20; pt++){

			laSpectra_gen[rpy]->SetBinContent(pt+1, genlaYield[rpy][pt]/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy] ));
				laSpectra_gen[rpy]->SetBinError(pt+1, sqrt( genlaYield[rpy][pt] )/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy]));

			laSpectra_gen_smear[rpy]->SetBinContent(pt+1, genlaSmearYield[rpy][pt]/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy] ));
				laSpectra_gen_smear[rpy]->SetBinError(pt+1, sqrt( genlaSmearYield[rpy][pt] )/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy] ));

		}
	}

	TLine* l1 = new TLine(0,1,9.0,1.0);
	l1->SetLineWidth(2);
	l1->SetLineColor(kRed);
	l1->SetLineStyle(2);

	TLatex* ratio[5];
    ratio[0] = new TLatex(1.65,1.06,"-2.87 < y < -1.8");
    ratio[1] = new TLatex(1.65,1.06,"-1.8 < y < -0.9");
    ratio[2] = new TLatex(1.65,1.06,"-0.9 < y < 0");
    ratio[3] = new TLatex(1.65,1.06,"0 < y < 0.93");
    ratio[4] = new TLatex(1.65,1.06,"0.93 < y < 1.93");

	TCanvas* c1 = new TCanvas();
	c1->Divide(3,2,0,0);
	for(rpy = 0; rpy < 5; rpy++){

		c1->cd(rpy+1);
		laSpectra_gen[rpy]->Divide( laSpectra_gen_smear[rpy] );
		laSpectra_gen[rpy]->SetMarkerStyle(20);
		laSpectra_gen[rpy]->SetLineColor(kBlack);
		laSpectra_gen[rpy]->GetYaxis()->SetRangeUser(0.9,1.1);
		laSpectra_gen[rpy]->SetTitle("#Lambda/#bar{#Lambda}");
		laSpectra_gen[rpy]->SetYTitle("noSmear/Smear");
		laSpectra_gen[rpy]->SetXTitle("pT(GeV/c)");
		laSpectra_gen[rpy]->SetStats(kFALSE);

		laSpectra_gen[rpy]->Draw("P");
		l1->Draw("same");
		ratio[rpy]->Draw("same");


	}


}