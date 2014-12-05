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

void testSmooth(){
/*
with partially smoothing efficiency, only smoothing after 2GeV:
 */

    TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Nov3_rapidity_v2_28ks_pTbins.root");
    TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/EffSmoothProducer/EffSmoothProducer_v1.root");

    TH1D* ks_eff_rpy[5];
    TH1D* la_eff_rpy[5];

    TH1D* ks_eff_rpy_new[5];
    TH1D* la_eff_rpy_new[5];

    double ks_eff[5][28];
    double la_eff[5][20];

    double ks_eff_err[5][28];
    double la_eff_err[5][20];

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

        ksName.str("");
        laName.str("");

        ksName << "ks_eff_rpy_new_";
        ksName << rpy+1;

        laName << "la_eff_rpy_new_";
        laName << rpy+1;

        ks_eff_rpy_new[rpy] = (TH1D*) file1->Get( ksName.str().c_str() );
        la_eff_rpy_new[rpy] = (TH1D*) file1->Get( laName.str().c_str() );

    }

    for(rpy = 0; rpy < 5; rpy++){

        for(int pt = 0; pt < 15; pt++){

            ks_eff[rpy][pt] = ks_eff_rpy[rpy]->GetBinContent( pt+1 );

                ks_eff_err[rpy][pt] = ks_eff_rpy[rpy]->GetBinError( pt+1 );
        }

        for(pt = 0; pt < 7; pt++){

            la_eff[rpy][pt] = la_eff_rpy[rpy]->GetBinContent( pt+1 );

                la_eff_err[rpy][pt] = la_eff_rpy[rpy]->GetBinError( pt+1 );  
        }

        /*
        adding smoothing after 2GeV:
         */
        
        for(pt = 15; pt < 28; pt++){

            ks_eff[rpy][pt] = ks_eff_rpy_new[rpy]->GetBinContent( pt+1 );

                ks_eff_err[rpy][pt] = ks_eff_rpy_new[rpy]->GetBinError( pt+1 );
        }

        for(pt = 7; pt < 20; pt++){

            la_eff[rpy][pt] = la_eff_rpy_new[rpy]->GetBinContent( pt+1 );

                la_eff_err[rpy][pt] = la_eff_rpy_new[rpy]->GetBinError( pt+1 );  
        }


    }
    double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    
    TH1D* temp = new TH1D("temp","temp",20,la_ptbins);

    for(pt = 0; pt < 20; pt++){

        temp->SetBinContent(pt+1,la_eff[0][pt]);
    }

    temp->Divide( la_eff_rpy[0] );
    temp->SetMarkerStyle(20);
    temp->Draw("P");




}