#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include <vector>
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include <sstream>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHist.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

float Mass_ks(float px_1,float py_1,float pz_1,float px_2,float py_2,float pz_2)
{
       
	float temp = 0.0;
        float E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.13957*0.13957));
        float E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.93827*0.93827));
        float E_tot = E1+E2;
	temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
	return sqrt(temp);
}

float Mass_la(float px_1,float py_1,float pz_1,float px_2,float py_2,float pz_2)
{
       
  float temp = 0.0;
        float E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.13957*0.13957));
        float E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.13957*0.13957));
        float E_tot = E1+E2;
  temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
  return sqrt(temp);
}

float Mass_e(float px_1,float py_1,float pz_1,float px_2,float py_2,float pz_2)
{
        float temp = 0.0;
        float E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.000511*0.000511));
        float E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.000511*0.000511));
        float E_tot = E1+E2;
        temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
        return sqrt(temp);
}

double errorCal(double yield_la, double yield_ks){

    double ks_temp = sqrt( yield_ks );
    double la_temp = sqrt( yield_la );

    double first = (la_temp*la_temp)/((yield_ks)*(yield_ks));
    double second = ( (yield_la)*(yield_la)*ks_temp*ks_temp )/((yield_ks)*(yield_ks)*(yield_ks)*(yield_ks));

    double error = sqrt( first + second );
    return error;

}

double ks_YieldCal( TH1D* inputHist ){

//Starting to do RooFit:

    RooRealVar x("x","mass",0.44,0.56);
    RooDataHist data("data","dataset",x, inputHist );
    RooPlot* xframe = x.frame(240);

    data.plotOn(xframe, Name("data"));

    RooRealVar mean("mean","mean",0.50,0.49,0.51);
    RooRealVar sigma1("sigma1","sigma1",0.003,0.001,0.01);
    RooRealVar sigma2("sigma2","sigma2",0.003,0.001,0.01);
    RooRealVar sig1("sig1","signal1",10,0,10000000);
    RooRealVar sig2("sig2","signal2",10,0,10000000);
    RooRealVar a("a","a",0,-100000,100000);
    RooRealVar b("b","b",0,-100000,100000);
    RooRealVar cp("cp","cp",0,-100000,100000);
    RooRealVar d("d","d",0,-100000,100000);

    RooRealVar f("f","f",0,-100000,100000);
    RooRealVar g("g","g",0,-100000,100000);
    RooRealVar h("h","h",0,-100000,100000);
    RooRealVar k("k","k",0,-100000,100000);

    RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
    RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
    RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
    RooRealVar polysig("polysig","polysig",10,0,10000000);
    RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

    x.setRange("cut",0.45,0.54);

    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));

    sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
    sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));

    xframe->Draw();

    double chi2  = xframe->chiSquare("sum","data");
    double meanf  = mean.getVal();
    double meanfe  = mean.getError();
    double sigmaf1  = sigma1.getVal();
    double sigmaf2  = sigma2.getVal();
    double bkgf  = polysig.getVal();
    double sigf1  = sig1.getVal();
    double sigf2  = sig2.getVal();
    double sigwf1  = sigf1 /(sigf1 +sigf2 );
    double sigwf2  = sigf2 /(sigf1 +sigf2 );
    double c1 = a.getVal();
    double c2 = b.getVal();

    double sigmaf  = sqrt(sigmaf1 **2*sigwf1  + sigmaf2 **2*sigwf2 );
    double massmin  = meanf  - 5*sigmaf ;
    double massmax  = meanf  + 5*sigmaf ;

    int nmin  =  inputHist  ->GetXaxis()->FindBin(massmin );
    int nmax  =  inputHist  ->GetXaxis()->FindBin(massmax );
    int anmin  =  inputHist  ->GetXaxis()->FindBin(0.44);
    int anmax  =  inputHist  ->GetXaxis()->FindBin(0.56);

    double awyh1  =  inputHist  ->Integral(anmin ,nmin );
    double awyh2  =  inputHist  ->Integral(nmax ,anmax );
    double awyh  = awyh1  + awyh2 ;
    double totyh  =  inputHist  ->Integral(nmin ,nmax );

    x.setRange("cut",massmin ,massmax );
    RooAbsReal* ibkg  = poly.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig1  = gaus1.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig2  = gaus2.createIntegral(x,NormSet(x),Range("cut"));
    double ibkgf  = ibkg ->getVal();
    double bkgfe  = polysig.getError();
    double isig1f  = isig1 ->getVal();
    double isig2f  = isig2 ->getVal();

    double bkgy  = ibkgf *bkgf ;
    double bkgye  = ibkgf *bkgfe ;
    double sigy1  = isig1f *sigf1 ;
    double sigy2  = isig2f *sigf2 ;
    double sigy  = sigy1  + sigy2 ;
    double toty  = bkgy  + sigy ;

    double abkgy  = (1-ibkgf )*bkgf ;
    double asigy1  = (1-isig1f )*sigf1 ;
    double asigy2  = (1-isig2f )*sigf2 ;
    double asigy  = asigy1  + asigy2 ;
    double awy  = abkgy  + asigy ;

    double sigfrac  = sigy /toty ;
    double bkgfrac  = bkgy /toty ;

    double sigyh  = totyh  - bkgy ;
    double sigfrach  = sigyh /totyh ;
    double bkgfrach  = bkgy /totyh ;

    double signif  = sigyh / sqrt( totyh );

    return sigyh;
    
}

//Similarly for Lambda;
double la_YieldCal( TH1D* inputHist ){

//starting to do RooFit:

    RooRealVar x("x","mass",1,1.2);
    RooDataHist data("data","dataset",x, inputHist );
    RooPlot* xframe = x.frame(160);

    data.plotOn(xframe, Name("data"));

    RooRealVar mean("mean","mean",1.115,1.11,1.12);
    RooRealVar sigma1("sigma1","sigma1",0.003,0.001,0.01);
    RooRealVar sigma2("sigma2","sigma2",0.003,0.001,0.01);
    RooRealVar sig1("sig1","signal1",10,0,10000000);
    RooRealVar sig2("sig2","signal2",10,0,10000000);
    RooRealVar a("a","a",0,-100000,100000);
    RooRealVar b("b","b",0,-100000,100000);
    RooRealVar cp("cp","cp",0,-100000,100000);
    RooRealVar d("d","d",0,-100000,100000);

    RooRealVar f("f","f",0,-100000,100000);
    RooRealVar g("g","g",0,-100000,100000);
    RooRealVar h("h","h",0,-100000,100000);
    RooRealVar k("k","k",0,-100000,100000);

    RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
    RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
    RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
    RooRealVar polysig("polysig","polysig",10,0,10000000);
    RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

    x.setRange("cut",1.10,1.14);

    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));

    sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
    sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));

    xframe->Draw();

    double chi2  = xframe->chiSquare("sum","data");
    double meanf  = mean.getVal();
    double meanfe  = mean.getError();
    double sigmaf1  = sigma1.getVal();
    double sigmaf2  = sigma2.getVal();
    double bkgf  = polysig.getVal();
    double sigf1  = sig1.getVal();
    double sigf2  = sig2.getVal();
    double sigwf1  = sigf1 /(sigf1 +sigf2 );
    double sigwf2  = sigf2 /(sigf1 +sigf2 );
    double c1 = a.getVal();
    double c2 = b.getVal();

    double sigmaf  = sqrt(sigmaf1 **2*sigwf1  + sigmaf2 **2*sigwf2 );
    double massmin  = meanf  - 5*sigmaf ;
    double massmax  = meanf  + 5*sigmaf ;

    int nmin  = inputHist->GetXaxis()->FindBin(massmin );
    int nmax  = inputHist->GetXaxis()->FindBin(massmax );
    int anmin  = inputHist->GetXaxis()->FindBin(1.0);
    int anmax  = inputHist->GetXaxis()->FindBin(1.2);

    double awyh1  = inputHist->Integral(anmin ,nmin );
    double awyh2  = inputHist->Integral(nmax ,anmax );
    double awyh  = awyh1  + awyh2 ;
    double totyh  = inputHist->Integral(nmin ,nmax );

    x.setRange("cut",massmin ,massmax );
    RooAbsReal* ibkg  = poly.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig1  = gaus1.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig2  = gaus2.createIntegral(x,NormSet(x),Range("cut"));
    double ibkgf  = ibkg ->getVal();
    double bkgfe  = polysig.getError();
    double isig1f  = isig1 ->getVal();
    double isig2f  = isig2 ->getVal();

    double bkgy  = ibkgf *bkgf ;
    double bkgye  = ibkgf *bkgfe ;
    double sigy1  = isig1f *sigf1 ;
    double sigy2  = isig2f *sigf2 ;
    double sigy  = sigy1  + sigy2 ;
    double toty  = bkgy  + sigy ;

    double abkgy  = (1-ibkgf )*bkgf ;
    double asigy1  = (1-isig1f )*sigf1 ;
    double asigy2  = (1-isig2f )*sigf2 ;
    double asigy  = asigy1  + asigy2 ;
    double awy  = abkgy  + asigy ;

    double sigfrac  = sigy /toty ;
    double bkgfrac  = bkgy /toty ;

    double sigyh  = totyh  - bkgy ;
    double sigfrach  = sigyh /totyh ;
    double bkgfrach  = bkgy /totyh ;

    //double signif  = sigyh / sqrt( totyh );

    return sigyh;

}


void MBhistoSpectraRatio(){

    gStyle->SetErrorX(0);
/*
Getting the 3D histograms, and store in a 1D 3dimentional histogram:
 */

    TH3D* ksHist[4];
    TH3D* laHist[4];

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/MBpPb_MB1_July23_2014.root");
    ksHist[0] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[0] = (TH3D*)file1->Get("ana/InvMass_la_underlying");

    TFile* file2 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/MBpPb_MB2_July23_2014.root");
    ksHist[1] = (TH3D*)file2->Get("ana/InvMass_ks_underlying");
    laHist[1] = (TH3D*)file2->Get("ana/InvMass_la_underlying");

    TFile* file3 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/MBpPb_MB3_July23_2014.root");
    ksHist[2] = (TH3D*)file3->Get("ana/InvMass_ks_underlying");
    laHist[2] = (TH3D*)file3->Get("ana/InvMass_la_underlying");

    TFile* file4 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/MBpPb_MB4_July23_2014.root");
    ksHist[3] = (TH3D*)file4->Get("ana/InvMass_ks_underlying");
    laHist[3] = (TH3D*)file4->Get("ana/InvMass_la_underlying");
 
   
   double ks_norm[4];
   double la_norm[4];

        ks_norm[0] = ksHist[0]->GetEntries();
        ks_norm[1] = ksHist[1]->GetEntries();
        ks_norm[2] = ksHist[2]->GetEntries();
        ks_norm[3] = ksHist[3]->GetEntries();

        la_norm[0] = laHist[0]->GetEntries();
        la_norm[1] = laHist[1]->GetEntries();
        la_norm[2] = laHist[2]->GetEntries();
        la_norm[3] = laHist[3]->GetEntries();


    TH1D* ks_MB[4][6][20];
    TH1D* la_MB[4][6][20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};

    stringstream ksHistName;
    stringstream laHistName;

    for (int mult = 0; mult < 4; mult++){

        for (int eta = 0; eta < 6; eta++){

            for (int pt = 0; pt < 20; pt++){

                ksHistName.str("");
                laHistName.str("");

                ksHistName << "ks1_";
                ksHistName << mult;
                ksHistName << "_";
                ksHistName << eta;
                ksHistName << "_";
                ksHistName << pt;

                laHistName << "la1_";
                laHistName << mult;
                laHistName << "_";
                laHistName << eta;
                laHistName << "_";
                laHistName << pt;

                ks_MB[mult][eta][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );
                la_MB[mult][eta][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );
            }
        }
    }


/*
****************************************
*/


/**
 * Getting efficiency from the table:
 */

    TFile* t1 = new TFile("~/Desktop/Efficiency2D_V0_all_smooth10.root");
    
    TH2D* hnew1 = t1->Get("ks2Dnew");
    TH2D* hnew2 = t1->Get("la2Dnew");

    double ks_eff[6][20];
    double la_eff[6][20];

    for (int i = 0;i < 6;i++){

        for (int r = 0; r < 20; r++){

            ks_eff[i][r] = hnew1->GetBinContent(i+1,r+2);
            la_eff[i][r] = hnew2->GetBinContent(i+1,r+2);

        }
    }


/*
****************************************
 */

/*
Start to fit all histograms to obtain the eff_corr yields:
 */

    double ks_MB_yield[4][6][20];
    double la_MB_yield[4][6][20];

    for (mult = 0; mult < 4; mult++){
        
        for (eta = 0; eta < 6; eta++){

            for (pt = 0; pt < 20; pt++){

                ks_MB_yield[mult][eta][pt] = ks_YieldCal( ks_MB[mult][eta][pt] );
                    ks_MB_yield[mult][eta][pt] = ks_MB_yield[mult][eta][pt]/ ks_eff[eta][pt];
                la_MB_yield[mult][eta][pt] = la_YieldCal( la_MB[mult][eta][pt] ); 
                    la_MB_yield[mult][eta][pt] = la_MB_yield[mult][eta][pt]/la_eff[eta][pt];


            }
        }
    }

/*
*****************************************
 */

    double ks_MB_pTyield[4][20];
    double la_MB_pTyield[4][20];

    for (mult = 0; mult < 4; mult++){

        for (pt = 0; pt < 20; pt++){

            ks_MB_pTyield[mult][pt] = 0;
            la_MB_pTyield[mult][pt] = 0;

            for (eta = 0; eta < 6; eta++){

                ks_MB_pTyield[mult][pt] = ks_MB_yield[mult][eta][pt] + ks_MB_pTyield[mult][pt];
                la_MB_pTyield[mult][pt] = la_MB_yield[mult][eta][pt] + la_MB_pTyield[mult][pt];
            }

            cout << "ks_MB_pTyield: " << ks_MB_pTyield[mult][pt] << endl;
            cout << "la_MB_pTyield: " << la_MB_pTyield[mult][pt] << endl;

        }
    }

    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};

    TH1D* ksSpectra[4];
    TH1D* laSpectra[4];
    TH1D* ratioHist[4];

    stringstream ratioHistName;

    for (mult = 0; mult < 4; mult++){
       
        ksHistName.str("");
        laHistName.str("");
        ratioHistName.str("");

            ksHistName << "ksSpectra_";
            ksHistName << mult+1;

            laHistName << "laSpectra_";
            laHistName << mult+1;

            ratioHistName << "ratioHist_";
            ratioHistName << mult+1;

        ksSpectra[mult] = new TH1D(ksHistName.str().c_str(),ksHistName.str().c_str(),20,ptbins);
        laSpectra[mult] = new TH1D(laHistName.str().c_str(),laHistName.str().c_str(),20,ptbins);
        ratioHist[mult] = new TH1D(ratioHistName.str().c_str(),ratioHistName.str().c_str(),20,ptbins);
        
        for (pt = 0; pt < 20; pt++){

            double ks_temp = (ks_MB_pTyield[mult][pt]/binwidth[pt])/(2*3.1415926*ptbins[pt+1]*4.8*ks_norm[mult]);
            double la_temp = (la_MB_pTyield[mult][pt]/binwidth[pt])/(2*3.1415926*ptbins[pt+1]*4.8*la_norm[mult]);

            ksSpectra[mult]->SetBinContent(pt+1, ks_temp );
            ksSpectra[mult]->SetBinError(pt+1, sqrt((ks_MB_pTyield[mult][pt]/binwidth[pt]))/(2*3.1415926*ptbins[pt+1]*4.8*ks_norm[mult]));

            laSpectra[mult]->SetBinContent(pt+1, la_temp );
            laSpectra[mult]->SetBinError(pt+1, sqrt((la_MB_pTyield[mult][pt]/binwidth[pt]))/(2*3.1415926*ptbins[pt+1]*4.8*la_norm[mult]));

            ratioHist[mult]->SetBinContent(pt+1, la_MB_pTyield[mult][pt]/(2*ks_MB_pTyield[mult][pt]));
            double err = errorCal( (la_MB_pTyield[mult][pt]/binwidth[pt]), (ks_MB_pTyield[mult][pt]/binwidth[pt]) );
            ratioHist[mult]->SetBinError(pt+1, err );
        }

        cout << "pt: " << pt << endl;

        cout << "last pT bins KS: " << ks_MB_pTyield[mult][pt-1] << endl;

        cout << "last pT bins LA: " << la_MB_pTyield[mult][pt-1] << endl;

        cout << "last pT bins ratio: " << la_MB_pTyield[mult][pt-1]/(2*ks_MB_pTyield[mult][pt-1]) << endl;
    }

            TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
            w1->SetLineColor(kWhite);
            w1->AddEntry(ksSpectra[0],"0 < N^{offline}_{trk} < 35");
            w1->AddEntry(ksSpectra[1],"35 < N^{offline}_{trk} < 60");
            w1->AddEntry(ksSpectra[2],"60 < N^{offline}_{trk} < 90");
            w1->AddEntry(ksSpectra[3],"90 < N^{offline}_{trk} < 120");

            TCanvas* c1 = new TCanvas();

            TFile f2("MBbins.root","new");
            for (int you = 0; you < 4; you++){

                ksSpectra[you]->Write();
                laSpectra[you]->Write();
                ratioHist[you]->Write();
            }

            ksSpectra[0]->SetMarkerColor(kRed);
            ksSpectra[0]->SetMarkerStyle(22);
            ksSpectra[0]->SetLineColor(kRed);
            ksSpectra[0]->Scale(1);

            ksSpectra[1]->SetMarkerColor(kYellow);
            ksSpectra[1]->SetLineColor(kYellow);
            ksSpectra[1]->SetMarkerStyle(22);
            ksSpectra[1]->Scale(2);

            ksSpectra[2]->SetMarkerColor(kGreen);
            ksSpectra[2]->SetLineColor(kGreen);
            ksSpectra[2]->SetMarkerStyle(22);
            ksSpectra[2]->Scale(4);

            ksSpectra[3]->SetMarkerColor(kBlue);
            ksSpectra[3]->SetLineColor(kBlue);
            ksSpectra[3]->SetMarkerStyle(22);
            ksSpectra[3]->SetXTitle("P^{}_{T,V0}(GeV/ck)");
            ksSpectra[3]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}d#eta) [(GeV/c)^{-2}]");
            ksSpectra[3]->SetStats(kFALSE);
            ksSpectra[3]->Scale(8);

            ksSpectra[3]->Draw("P");
            ksSpectra[2]->Draw("Psame");
            ksSpectra[1]->Draw("Psame");
            ksSpectra[0]->Draw("Psame");

            w1->Draw("same");

            ksSpectra[3]->GetYaxis()->SetRangeUser(0.0000001,150);
            ksSpectra[3]->GetXaxis()->SetRangeUser(0,9);

            TCanvas* c2 = new TCanvas();

            laSpectra[0]->SetMarkerColor(kRed);
            laSpectra[0]->SetLineColor(kRed);
            laSpectra[0]->SetMarkerStyle(22);
            laSpectra[0]->Scale(1);

            laSpectra[1]->SetMarkerColor(kYellow);
            laSpectra[1]->SetLineColor(kYellow);
            laSpectra[1]->SetMarkerStyle(22);
            laSpectra[1]->Scale(2);

            laSpectra[2]->SetMarkerColor(kGreen);
            laSpectra[2]->SetLineColor(kGreen);
            laSpectra[2]->SetMarkerStyle(22);
            laSpectra[2]->Scale(4);

            laSpectra[3]->SetMarkerColor(kBlue);
            laSpectra[3]->SetLineColor(kBlue);
            laSpectra[3]->SetMarkerStyle(22);
            laSpectra[3]->SetStats(kFALSE);
            laSpectra[3]->SetXTitle("P^{}_{T,V0}(GeV/c)");
            laSpectra[3]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}d#eta) [(GeV/c)^{-2}]");
            laSpectra[3]->Scale(8);

            laSpectra[3]->Draw("P");
            laSpectra[2]->Draw("Psame");
            laSpectra[1]->Draw("Psame");
            laSpectra[0]->Draw("Psame");

            w1->Draw("same");

            laSpectra[3]->GetYaxis()->SetRangeUser(0.0000001,150);
            laSpectra[3]->GetXaxis()->SetRangeUser(0,9);


            TCanvas* c3 = new TCanvas();

            ratioHist[0]->SetMarkerColor(kRed);
            ratioHist[0]->SetMarkerStyle(22);
            ratioHist[0]->SetAxisRange(0,1.5,"Y");
            ratioHist[0]->SetYTitle("#Lambda+#bar{#Lambda}/2K^{0}_{s}");
            ratioHist[0]->SetXTitle("P^{}_{T,V0}(GeV/c)");
            ratioHist[0]->SetStats(kFALSE);
            ratioHist[0]->Draw("P");

            ratioHist[1]->SetMarkerColor(kYellow);
            ratioHist[1]->SetMarkerStyle(22);
            ratioHist[1]->Draw("Psame");

            ratioHist[2]->SetMarkerColor(kGreen);
            ratioHist[2]->SetMarkerStyle(22);
            ratioHist[2]->Draw("Psame");

            ratioHist[3]->SetMarkerColor(kBlue);
            ratioHist[3]->SetMarkerStyle(22);
            ratioHist[3]->Draw("Psame");


            w1->Draw("same");
        


}