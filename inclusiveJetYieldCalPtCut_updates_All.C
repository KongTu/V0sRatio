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
    double massmin  = meanf  - 2*sigmaf ;
    double massmax  = meanf  + 2*sigmaf ;

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
    double massmin  = meanf  - 2*sigmaf ;
    double massmax  = meanf  + 2*sigmaf ;

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

void inclusiveJetYieldCalPtCut_updates_All(){

/**
 * MB sample, with Ntrk > 50; 4 bins:
 */

    TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_histo/MBpPb_histo_20pTBins_total.root");

    TH1D* ks_0_underlying_mass[20];
    TH1D* la_0_underlying_mass[20];

    TH1D* ks_1_underlying_mass[20];
    TH1D* la_1_underlying_mass[20];

    TH1D* ks_2_underlying_mass[20];
    TH1D* la_2_underlying_mass[20];

    TH1D* ks_3_underlying_mass[20];
    TH1D* la_3_underlying_mass[20];

    for (int value = 0; value < 20; value++){

        stringstream strr0;
        strr0 << "ana/InvMass_ks_0_underlying";
        strr0 << value;

        stringstream strr_la;
        strr_la << "ana/InvMass_la_0_underlying";
        strr_la << value;

        stringstream strr;
        strr << "ana/InvMass_ks_1_underlying";
        strr << value;

        stringstream strr1;
        strr1 << "ana/InvMass_la_1_underlying";
        strr1 << value;

        stringstream strr2;
        strr2 << "ana/InvMass_ks_2_underlying";
        strr2 << value;

        stringstream strr3;
        strr3 << "ana/InvMass_la_2_underlying";
        strr3 << value;

        stringstream strr4;
        strr4 << "ana/InvMass_ks_3_underlying";
        strr4 << value;

        stringstream strr5;
        strr5 << "ana/InvMass_la_3_underlying";
        strr5 << value;

        file->GetObject(strr0.str().c_str(), ks_0_underlying_mass[value] );
        file->GetObject(strr_la.str().c_str(), la_0_underlying_mass[value] );
        file->GetObject(strr.str().c_str(), ks_1_underlying_mass[value] );
        file->GetObject(strr1.str().c_str(), la_1_underlying_mass[value] );
        file->GetObject(strr2.str().c_str(), ks_2_underlying_mass[value] );
        file->GetObject(strr3.str().c_str(), la_2_underlying_mass[value] );
        file->GetObject(strr4.str().c_str(), ks_3_underlying_mass[value] );
        file->GetObject(strr5.str().c_str(), la_3_underlying_mass[value] );
        
    }


/**
 * HM trigger data: Ntrk > 120; 4 bins
 */

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_histo/HMpPb_histo_20pTBins_total.root");

    TH1D* ks_4_underlying_mass[20];
    TH1D* la_4_underlying_mass[20];

    TH1D* ks_5_underlying_mass[20];
    TH1D* la_5_underlying_mass[20];

    TH1D* ks_6_underlying_mass[20];
    TH1D* la_6_underlying_mass[20];

    TH1D* ks_7_underlying_mass[20];
    TH1D* la_7_underlying_mass[20];

    TH1D* ks_8_underlying_mass[20];
    TH1D* la_8_underlying_mass[20];

    for (int value = 0; value < 20; value++){

        stringstream HMstrr0;
        HMstrr0 << "ana/InvMass_ks_0_underlying";
        HMstrr0 << value;

        stringstream HMstrr_la;
        HMstrr_la << "ana/InvMass_la_0_underlying";
        HMstrr_la << value;

        stringstream HMstrr;
        HMstrr << "ana/InvMass_ks_1_underlying";
        HMstrr << value;

        stringstream HMstrr1;
        HMstrr1 << "ana/InvMass_la_1_underlying";
        HMstrr1 << value;

        stringstream HMstrr2;
        HMstrr2 << "ana/InvMass_ks_2_underlying";
        HMstrr2 << value;

        stringstream HMstrr3;
        HMstrr3 << "ana/InvMass_la_2_underlying";
        HMstrr3 << value;

        stringstream HMstrr4;
        HMstrr4 << "ana/InvMass_ks_3_underlying";
        HMstrr4 << value;

        stringstream HMstrr5;
        HMstrr5 << "ana/InvMass_la_3_underlying";
        HMstrr5 << value;

        stringstream HMstrr6;
        HMstrr6 << "ana/InvMass_ks_260_underlying";
        HMstrr6 << value;

        stringstream HMstrr7;
        HMstrr7 << "ana/InvMass_la_260_underlying";
        HMstrr7 << value;

        file1->GetObject(HMstrr0.str().c_str(), ks_4_underlying_mass[value] );
        file1->GetObject(HMstrr_la.str().c_str(), la_4_underlying_mass[value] );
        file1->GetObject(HMstrr.str().c_str(), ks_5_underlying_mass[value] );
        file1->GetObject(HMstrr1.str().c_str(), la_5_underlying_mass[value] );
        file1->GetObject(HMstrr2.str().c_str(), ks_6_underlying_mass[value] );
        file1->GetObject(HMstrr3.str().c_str(), la_6_underlying_mass[value] );
        file1->GetObject(HMstrr4.str().c_str(), ks_7_underlying_mass[value] );
        file1->GetObject(HMstrr5.str().c_str(), la_7_underlying_mass[value] );
        file1->GetObject(HMstrr6.str().c_str(), ks_8_underlying_mass[value] );
        file1->GetObject(HMstrr7.str().c_str(), la_8_underlying_mass[value] );

    }

    TH1D* ks_underlying_mass[9][20];
    TH1D* la_underlying_mass[9][20];

   
    for (int y = 0; y < 20; y++){

        ks_underlying_mass[0][y] = ks_0_underlying_mass[y];
        ks_underlying_mass[1][y] = ks_1_underlying_mass[y];
        ks_underlying_mass[2][y] = ks_2_underlying_mass[y];
        ks_underlying_mass[3][y] = ks_3_underlying_mass[y];
        ks_underlying_mass[4][y] = ks_4_underlying_mass[y];
        ks_underlying_mass[5][y] = ks_5_underlying_mass[y];
        ks_underlying_mass[6][y] = ks_6_underlying_mass[y];
        ks_underlying_mass[7][y] = ks_7_underlying_mass[y];
        ks_underlying_mass[8][y] = ks_8_underlying_mass[y];

        la_underlying_mass[0][y] = la_0_underlying_mass[y];
        la_underlying_mass[1][y] = la_1_underlying_mass[y];
        la_underlying_mass[2][y] = la_2_underlying_mass[y];
        la_underlying_mass[3][y] = la_3_underlying_mass[y];
        la_underlying_mass[4][y] = la_4_underlying_mass[y];
        la_underlying_mass[5][y] = la_5_underlying_mass[y];
        la_underlying_mass[6][y] = la_6_underlying_mass[y];
        la_underlying_mass[7][y] = la_7_underlying_mass[y];
        la_underlying_mass[8][y] = la_8_underlying_mass[y];

    }
  

    float ks_underlying_yield[9][20];
    float la_underlying_yield[9][20];

    TFile* t1 = new TFile("~/Desktop/Efficiency2D_V0_pPbHijing_counting.root");
    
    TH1D* hnew = Eff2D_ks->ProjectionY("ks",1,1000);
    TH1D* hnew2 = Eff2D_la->ProjectionY("la",1,1000);

    float ks_eff[20];
    float la_eff[20];

    ks_eff[0] = (0.0167*(hnew->GetBinContent(60)));
    ks_eff[1] = (0.0167*(hnew->GetBinContent(80)));
    ks_eff[2] = (0.0167*(hnew->GetBinContent(100)));
    ks_eff[3] = (0.0167*(hnew->GetBinContent(120)));
    ks_eff[4] = (0.0167*(hnew->GetBinContent(140)));
    ks_eff[5] = (0.0167*(hnew->GetBinContent(160)));
    ks_eff[6] = (0.0167*(hnew->GetBinContent(180)));
    ks_eff[7] = (0.0167*(hnew->GetBinContent(200)));
    ks_eff[8] = (0.0167*(hnew->GetBinContent(220)));
    ks_eff[9] = (0.0167*(hnew->GetBinContent(240)));
    ks_eff[10] = (0.0167*(hnew->GetBinContent(260)));
    ks_eff[11] = (0.0167*(hnew->GetBinContent(280)));
    ks_eff[12] = (0.0167*(hnew->GetBinContent(300)));
    ks_eff[13] = (0.0167*(hnew->GetBinContent(340)));
    ks_eff[14] = (0.0167*(hnew->GetBinContent(380)));
    ks_eff[15] = (0.0167*(hnew->GetBinContent(420)));
    ks_eff[16] = (0.0167*(hnew->GetBinContent(460)));
    ks_eff[17] = (0.0167*(hnew->GetBinContent(500)));
    ks_eff[18] = (0.0167*(hnew->GetBinContent(560)));
    ks_eff[19] = (0.0167*(hnew->GetBinContent(570)));

    la_eff[0] = (0.0167*(hnew2->GetBinContent(60)));
    la_eff[1] = (0.0167*(hnew2->GetBinContent(80)));
    la_eff[2] = (0.0167*(hnew2->GetBinContent(100)));
    la_eff[3] = (0.0167*(hnew2->GetBinContent(120)));
    la_eff[4] = (0.0167*(hnew2->GetBinContent(140)));
    la_eff[5] = (0.0167*(hnew2->GetBinContent(160)));
    la_eff[6] = (0.0167*(hnew2->GetBinContent(180)));
    la_eff[7] = (0.0167*(hnew2->GetBinContent(200)));
    la_eff[8] = (0.0167*(hnew2->GetBinContent(220)));
    la_eff[9] = (0.0167*(hnew2->GetBinContent(240)));
    la_eff[10] = (0.0167*(hnew2->GetBinContent(260)));
    la_eff[11] = (0.0167*(hnew2->GetBinContent(280)));
    la_eff[12] = (0.0167*(hnew2->GetBinContent(300)));
    la_eff[13] = (0.0167*(hnew2->GetBinContent(340)));
    la_eff[14] = (0.0167*(hnew2->GetBinContent(380)));
    la_eff[15] = (0.0167*(hnew2->GetBinContent(420)));
    la_eff[16] = (0.0167*(hnew2->GetBinContent(460)));
    la_eff[17] = (0.0167*(hnew2->GetBinContent(500)));
    la_eff[18] = (0.0167*(hnew2->GetBinContent(560)));
    la_eff[19] = (0.0167*(hnew2->GetBinContent(570)));


    TCanvas* c1[9];
    TCanvas* c2[9];

    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};

    TH1F* h1[9];
    TH1F* h2[9]; 
    TH1F* h3[9];

    int it = 0;
    stringstream ks_name;
    stringstream la_name;


/*
start to fit and print the histos:
 */

    for ( it = 0; it < 9; it++){

        c1[it] = new TCanvas();
        h1[it] = new TH1F(Form("h1%d",it),"h1",20,ptbins);

        ks_name.str("");
        ks_name << "ksHist_July18_underlying_";
        ks_name << it;
        ks_name << ".pdf[";

            c1[it]->Print( ks_name.str().c_str() );

        ks_name.str("");
        ks_name << "ksHist_July18_underlying_";
        ks_name << it;
        ks_name << ".pdf";

        for (int is = 0; is < 20; is++){

            ks_underlying_yield[it][is] = ks_YieldCal( ks_underlying_mass[it][is] );
            c1[it]->Print( ks_name.str().c_str() );
            ks_underlying_yield[it][is] = ks_underlying_yield[it][is]/ks_eff[is];
            h1[it]->SetBinContent(is+1, (ks_underlying_yield[it][is]/binwidth[is]) );
            h1[it]->SetBinError(is+1, sqrt( ks_underlying_yield[it][is] ) );

        }

        ks_name.str("");
        ks_name << "ksHist_July18_underlying_";
        ks_name << it;
        ks_name << ".pdf]";

        c1[it]->Print( ks_name.str().c_str() );

    }


    for ( it = 0; it < 9; it++){

        c2[it] = new TCanvas();
        h2[it] = new TH1F(Form("h2%d",it),"h2",20,ptbins);

        la_name.str("");
        la_name << "laHist_July18_underlying_";
        la_name << it;
        la_name << ".pdf[";

            c2[it]->Print( la_name.str().c_str() );

        la_name.str("");
        la_name << "laHist_July18_underlying_";
        la_name << it;
        la_name << ".pdf";

        for (int is = 0; is < 20; is++){

            la_underlying_yield[it][is] = la_YieldCal( la_underlying_mass[it][is] );
            c2[it]->Print( la_name.str().c_str() );
            la_underlying_yield[it][is] = la_underlying_yield[it][is]/la_eff[is];
            h2[it]->SetBinContent(is+1, (la_underlying_yield[it][is]/binwidth[is]) );
            h2[it]->SetBinError(is+1, sqrt( la_underlying_yield[it][is] ) );

        }

        la_name.str("");
        la_name << "laHist_July18_underlying_";
        la_name << it;
        la_name << ".pdf]";

        c2[it]->Print( la_name.str().c_str() );

    }

/*
Calculate the Lambda/K0short ratio as well as the errors:
 */

    for (it = 0; it < 9; it++){

        h3[it] =  new TH1F(Form("h3%d",it),"#Lambda/K^_{0}_{s} ratio in MB+HM pPb data (full stats)",20,ptbins);

        for (is = 0; is < 20; is++){

            double err = errorCal( la_underlying_yield[it][is], ks_underlying_yield[it][is] );
            double temp = la_underlying_yield[it][is]/( 2 * ks_underlying_yield[it][is] );
            h3[it]->SetBinContent(is+1, temp );
            h3[it]->SetBinError(is+1, err );
        }

    }

    TCanvas* c3 = new TCanvas();

    for(it = 0; it < 9; it++){

        h3[it]->SetMarkerStyle(22);
        h3[it]->SetAxisRange(0,1,"Y");
        h3[it]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        h3[it]->SetYTitle("#Lambda+#bar{#Lambda}/2K^{0}_{s}");
        h3[it]->SetStats(kFALSE);
    }

    h3[0]->SetMarkerColor(kRed);
    h3[0]->SetLineColor(kRed);

    h3[1]->SetMarkerColor(kOrange);
    h3[1]->SetLineColor(kOrange);

    h3[2]->SetMarkerColor(kYellow);
    h3[2]->SetLineColor(kYellow);

    h3[3]->SetMarkerColor(kGreen);
    h3[3]->SetLineColor(kGreen);

    h3[4]->SetMarkerColor(kCyan);
    h3[4]->SetLineColor(kCyan);

    h3[5]->SetMarkerColor(kBlue);
    h3[5]->SetLineColor(kBlue);

    h3[6]->SetMarkerColor(kMagenta);
    h3[6]->SetLineColor(kMagenta);

    h3[7]->SetMarkerColor(kBlack);
    h3[7]->SetLineColor(kBlack);

    h3[8]->SetMarkerColor(kPink);
    h3[8]->SetLineColor(kPink);


    TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
    w1->AddEntry(h3[0],"0 < N^{offline}_{trk} < 35");
    w1->AddEntry(h3[1],"35 < N^{offline}_{trk} < 60");
    w1->AddEntry(h3[2],"60 < N^{offline}_{trk} < 90");
    w1->AddEntry(h3[3],"90 < N^{offline}_{trk} < 120");
    w1->AddEntry(h3[4],"120 < N^{offline}_{trk} < 150");
    w1->AddEntry(h3[5],"150 < N^{offline}_{trk} < 185");
    w1->AddEntry(h3[6],"185 < N^{offline}_{trk} < 220");
    w1->AddEntry(h3[7],"220 < N^{offline}_{trk} < 260");
    w1->AddEntry(h3[8],"N^{offline}_{trk} > 260");
    
    h3[0]->Draw("P");
    for (int t = 1; t < 9; t++){

        h3[t]->Draw("Psame");
    }

    w1->Draw("same");


    for (int u = 0; u < 9; u++){

        h1[u]->SetStats(kFALSE);
        h1[u]->SetMarkerStyle(22);
        h1[u]->SetAxisRange(1000,1000000000,"Y");
        h1[u]->SetYTitle("#Events");
        h1[u]->SetXTitle("P^_{}_{T,V0} (GeV/c)");

        h2[u]->SetStats(kFALSE);
        h2[u]->SetMarkerStyle(22);
        h2[u]->SetAxisRange(1000,1000000000,"Y");
        h2[u]->SetYTitle("#Events");
        h2[u]->SetXTitle("P^_{}_{T,V0} (GeV/c)");

    }

    h1[0]->SetMarkerColor(kRed);
    h1[1]->SetMarkerColor(kOrange);
    h1[2]->SetMarkerColor(kYellow);
    h1[3]->SetMarkerColor(kGreen);
    h1[4]->SetMarkerColor(kCyan);
    h1[5]->SetMarkerColor(kBlue);
    h1[6]->SetMarkerColor(kMagenta);
    h1[7]->SetMarkerColor(kBlack);
    h1[8]->SetMarkerColor(kPink);
    h1[0]->SetLineColor(kRed);
    h1[1]->SetLineColor(kOrange);
    h1[2]->SetLineColor(kYellow);
    h1[3]->SetLineColor(kGreen);
    h1[4]->SetLineColor(kCyan);
    h1[5]->SetLineColor(kBlue);
    h1[6]->SetLineColor(kMagenta);
    h1[7]->SetLineColor(kBlack);
    h1[8]->SetLineColor(kPink);

    h2[0]->SetMarkerColor(kRed);
    h2[1]->SetMarkerColor(kOrange);
    h2[2]->SetMarkerColor(kYellow);
    h2[3]->SetMarkerColor(kGreen);
    h2[4]->SetMarkerColor(kCyan);
    h2[5]->SetMarkerColor(kBlue);
    h2[6]->SetMarkerColor(kMagenta);
    h2[7]->SetMarkerColor(kBlack);
    h2[8]->SetMarkerColor(kPink);
    h2[0]->SetLineColor(kRed);
    h2[1]->SetLineColor(kOrange);
    h2[2]->SetLineColor(kYellow);
    h2[3]->SetLineColor(kGreen);
    h2[4]->SetLineColor(kCyan);
    h2[5]->SetLineColor(kBlue);
    h2[6]->SetLineColor(kMagenta);
    h2[7]->SetLineColor(kBlack);
    h2[8]->SetLineColor(kPink);


    TCanvas* c4 = new TCanvas();

    TLegend *w2 = new TLegend(0.25,0.4,0.5,0.5);
    w2->AddEntry(h1[0],"0 < N^{offline}_{trk} < 35");
    w2->AddEntry(h1[1],"35 < N^{offline}_{trk} < 60");
    w2->AddEntry(h1[2],"60 < N^{offline}_{trk} < 90");
    w2->AddEntry(h1[3],"90 < N^{offline}_{trk} < 120");
    w2->AddEntry(h1[4],"120 < N^{offline}_{trk} < 150");
    w2->AddEntry(h1[5],"150 < N^{offline}_{trk} < 185");
    w2->AddEntry(h1[6],"185 < N^{offline}_{trk} < 220");
    w2->AddEntry(h1[7],"220 < N^{offline}_{trk} < 260");
    w2->AddEntry(h1[8],"N^{offline}_{trk} > 260");
    
    h1[0]->Draw("P");

    for (u = 1; u < 9; u++){

        h1[u]->Draw("Psame");
    }

    w2->Draw("same");

    TCanvas* c5 = new TCanvas();

    TLegend *w3 = new TLegend(0.25,0.4,0.5,0.5);
    w3->AddEntry(h2[0],"0 < N^{offline}_{trk} < 35");
    w3->AddEntry(h2[1],"35 < N^{offline}_{trk} < 60");
    w3->AddEntry(h2[2],"60 < N^{offline}_{trk} < 90");
    w3->AddEntry(h2[3],"90 < N^{offline}_{trk} < 120");
    w3->AddEntry(h2[4],"120 < N^{offline}_{trk} < 150");
    w3->AddEntry(h2[5],"150 < N^{offline}_{trk} < 185");
    w3->AddEntry(h2[6],"185 < N^{offline}_{trk} < 220");
    w3->AddEntry(h2[7],"220 < N^{offline}_{trk} < 260");
    w3->AddEntry(h2[8],"N^{offline}_{trk} > 260");
    
    h2[0]->Draw("P");

    for (u = 1; u < 9; u++){

        h2[u]->Draw("Psame");
    }

    w3->Draw("same");
    








}