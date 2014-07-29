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
#include "RooFitResult.h"


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
//
    //RooFitResult* rForward = forward.fitTo( *forwardBS, RooFit::PrintLevel(-1) );
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

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
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    //RooFitResult* rForward = forward.fitTo( *forwardBS, RooFit::PrintLevel(-1) );

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