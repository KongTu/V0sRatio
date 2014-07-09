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

    TFile* file = new TFile("~/Desktop/HMTripPb_histo_jetAnalysis_July7_2014.root");

    TH1D* ks_mass[10];
    TH1D* la_mass[10];
    TH1D* ks_underlying_mass[10];
    TH1D* la_underlying_mass[10];



    for (int value = 0; value < 10; value++){

        stringstream strr;
        strr << "ana/InvMass_ks";
        strr << value;

        stringstream strr1;
        strr1 << "ana/InvMass_la";
        strr1 << value;

        stringstream strr2;
        strr2 << "ana/InvMass_ks_underlying";
        strr2 << value;

        stringstream strr3;
        strr3 << "ana/InvMass_la_underlying";
        strr3 << value;

        file->GetObject(strr.str().c_str(),ks_mass[value] );
        file->GetObject(strr1.str().c_str(),la_mass[value] );
        file->GetObject(strr2.str().c_str(), ks_underlying_mass[value] );
        file->GetObject(strr3.str().c_str(), la_underlying_mass[value] );
    
    }

    float ks_yield[10];
    float la_yield[10];
    float ks_underlying_yield[10];
    float la_underlying_yield[10];

    TFile* t1 = new TFile("~/Desktop/Efficiency2D_V0_pPbHijing_counting.root");
    
    TH1D* hnew = Eff2D_ks->ProjectionY("ks",1,1000);
    TH1D* hnew2 = Eff2D_la->ProjectionY("la",1,1000);

//K0s in jet:

    TCanvas* c1 = new TCanvas();
    TH1F* h1 = new TH1F();
    
    c1->Print("ksHist_July8.pdf[");

    for (int it = 0; it < 10; it++){

        ks_yield[it] = ks_YieldCal( ks_mass[it] );
        h1->Fill( ks_yield[it] );
        c1->Print("ksHist_July8.pdf");
    }
    c1->Print("ksHist_July8.pdf]");

//K0s underlying:

    TCanvas* c1_underlying = new TCanvas();
    TH1F* h1_underlying = new TH1F();
    
    c1_underlying->Print("ksHist_July8_underlying.pdf[");

    for (int ip = 0; ip < 10; ip++){

        ks_underlying_yield[ip] = ks_YieldCal( ks_underlying_mass[ip] );
        h1_underlying->Fill( ks_underlying_yield[ip] );
        c1_underlying->Print("ksHist_July8_underlying.pdf");
    }
    c1_underlying->Print("ksHist_July8_underlying.pdf]");

//Lambda in jet:

    TCanvas* c2 = new TCanvas();
    TH1F* h2 = new TH1F();

    c2->Print("laHist_July8.pdf[");

    for (int is = 0; is < 10; is++){

        la_yield[is] = la_YieldCal ( la_mass[is] );
        h2->Fill( la_yield[is] );
        c2->Print("laHist_July8.pdf");

    }

    c2->Print("laHist_July8.pdf]");

//Lambda underlying:

    TCanvas* c2_underlying = new TCanvas();
    TH1F* h2_underlying = new TH1F();

    c2_underlying->Print("laHist_July8_underlying.pdf[");

    for (int iu = 0; iu < 10; iu++){

        la_underlying_yield[iu] = la_YieldCal ( la_underlying_mass[iu] );
        h2_underlying->Fill( la_underlying_yield[iu] );
        c2_underlying->Print("laHist_July8_underlying.pdf");

    }

    c2_underlying->Print("laHist_July8_underlying.pdf]");

    ks_yield[0] = ks_yield[0]/(0.0167*(hnew->GetBinContent(105)));
    ks_yield[1] = ks_yield[1]/(0.0167*(hnew->GetBinContent(180)));
    ks_yield[2] = ks_yield[2]/(0.0167*(hnew->GetBinContent(250)));
    ks_yield[3] = ks_yield[3]/(0.0167*(hnew->GetBinContent(330)));
    ks_yield[4] = ks_yield[4]/(0.0167*(hnew->GetBinContent(440)));
    ks_yield[5] = ks_yield[5]/(0.0167*(hnew->GetBinContent(570)));
    ks_yield[6] = ks_yield[6]/(0.0167*(hnew->GetBinContent(570)));
    ks_yield[7] = ks_yield[7]/(0.0167*(hnew->GetBinContent(570)));
    ks_yield[8] = ks_yield[8]/(0.0167*(hnew->GetBinContent(570)));
    ks_yield[9] = ks_yield[9]/(0.0167*(hnew->GetBinContent(570)));
    
    la_yield[0] = la_yield[0]/(0.0167*(hnew2->GetBinContent(105)));
    la_yield[1] = la_yield[1]/(0.0167*(hnew2->GetBinContent(180)));
    la_yield[2] = la_yield[2]/(0.0167*(hnew2->GetBinContent(250)));
    la_yield[3] = la_yield[3]/(0.0167*(hnew2->GetBinContent(330)));
    la_yield[4] = la_yield[4]/(0.0167*(hnew2->GetBinContent(440)));
    la_yield[5] = la_yield[5]/(0.0167*(hnew2->GetBinContent(570)));
    la_yield[6] = la_yield[6]/(0.0167*(hnew2->GetBinContent(570)));
    la_yield[7] = la_yield[7]/(0.0167*(hnew2->GetBinContent(570)));
    la_yield[8] = la_yield[8]/(0.0167*(hnew2->GetBinContent(570)));
    la_yield[9] = la_yield[9]/(0.0167*(hnew2->GetBinContent(570)));


    ks_underlying_yield[0] = ks_underlying_yield[0]/(0.0167*(hnew->GetBinContent(105)));
    ks_underlying_yield[1] = ks_underlying_yield[1]/(0.0167*(hnew->GetBinContent(180)));
    ks_underlying_yield[2] = ks_underlying_yield[2]/(0.0167*(hnew->GetBinContent(250)));
    ks_underlying_yield[3] = ks_underlying_yield[3]/(0.0167*(hnew->GetBinContent(330)));
    ks_underlying_yield[4] = ks_underlying_yield[4]/(0.0167*(hnew->GetBinContent(440)));
    ks_underlying_yield[5] = ks_underlying_yield[5]/(0.0167*(hnew->GetBinContent(570)));
    ks_underlying_yield[6] = ks_underlying_yield[6]/(0.0167*(hnew->GetBinContent(570)));
    ks_underlying_yield[7] = ks_underlying_yield[7]/(0.0167*(hnew->GetBinContent(570)));
    ks_underlying_yield[8] = ks_underlying_yield[8]/(0.0167*(hnew->GetBinContent(570)));
    ks_underlying_yield[9] = ks_underlying_yield[9]/(0.0167*(hnew->GetBinContent(570)));
    
    la_underlying_yield[0] = la_underlying_yield[0]/(0.0167*(hnew2->GetBinContent(105)));
    la_underlying_yield[1] = la_underlying_yield[1]/(0.0167*(hnew2->GetBinContent(180)));
    la_underlying_yield[2] = la_underlying_yield[2]/(0.0167*(hnew2->GetBinContent(250)));
    la_underlying_yield[3] = la_underlying_yield[3]/(0.0167*(hnew2->GetBinContent(330)));
    la_underlying_yield[4] = la_underlying_yield[4]/(0.0167*(hnew2->GetBinContent(440)));
    la_underlying_yield[5] = la_underlying_yield[5]/(0.0167*(hnew2->GetBinContent(570)));
    la_underlying_yield[6] = la_underlying_yield[6]/(0.0167*(hnew2->GetBinContent(570)));
    la_underlying_yield[7] = la_underlying_yield[7]/(0.0167*(hnew2->GetBinContent(570)));
    la_underlying_yield[8] = la_underlying_yield[8]/(0.0167*(hnew2->GetBinContent(570)));
    la_underlying_yield[9] = la_underlying_yield[9]/(0.0167*(hnew2->GetBinContent(570)));

//Lambda/K0s ratio

    for (int i = 0; i < 10; i++){

        cout << "let's see K0short yield " << i << ": " << ks_yield[i] << endl;
        cout << "let's ses K0short underlying yield" << i << ": " << ks_underlying_yield[i] << endl;
        cout << "let's see Lambda yield " << i << ": " << la_yield[i] << endl;
        cout << "let's see Lambda underlying yield" << i << ": " << la_underlying_yield[i] << endl;
    
    }

    double ptbins[] = {0.7,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0};

    TH1F* h3 = new TH1F("h3","h3",10,ptbins);
    TH1F* h4 = new TH1F("h4","h4",10,ptbins);

    for (int y = 0; y < 10; y++){

        float temp = la_yield[y]/( 2 * ks_yield[y] );
        float temp1 = la_underlying_yield[y]/( 2 * ks_underlying_yield[y] );
        h3->SetBinContent(y+1,temp );
        h4->SetBinContent(y+1,temp1 );

    }

    TCanvas* c3 = new TCanvas();
    h3->SetMarkerStyle(21);
    h3->SetMarkerColor(kRed);
    h3->SetAxisRange(0,1,"Y");
    h3->SetXTitle("P^{}_{T,V0}(GeV/c)");
    h3->SetYTitle("#Events");
    h3->SetStats(kFALSE);
    
    h4->SetMarkerStyle(21);
    h4->SetMarkerColor(kBlue);
    h4->SetAxisRange(0,1,"Y");

    TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
    w1->AddEntry(h3,"Inside Jet V0s");
    w1->AddEntry(h4,"All V0s");
    
    h3->Draw("P");
    h4->Draw("Psame");
    w1->Draw("same");


}