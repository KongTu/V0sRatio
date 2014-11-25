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

vector<double> x,ex, y,ey,x1,ex1,y1,ey1;

double getChi2;
int getNdf;

const double R = 1.0;
const double dr = 1e-2; // FIXME

double integral(double beta_s, double T, double n, double pt, double mt)
{
  
  double s = 0.;
  for(double r = dr/2; r < R; r += dr)
  {
    double beta_r = beta_s * TMath::Power((r/R),n);
    double rho = TMath::ATanH(beta_r);

    s += r * dr * TMath::BesselK1((mt*cosh(rho))/T) * TMath::BesselI0((pt*sinh(rho))/T);
  }

  return s;
}

void function(int &npar, double *gin, double &f, double *par, int flag)
{
  double aka    = par[0];
  double aka1   = par[4];
  double T      = par[2];
  double beta_s = par[1];
  double n = par[3];

  double chi2 = 0.;

  for(int i = 0; i < int(x.size()); i++)
  {
    double m = 0.497;

    double pt = x[i];
    double mt = sqrt(m*m + pt*pt);

    double a = 0.;
    a = aka;
  
    double theo = a * mt * integral(beta_s, T, n, pt, mt);
    double q = (y[i] - theo)/ey[i];

    chi2 += q*q;
  }

  for(i = 0; i < int(x1.size()); i++){

  	double m1 = 1.115;

  	double pt1 = x1[i];
  	double mt1 = sqrt(m1*m1 + pt1*pt1);

  	double a1 = 0.;
  	a1 = aka1;

  	double theo1 = a1 * mt1 * integral(beta_s,T,n,pt1,mt1);
  	double q1 = (y1[i]-theo1)/ey1[i];

  	chi2 += q1*q1;
  }

  f = chi2;

  getChi2 = f;
  getNdf  = x.size() + x1.size() - 5 - 1;
}

double MyFunc( double *pt, double *p){

  double mass = 0.497;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}

double MyFunc1( double *pt, double *p){

  double mass = 1.115;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}

double beta_T(double beta_s,double n){

  return (beta_s*2)/(n+2);

}


void MinuitBlastWaveCMSdataFitRapidity(){

  gStyle->SetErrorX(0);


  TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/SpectraRapidityProducer/SpectraRapidityProducer220plus_vtxReweighEPOS_26ksPTbins_v4.root");
  
  TH1D* ksSpectra[5];
  TH1D* laSpectra[5];

  stringstream ksHistName;
  stringstream laHistName;

  for (int rpy = 0; rpy < 5; rpy++){

    ksHistName.str("");
    laHistName.str("");

    ksHistName << "ksSpectra_";
    ksHistName << rpy+1;

    laHistName << "laSpectra_";
    laHistName << rpy+1;

    ksSpectra[rpy] = (TH1D*)file->Get(ksHistName.str().c_str());
    laSpectra[rpy] = (TH1D*)file->Get(laHistName.str().c_str());
    

  }

  double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_ptbincenter[26] = {0.15,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
  double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

  TCanvas* c1 = new TCanvas();
    c1->Divide(2,3,0,0);

  double beta_s[5],Tkin[5],aka[5],aka1[5],n[5];
  double ebeta_s[5],eTkin[5],eaka[5],eaka1[5],en[5];

  stringstream f1Name;
  stringstream f2Name;

  TF1* f1[5];
  TF1* f2[5];

  TLatex* ratio[5];
  ratio[0] = new TLatex(0.57,0.17,"-2.87 < y < -1.8");
  ratio[1] = new TLatex(0.57,0.17,"-1.8 < y < -0.9");
  ratio[2] = new TLatex(0.57,0.17,"-0.9 < y < 0");
  ratio[3] = new TLatex(0.57,0.17,"0 < y < 0.93");
  ratio[4] = new TLatex(0.57,0.17,"0.93 < y < 1.93");

  TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
  w1->SetLineColor(kWhite);
  w1->SetFillColor(0);

  w1->AddEntry(ksSpectra[4],"K^{0}_{s}");
  w1->AddEntry(laSpectra[4],"#Lambda/#bar{#Lambda}");
  
  TGraph* test[8];

  for(int rpy = 0; rpy < 5; rpy++){

      x.clear();
      ex.clear();
      y.clear();
      ey.clear();

      x1.clear();
      ex1.clear();
      y1.clear();
      ey1.clear();

      for (int pt = 1; pt < 13 ; pt++){

        x.push_back( ks_ptbincenter[pt] );
        ex.push_back(0.0);
        y.push_back( ksSpectra[rpy]->GetBinContent(pt+1) );
        ey.push_back( ksSpectra[rpy]->GetBinError(pt+1) );
        
      }

      for(pt = 1; pt < 13; pt++){

        x1.push_back( la_ptbincenter[pt] );
        ex1.push_back(0.0);
        y1.push_back( laSpectra[rpy]->GetBinContent(pt+1) );
        ey1.push_back( laSpectra[rpy]->GetBinError(pt+1) );
        
      }

  /**
   * simultaneous fit for K0s and Lambda:
   */

      TMinuit * gMinuit[5];
      gMinuit[rpy] = new TMinuit(5);

      //set the function to minimize with minuit;
      gMinuit[rpy]->SetFCN(function);

      double arglist[10];
      arglist[4] = 0.001;
      gMinuit[rpy]->mnexcm("SET ERR", arglist, 1, 0);

      gMinuit[rpy]->mnparm(0, "aka",    10,   0.1,  1,    100000, 0);
      gMinuit[rpy]->mnparm(1, "beta_s", 0.70, 0.01, 0.15, 1.0,   0);
      gMinuit[rpy]->mnparm(2, "Tkin",   0.15, 0.01, 0.05, 1.0,   0);
      gMinuit[rpy]->mnparm(3, "n",      1.0,  0.01, 0.1,  10.0,  0);
      gMinuit[rpy]->mnparm(4, "aka1",   10,   0.1,  1,    100000, 0);

      gMinuit[rpy]->mnexcm("MIGRAD",    arglist,  1,   0);
      gMinuit[rpy]->mnexcm("CALL FCN",  arglist,  1,   0);
      //gMinuit[rpy]->mnexcm("HESSE",     arglist,  1,   0);

      gMinuit[rpy]->GetParameter(0, aka[rpy],    eaka[rpy]);
      gMinuit[rpy]->GetParameter(1, beta_s[rpy], ebeta_s[rpy]);
      gMinuit[rpy]->GetParameter(2, Tkin[rpy],   eTkin[rpy]);
      gMinuit[rpy]->GetParameter(3, n[rpy],      en[rpy]);
      gMinuit[rpy]->GetParameter(4, aka1[rpy],   eaka1[rpy]);

      gMinuit[rpy]->SetErrorDef(4);
      test[rpy] = (TGraph*)gMinuit[rpy]->Contour(100,1,2);

  /**
   * define the fit function by using the parameters from fit:
   */

      double xmin = 0.35;
      double xmax = 3.0;

      double xmin1 = 0.8;
      double xmax1 = 3.3;

      f1Name.str("");
      f2Name.str("");

      f1Name << "f1_";
      f1Name << rpy+1;

      f2Name << "f2_";
      f2Name << rpy+1;
  
      f1[rpy] = new TF1(f1Name.str().c_str(),MyFunc,xmin,xmax,4);
      f1[rpy]->SetParameter(0,beta_s[rpy]);
      f1[rpy]->SetParameter(1,Tkin[rpy]);
      f1[rpy]->SetParameter(2,aka[rpy]);
      f1[rpy]->SetParameter(3,n[rpy]);

      f2[rpy] = new TF1(f2Name.str().c_str(),MyFunc1,xmin1,xmax1,4);
      f2[rpy]->SetParameter(0,beta_s[rpy]);
      f2[rpy]->SetParameter(1,Tkin[rpy]);
      f2[rpy]->SetParameter(2,aka1[rpy]);
      f2[rpy]->SetParameter(3,n[rpy]);

      c1->cd(rpy+1);
      gPad->SetLogy(1);

      ksSpectra[rpy]->SetMarkerSize(1.3);
      ksSpectra[rpy]->SetMarkerStyle(20);
      ksSpectra[rpy]->SetMarkerColor(kBlue);
      ksSpectra[rpy]->SetStats(kFALSE);
      ksSpectra[rpy]->SetTitle("");
      ksSpectra[rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
      ksSpectra[rpy]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");

      ksSpectra[rpy]->GetYaxis()->SetRangeUser(0.0000001,1000);

      laSpectra[rpy]->SetMarkerSize(1.3);
      laSpectra[rpy]->SetMarkerStyle(20);
      laSpectra[rpy]->SetMarkerColor(kBlack);

      ksSpectra[rpy]->Draw("P");
      laSpectra[rpy]->Draw("Psame");
      ratio[rpy]->Draw("same");

      f2[rpy]->Draw("same");
      f1[rpy]->Draw("same");

  }

  w1->Draw("same");

  TCanvas* s1 = new TCanvas();

  TH1D* hist = new TH1D("h1","h1",100,0.2,0.7);
  hist->GetYaxis()->SetRangeUser(0.02,0.25);
  hist->GetXaxis()->SetRangeUser(0.2,0.7);
  hist->SetYTitle("Tkin(GeV)");
  hist->SetXTitle("<#beta^{}_{T}>");

  TGraph* g = new TGraph(5);

  g->SetMarkerStyle(20);

  for(rpy = 0; rpy < 5; rpy++){

    double temp = 0.0;
    temp = beta_T(beta_s[rpy],n[rpy]);

    g->SetPoint(rpy,temp,Tkin[rpy]);


  }

  hist->Draw("P");
  g->Draw("Psame");

  double ellipse_x[5][100];
  double ellipse_y[5][100];

  for(rpy = 0; rpy < 5; rpy++){

    double temp = beta_T(beta_s[rpy],n[rpy]);
    double shift = beta_s[rpy]-temp;

    for(int i = 0; i < 100; i++){

        test[rpy]->GetPoint(i,ellipse_x[rpy][i],ellipse_y[rpy][i]);
        ellipse_x[rpy][i] = ellipse_x[rpy][i] - shift;
    }   
  }

  TGraph* g1[5];

  for(rpy = 0; rpy < 5; rpy++){

    g1[rpy] = new TGraph(100);
      
      for(i = 0; i < 100; i++){
          
          g1[rpy]->SetPoint(i,ellipse_x[rpy][i],ellipse_y[rpy][i]);
          g1[rpy]->SetLineColor(kRed);
          g1[rpy]->SetMarkerColor(kRed);

      }

      g1[rpy]->Draw("same");
  }





}