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

vector<double> x,ex, y,ey;

double getChi2;
int getNdf;

double integral(double beta_s, double T, double pt, double mt)
{
  const double R = 1.;
  const double dr = 1e-2; // FIXME

  double s = 0.;
  for(double r = dr/2; r < R; r += dr)
  {
    double beta_r = beta_s * (r/R);
    double rho = TMath::ATanH(beta_r);

    s += r * dr * TMath::BesselK1(mt*cosh(rho)/T) * TMath::BesselI0(pt*sinh(rho)/T);
  }

  return s;
}

double MyFunc( double *pt, double *p){

  double mass = 0.497;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+ mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],pt[0],mt);

    return temp;
  
}

void BlastWaveTF1(){

  gStyle->SetErrorX(0);


	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy/AllMultbins_20pTBins_wSmooth10Eff_FullStats_v2_mt.root");
	
  TH1D* ksSpectra[9];
  TH1D* laSpectra[9];

  stringstream ksHistName;
  stringstream laHistName;

  for (int mult = 0; mult < 9; mult++){

    ksHistName.str("");
    laHistName.str("");

    ksHistName << "ksSpectra_";
    ksHistName << mult+1;

    laHistName << "laSpectra_";
    laHistName << mult+1;

    ksSpectra[mult] = (TH1D*)file->Get(ksHistName.str().c_str());
    laSpectra[mult] = (TH1D*)file->Get(laHistName.str().c_str());

  }
	

  double xmin = 0.2; 
  double xmax = 1.5;

  TCanvas* c1[9]; 

  TF1* f1 = new TF1("f1",MyFunc,xmin,xmax,3);

  f1->SetParameter(0,0.70);
  f1->SetParameter(1,0.15);
  f1->SetParameter(2,1000);

  f1->SetParLimits(0,0.1,1.0);
  f1->SetParLimits(1,0.05,1.0);
  f1->SetParLimits(2,1,1000000);

  TLatex* r3 = new TLatex(0.8,0.000001,"K^{0}_{s}");
  r3->SetTextSize(0.07);

  ksSpectra[0]->SetMarkerColor(kYellow-3);

  ksSpectra[1]->Scale(2);
  ksSpectra[1]->SetMarkerColor(1);
  ksSpectra[2]->Scale(4);
  ksSpectra[2]->SetMarkerColor(2);
  ksSpectra[3]->Scale(8);
  ksSpectra[3]->SetMarkerColor(3);
  ksSpectra[4]->Scale(16);
  ksSpectra[4]->SetMarkerColor(4);
  ksSpectra[5]->Scale(32);
  ksSpectra[5]->SetMarkerColor(5);
  ksSpectra[6]->Scale(64);
  ksSpectra[6]->SetMarkerColor(6);
  ksSpectra[7]->Scale(128);
  ksSpectra[7]->SetMarkerColor(7);
  ksSpectra[8]->Scale(256);
  ksSpectra[8]->SetMarkerColor(8);

  double Tkin[9];
  double beta_s[9];

  for(mult = 0; mult < 9; mult++){

    c1[mult] = new TCanvas();
    c1[mult]->SetLogy(1);

    ksSpectra[mult]->Fit("f1","Rsame");
    ksSpectra[mult]->GetXaxis()->SetRangeUser(0.0,8.0);
    ksSpectra[mult]->GetYaxis()->SetRangeUser(0.0000001,1000);
    ksSpectra[mult]->Draw("same");

    beta_s[mult] = f1->GetParameter(0);
    Tkin[mult] = f1->GetParameter(1);

  }

    ksSpectra[8]->SetMarkerStyle(20);
    ksSpectra[8]->SetStats(kFALSE);
    ksSpectra[8]->SetXTitle("m^{}_{T}-m^{}_{0} (GeV/c^{2})");
    ksSpectra[8]->SetYTitle("1/N^{}_{ev}1/(2#Pim^{}_{T}d^{2}N/(dm^{}_{T}dy) [(GeV/c)^{-2}]");

  for(mult = 0; mult < 8; mult++){

    ksSpectra[mult]->SetMarkerStyle(20);
    ksSpectra[mult]->SetStats(kFALSE);
    ksSpectra[mult]->Draw("Psame");

  }

  r3->Draw("same");

  TGraph* g1 = new TGraph(9);
  g1->GetXaxis()->SetTitle("beta_s");
  g1->GetYaxis()->SetTitle("Tkin (GeV)");
  for(mult = 0; mult < 9; mult++){

    g1->SetPoint(mult,beta_s[mult],Tkin[mult]);
    //g1->SetMarkerColor(mult+1);
    g1->SetMarkerStyle(20);

  }

  TCanvas* y1 = new TCanvas();
  g1->Draw("AP");
  

  /*cout << "beta_s: " << f1->GetParameter(0) << endl;
  cout << "T: " << f1->GetParameter(1) << endl;
  cout << "aK0: " << f1->GetParameter(2) << endl;
*/


}