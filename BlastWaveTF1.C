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

const double R = 1.0;
const double dr = 1e-3; // FIXME

double integral(double beta_s, double T, double n, double pt, double mt)
{
  
  double s = 0.;
  for(double r = dr/2; r < R; r += dr)
  {
    double beta_r = beta_s * TMath::Power((r/R),n);
    double rho = TMath::ATanH(beta_r);

    s += r * dr * TMath::BesselK1(mt*cosh(rho)/T) * TMath::BesselI0(pt*sinh(rho)/T);
  }

  return s;
}

double MyFunc( double *pt, double *p){

  double mass = 0.497;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}

//fit mt-m0:
/*
double MyFunc( double *x, double *p){

  double mass = 0.497;

  double temp = 0.;
    double mt = (x[0]+mass);
    double pt = sqrt(mt*mt - mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt,mt);

    return temp;
  
}*/


double beta_T(double beta_s, double n){

  double temp = 0.;
  double s = 0.;

    temp = 2 * 3.14 * beta_s;

    for(double r = dr/2; r < R; r += dr)
  {
    double value = TMath::Power((r/R),n);

    s += r * dr * value;
  }

 return (s * temp);

}

void BlastWaveTF1(){

  gStyle->SetErrorX(0);

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy/AllMultbins_23ks_pTBins_wSmooth10Eff_FullStats_v1_pt.root");
	
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
	

  double xmin = 0.5; 
  double xmax = 6.0;

  TCanvas* c1[9]; 

  TF1* f1 = new TF1("f1",MyFunc,xmin,xmax,4);

  f1->SetParameter(0,0.7);
  f1->SetParameter(1,0.15);
  f1->SetParameter(2,1000);
  f1->SetParameter(3,1);

  f1->SetParLimits(0,0.15,1.0);
  f1->SetParLimits(1,0.05,1.0);
  f1->SetParLimits(2,1,10000);
  f1->SetParLimits(3,0,10.0);

  TLatex* r1 = new TLatex(3.7,152,"CMS,p-Pb 2013,#sqrt{S^{}_{NN}} = 5.02 TeV");
  r1->SetTextSize(0.04);
  TLatex* r5 = new TLatex(3.7,24,"Blast Wave Fit (individual fits)");
  r5->SetTextSize(0.04);
  TLatex* r2 = new TLatex(6,5,"-2.4 < y < 2.4");
  r2->SetTextSize(0.04);
  TLatex* r3 = new TLatex(0.8,0.000001,"K^{0}_{s}");
  r3->SetTextSize(0.07);
  TLatex* r4 = new TLatex(0.8,0.000001,"#Lambda+#bar{#Lambda}");

  ksSpectra[0]->SetMarkerColor(kYellow-3);
  ksSpectra[0]->SetMarkerSize(1.3);

  ksSpectra[1]->Scale(2);
  ksSpectra[1]->SetMarkerColor(1);
  ksSpectra[1]->SetMarkerSize(1.3);
  ksSpectra[2]->Scale(4);
  ksSpectra[2]->SetMarkerColor(2);
  ksSpectra[2]->SetMarkerSize(1.3);
  ksSpectra[3]->Scale(8);
  ksSpectra[3]->SetMarkerColor(3);
  ksSpectra[3]->SetMarkerSize(1.3);
  ksSpectra[4]->Scale(16);
  ksSpectra[4]->SetMarkerColor(4);
  ksSpectra[4]->SetMarkerSize(1.3);
  ksSpectra[5]->Scale(32);
  ksSpectra[5]->SetMarkerColor(5);
  ksSpectra[5]->SetMarkerSize(1.3);
  ksSpectra[6]->Scale(64);
  ksSpectra[6]->SetMarkerColor(6);
  ksSpectra[6]->SetMarkerSize(1.3);
  ksSpectra[7]->Scale(128);
  ksSpectra[7]->SetMarkerColor(7);
  ksSpectra[7]->SetMarkerSize(1.3);
  ksSpectra[8]->Scale(256);
  ksSpectra[8]->SetMarkerColor(8);
  ksSpectra[8]->SetMarkerSize(1.3);

  double Tkin[9];
  double beta_s[9];
  double n[9];
  double beta_T[9];

  for(mult = 0; mult < 9; mult++){

    c1[mult] = new TCanvas();
    c1[mult]->SetLogy(1);

    ksSpectra[mult]->Fit("f1","R");
    ksSpectra[mult]->GetXaxis()->SetRangeUser(0.0,8.0);
    ksSpectra[mult]->GetYaxis()->SetRangeUser(0.0000001,1000);
    ksSpectra[mult]->Draw("Psame");

    beta_s[mult] = f1->GetParameter(0);
    Tkin[mult] = f1->GetParameter(1);
    n[mult] = f1->GetParameter(3);

  }

    ksSpectra[8]->SetMarkerStyle(20);
    ksSpectra[8]->SetStats(kFALSE);
    ksSpectra[8]->SetXTitle("P^{}_{T}(GeV/c)");
    ksSpectra[8]->SetYTitle("1/N^{}_{ev}1/(2#Pip^{}_{T}d^{2}N/(dp^{}_{T}dy) [(GeV/c)^{-2}]");

  for(mult = 0; mult < 8; mult++){

    ksSpectra[mult]->SetMarkerStyle(20);
    ksSpectra[mult]->SetStats(kFALSE);
    ksSpectra[mult]->Draw("Psame");

  }

  r1->Draw("same");
  r2->Draw("same");
  r3->Draw("same");
  r5->Draw("same");

/**
 * Below is to calculate the <beta_T> and the plot of Tkin vs <beta_T>
 */

  for(mult = 0; mult < 9; mult++){

    beta_T[mult] = beta_T(beta_s[mult],n[mult]);

  }



  TGraph* g1[9];

  for(mult = 0; mult < 9; mult++){

    g1[mult] = new TGraph(1);

    g1[mult]->SetPoint(mult,beta_s[mult],Tkin[mult]);
    g1[mult]->SetMarkerStyle(20);
    g1[mult]->SetMarkerColor(mult);

  }

  TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->AddEntry(ksSpectra[0],"0  < N^{offline}_{trk} < 35");
    w1->AddEntry(ksSpectra[1],"35  < N^{offline}_{trk} < 60");
    w1->AddEntry(ksSpectra[2],"60  < N^{offline}_{trk} < 90");
    w1->AddEntry(ksSpectra[3],"90  < N^{offline}_{trk} < 120");
    w1->AddEntry(ksSpectra[4],"120 < N^{offline}_{trk} < 150");
    w1->AddEntry(ksSpectra[5],"150 < N^{offline}_{trk} < 185");
    w1->AddEntry(ksSpectra[6],"185 < N^{offline}_{trk} < 220");
    w1->AddEntry(ksSpectra[7],"220 < N^{offline}_{trk} < 260");
    w1->AddEntry(ksSpectra[8],"     N^{offline}_{trk} > 260");


  TCanvas* y1 = new TCanvas();
  g1[0]->GetXaxis()->SetRangeUser(0.5,1.5);
  g1[0]->GetYaxis()->SetRangeUser(0,0.4);
  g1[0]->SetMarkerColor(kYellow-3);
  g1[0]->Draw("AP");

  for(mult = 0; mult < 9; mult++){

    g1[mult]->Draw("Psame");
  }
  
  w1->Draw("same");

  /*cout << "beta_s: " << f1->GetParameter(0) << endl;
  cout << "T: " << f1->GetParameter(1) << endl;
  cout << "aK0: " << f1->GetParameter(2) << endl;
*/


}