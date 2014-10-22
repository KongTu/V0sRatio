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
const double dr = 1e-2; // FIXME

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
    double mt = sqrt(pt[0]*pt[0] + mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}


double beta_T(double beta_s, double n){

  double temp = 0.;
  double s = 0.;

    temp = beta_s;

    for(double r = dr/2; r < R; r += dr)
  {
    double value = TMath::Power((r/R),n);

    s += r * dr * value;
  }

 return (s * temp);

}

void ALICEblastwaveFit(){

  gStyle->SetErrorX(0);

  TH1D* ksSpectra[9];
  TH1D* laSpectra[9];

  double ptbins[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.3,3.6,3.9,4.2,4.6,5.0,5.5,6.0,8.0};
  double yields[] = { 0.8664, 0.8167, 0.7138, 0.5872, 0.4752, 0.3726, 0.2879, 0.2252, 0.172, 
    0.1339, 0.1043, 0.0819, 0.06291, 0.04975, 0.03974, 0.03109, 0.02512, 0.02004, 0.01625, 
    0.01301, 0.009402, 0.006303, 0.004398, 0.003022, 0.002128, 0.001354, 8.566E-4, 5.426E-4, 3.412E-4, 
    2.217E-4, 1.309E-4, 7.942E-5, 4.63E-5, 1.598E-5 };

  ksSpectra[0] = new TH1D("ksSpectra_1","ksSpectra_1",34,ptbins);

  for(int pt = 0; pt < 34; pt++){
  
    ksSpectra[0]->SetBinContent(pt+1,yields[pt]);

  }

  double xmin = 0.00; 
  double xmax = 1.50;

  TF1* f1 = new TF1("f1",MyFunc,xmin,xmax,4);

  f1->SetParameter(0,0.70);
  f1->SetParameter(1,0.15);
  f1->SetParameter(2,1000);
  f1->SetParameter(3,7);

  f1->SetParLimits(0,0.15,1.0);
  f1->SetParLimits(1,0.05,1.0);
  f1->SetParLimits(2,1,1000000);
  f1->SetParLimits(3,1,20);

  TCanvas* c1 = new TCanvas();
  c1->SetLogy(1);

  ksSpectra[0]->Fit("f1","Rsame");
  ksSpectra[0]->GetXaxis()->SetRangeUser(0.0,8.0);
  ksSpectra[0]->GetYaxis()->SetRangeUser(0.0000001,1000);
  ksSpectra[0]->SetMarkerStyle(20);
  ksSpectra[0]->SetMarkerSize(1.3);

  ksSpectra[0]->Draw("Psame");

  double Tkin = f1->GetParameter(1);
  double beta_s = f1->GetParameter(0);
  double n = f1->GetParameter(3);


  cout << "Tkin: " << f1->GetParameter(1) << endl;
  cout << "beta_s: " << f1->GetParameter(0) << endl;
  cout << "n: " << f1->GetParameter(3) << endl;
  cout << "beta_T: " << beta_T(beta_s,n) << endl;





 
	



}