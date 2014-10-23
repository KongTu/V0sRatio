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

    s += r * dr * TMath::BesselK1((mt*cosh(rho))/T) * TMath::BesselI0((pt*sinh(rho))/T);
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


double beta_T(double beta_s, double n){

  double temp = 0.;
  double s = 0.;

    temp = 2*3.14*beta_s;

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
  
  double yields_1[] = { 1.532, 1.373, 1.22, 1.026, 0.8654, 0.7002, 0.5571, 0.4442, 0.3476, 
    0.2745, 0.2161, 0.1731, 0.1342, 0.1068, 0.0867, 0.06848, 0.0553, 0.04467, 0.03537, 
    0.02874, 0.02089, 0.01407, 0.009662, 0.006504, 0.004485, 0.002855, 0.001833, 0.00112, 7.169E-4, 
    4.67E-4, 2.456E-4, 1.467E-4, 9.063E-5, 2.868E-5 };

  double yields_2[] = { 1.282, 1.125, 0.9841, 0.8456, 0.6969, 0.5672, 0.447, 0.3563, 0.2742, 
    0.2176, 0.1696, 0.1349, 0.1055, 0.08433, 0.06682, 0.05295, 0.04317, 0.03504, 0.02797, 
    0.02255, 0.01616, 0.01112, 0.007565, 0.005198, 0.00362, 0.002319, 0.001436, 9.019E-4, 5.672E-4, 
    3.625E-4, 2.299E-4, 1.311E-4, 6.245E-5, 2.493E-5 };

  double yields_3[] = { 1.032, 0.9984, 0.8801, 0.7261, 0.5949, 0.4825, 0.3799, 0.2988, 0.2292, 
    0.1805, 0.141, 0.1106, 0.0862, 0.06891, 0.05424, 0.04371, 0.03475, 0.02833, 0.02264, 
    0.01838, 0.01337, 0.009004, 0.006065, 0.004339, 0.003005, 0.001894, 0.001185, 7.542E-4, 4.686E-4, 
    3.037E-4, 1.802E-4, 1.081E-4, 6.501E-5, 2.173E-5 };

  double yields_4[] = { 0.8664, 0.8167, 0.7138, 0.5872, 0.4752, 0.3726, 0.2879, 0.2252, 0.172, 
    0.1339, 0.1043, 0.0819, 0.06291, 0.04975, 0.03974, 0.03109, 0.02512, 0.02004, 0.01625, 
    0.01301, 0.009402, 0.006303, 0.004398, 0.003022, 0.002128, 0.001354, 8.566E-4, 5.426E-4, 3.412E-4, 
    2.217E-4, 1.309E-4, 7.942E-5, 4.63E-5, 1.598E-5 };

  double yields_5[] = { 0.6602, 0.6161, 0.529, 0.4355, 0.3441, 0.2617, 0.2011, 0.1554, 0.1174, 
    0.08973, 0.06793, 0.05317, 0.04105, 0.03174, 0.0248, 0.01966, 0.01595, 0.01258, 0.009976, 
    0.008011, 0.005954, 0.003914, 0.002692, 0.001828, 0.001312, 8.455E-4, 5.218E-4, 3.186E-4, 2.097E-4, 
    1.396E-4, 8.241E-5, 4.926E-5, 2.905E-5, 9.493E-6 };

  double yields_6[] = { 0.5099, 0.4278, 0.3651, 0.2804, 0.2207, 0.1663, 0.1246, 0.09248, 0.06821, 
    0.0513, 0.03895, 0.02908, 0.02266, 0.01756, 0.01333, 0.01052, 0.008219, 0.006613, 0.005288, 
    0.004173, 0.003059, 0.002031, 0.001332, 9.469E-4, 6.523E-4, 4.37E-4, 2.7E-4, 1.664E-4, 1.1E-4, 
    6.73E-5, 4.216E-5, 2.427E-5, 1.468E-5, 4.802E-6 };

  double yields_7[] = { 0.2596, 0.2172, 0.1874, 0.1448, 0.1059, 0.07555, 0.05461, 0.03886, 0.0274, 
    0.02033, 0.01438, 0.01057, 0.008109, 0.006127, 0.004504, 0.003565, 0.002764, 0.002139, 0.001738, 
    0.001358, 9.955E-4, 6.173E-4, 4.244E-4, 2.994E-4, 1.98E-4, 1.296E-4, 7.599E-5, 5.017E-5, 2.996E-5, 
    1.953E-5, 1.094E-5, 8.026E-6, 4.519E-6, 1.483E-6 };

  ksSpectra[0] = new TH1D("ksSpectra_1","ksSpectra_1",34,ptbins);
  ksSpectra[1] = new TH1D("ksSpectra_2","ksSpectra_2",34,ptbins);
  ksSpectra[2] = new TH1D("ksSpectra_3","ksSpectra_3",34,ptbins);
  ksSpectra[3] = new TH1D("ksSpectra_4","ksSpectra_4",34,ptbins);
  ksSpectra[4] = new TH1D("ksSpectra_5","ksSpectra_5",34,ptbins);
  ksSpectra[5] = new TH1D("ksSpectra_6","ksSpectra_6",34,ptbins);
  ksSpectra[6] = new TH1D("ksSpectra_7","ksSpectra_7",34,ptbins);

  for(int pt = 0; pt < 34; pt++){
  
    ksSpectra[0]->SetBinContent(pt+1,yields_1[pt]);
    ksSpectra[0]->SetBinError(pt+1,sqrt(yields_1[pt]));   
    ksSpectra[1]->SetBinContent(pt+1,yields_2[pt]);
    ksSpectra[1]->SetBinError(pt+1,sqrt(yields_2[pt]));
    ksSpectra[2]->SetBinContent(pt+1,yields_3[pt]);
    ksSpectra[2]->SetBinError(pt+1,sqrt(yields_3[pt]));
    ksSpectra[3]->SetBinContent(pt+1,yields_4[pt]);
    ksSpectra[3]->SetBinError(pt+1,sqrt(yields_4[pt]));
    ksSpectra[4]->SetBinContent(pt+1,yields_5[pt]);
    ksSpectra[4]->SetBinError(pt+1,sqrt(yields_5[pt]));
    ksSpectra[5]->SetBinContent(pt+1,yields_6[pt]);
    ksSpectra[5]->SetBinError(pt+1,sqrt(yields_6[pt]));
    ksSpectra[6]->SetBinContent(pt+1,yields_7[pt]);
    ksSpectra[6]->SetBinError(pt+1,sqrt(yields_7[pt]));

  }

  double xmin = 0.0; 
  double xmax = 1.5;

  TF1* f1 = new TF1("f1",MyFunc,xmin,xmax,4);

  f1->SetParameter(0,0.8);
  f1->SetParameter(1,0.14);
  f1->SetParameter(2,1000);
  f1->SetParameter(3,1);

  f1->SetParLimits(0,0.5,1.0);
  f1->SetParLimits(1,0.05,0.25);
  f1->SetParLimits(2,1,10000);
  f1->SetParLimits(3,0.5,5);

  ksSpectra[0]->SetMarkerColor(kYellow-3);
  ksSpectra[0]->SetMarkerSize(1.3);
  ksSpectra[0]->Scale(64);

  ksSpectra[1]->Scale(32);
  ksSpectra[1]->SetMarkerColor(1);
  ksSpectra[1]->SetMarkerSize(1.3);
  ksSpectra[2]->Scale(16);
  ksSpectra[2]->SetMarkerColor(2);
  ksSpectra[2]->SetMarkerSize(1.3);
  ksSpectra[3]->Scale(8);
  ksSpectra[3]->SetMarkerColor(3);
  ksSpectra[3]->SetMarkerSize(1.3);
  ksSpectra[4]->Scale(4);
  ksSpectra[4]->SetMarkerColor(4);
  ksSpectra[4]->SetMarkerSize(1.3);
  ksSpectra[5]->Scale(2);
  ksSpectra[5]->SetMarkerColor(5);
  ksSpectra[5]->SetMarkerSize(1.3);
  ksSpectra[6]->Scale(1);
  ksSpectra[6]->SetMarkerColor(6);
  ksSpectra[6]->SetMarkerSize(1.3);

  double Tkin[7];
  double beta_s[7];
  double n[7];
  double beta_T[7];

  TCanvas* c1[7];

  for(int mult = 0; mult < 7; mult++){

    c1[mult] = new TCanvas();
    c1[mult]->SetLogy(1);

    ksSpectra[mult]->Fit("f1","LLRsame");
    ksSpectra[mult]->GetXaxis()->SetRangeUser(0.0,8.0);
    ksSpectra[mult]->GetYaxis()->SetRangeUser(0.0000001,1000);
    ksSpectra[mult]->Draw("same");

    beta_s[mult] = f1->GetParameter(0);
    Tkin[mult] = f1->GetParameter(1);
    n[mult] = f1->GetParameter(3);

  }

  cout << "STATUS: " << gMinuit->fStatus << endl;
  cout << "fCstatus:" << gMinuit->fCstatu << endl;

  ksSpectra[6]->SetMarkerStyle(20);
  ksSpectra[6]->SetStats(kFALSE);
  ksSpectra[6]->SetXTitle("m^{}_{T}-m^{}_{0} (GeV/c^{2})");
  ksSpectra[6]->SetYTitle("1/N^{}_{ev}1/(2#Pim^{}_{T}d^{2}N/(dm^{}_{T}dy) [(GeV/c)^{-2}]");
  ksSpectra[6]->Draw("Psame");

  for(mult = 0; mult < 6; mult++){

    ksSpectra[mult]->SetMarkerStyle(20);
    ksSpectra[mult]->SetStats(kFALSE);
    ksSpectra[mult]->Draw("Psame");

  }

  /**
 * Below is to calculate the <beta_T> and the plot of Tkin vs <beta_T>
 */

  for(mult = 0; mult < 7; mult++){

    beta_T[mult] = beta_T(beta_s[mult],n[mult]);

  }

  TGraph* g1 = new TGraph(7);

  for(mult = 0; mult < 7; mult++){

    g1->SetPoint(mult,beta_s[mult],Tkin[mult]);
    //g1->SetMarkerColor(mult+1);
    g1->SetMarkerStyle(20);

    cout << "Tkin: " << Tkin[mult] << endl;
    cout << "n: " << n[mult] << endl;

  }

  TCanvas* y1 = new TCanvas();
  g1->Draw("AP");




 
	



}