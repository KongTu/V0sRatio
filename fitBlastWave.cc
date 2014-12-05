#include <iostream>
#include <fstream>
#include <cmath>

//#include "../interface/Enumerators.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"

using namespace std;

vector<double> x,ex, y,ey;
vector<int> pid;

char energy[256];
bool plot = false;

ofstream fileOut;

double getChi2;
int getNdf;

TF1 * I0;
TF1 * K1;

/*****************************************************************************/
double integral(double beta_s, double pt, double mt, double T)
{
  const double R = 1.;
  const double dr = 1e-2; // FIXME

  double s = 0.;
  for(double r = dr/2; r < R; r += dr)
  {
    double beta_r = beta_s * (r/R);
    double rho = atanh(beta_r);

    s += r * dr * K1->Eval(mt*cosh(rho)/T) * I0->Eval(pt*sinh(rho)/T);
  }

  return s;
}

/*****************************************************************************/
double getMass(int i)
{
  double m = 0;

  if(i == pip || i == pim) m = mass[pion];
  if(i == kap || i == kam) m = mass[kaon];
  if(i == prp || i == prm) m = mass[prot];
  if(i == elp || i == elm) m = mass[elec];

  return m;
}

/*****************************************************************************/
void function(int &npar, double *gin, double &f, double *par, int flag)
{
  double api    = par[0];
  double aka    = par[1];
  double aprp   = par[2];
  double aprm   = par[3];

  double T      = par[4];
  double beta_s = par[5];

  double chi2 = 0.;

  for(int i = 0; i < int(x.size()); i++)
  {
    double m = getMass(pid[i]);

    double pt = x[i];
    double mt = sqrt(m*m + pt*pt);

    double a = 0.;
    if(pid[i] == pip || pid[i] == pim) a = api;
    if(pid[i] == kap || pid[i] == kam) a = aka;
    if(pid[i] == prp) a = aprp;
    if(pid[i] == prm) a = aprm;
  
    double theo = a * pt * mt * integral(beta_s, pt, mt, T);

    if(plot)
    {
      if(i >= 1 && pid[i] != pid[i-1])
         fileOut << endl << endl;

      fileOut << " " << mt - m 
              << " " << y[i] << " " << ey[i] << " " << theo << endl;
    }
   
    double q = (y[i] - theo)/ey[i];

    chi2 += q*q;
  }

  f = chi2;

  getChi2 = f;
  getNdf  = x.size() - 6 - 1;
}

/*****************************************************************************/
void options(int arg, char **arc)
{
  int i = 1;

  do
  {
    if(strcmp(arc[i],"-energy") == 0) sprintf(energy,"%s",arc[++i]);

    i++;
  }
  while(i < arg);
}

/*****************************************************************************/
void fitBlastWave(int arg, char **arc)
{
  options(arg,arc);

  double mc = false;
  if(strcmp(energy,"AMPT")    == 0) mc = true;
  if(strcmp(energy,"EPOS")    == 0) mc = true;
  if(strcmp(energy,"HIJI")    == 0) mc = true;

  char name[256];
  if(!mc) sprintf(name, "../%s/res/blast_spectra.dat", energy);
     else sprintf(name, "../5TeV/res/blast_%s.dat", energy);

  fileOut.open(name);

  ofstream fileRes;

  if(!mc) sprintf(name, "../%s/res/blast_results.dat", energy);
     else sprintf(name, "../5TeV/res/blast_%s.res", energy);

  fileRes.open(name);

  int maxMult;
  if(strcmp(energy,"900GeV" ) == 0) maxMult =  7;
  if(strcmp(energy,"2760GeV") == 0) maxMult =  9;
  if(strcmp(energy,"7TeV")    == 0) maxMult = 12;
  if(strcmp(energy,"5TeV")    == 0) maxMult = 20;

  if(mc) maxMult = 0;

  cerr << " " << energy << endl;

  I0 = new TF1("besseI0","ROOT::Math::cyl_bessel_i(0.,x)",0.,10.);
  K1 = new TF1("besseK1","ROOT::Math::cyl_bessel_k(1.,x)",0.,10.);

  // Fit
  TVirtualFitter * func = TVirtualFitter::Fitter(0, 6);

  // Take all
  for(int mult = 0; mult <= maxMult; mult++)
  {
    cerr << "  mult: " << mult << " |";

    x.clear(); ex.clear();
    y.clear(); ey.clear();
    pid.clear();

    for(int i = pip; i <= prm; i++)
    {
      plot = false;
   
      if(mc)
      {
        sprintf(name,"../mc/dndpt_5TeV_%s.dat",energy);
        ifstream file(name);

        for(int j = pip; j <= prm; j++)
        for(int ipt = 0; ipt < 40; ipt++)
        {
          double pt, val, err;

          file >> pt >> val >> err;
           x.push_back(pt);
          ex.push_back(0.);

         // 1/mt dn/dmt 
           y.push_back(val / pt);
          ey.push_back(err / pt);

        pid.push_back(j);
        }
      }
      else
      {
      char name[256];

      if(mult == 0) sprintf(name,"../%s/res/dndpt_unknown.dat",energy);
               else sprintf(name,"../%s/res/dndpt_mult%d.dat" ,energy,mult);
      ifstream file(name);

      while(!file.eof()) 
      {
        double rap; file >> rap;

        if(!file.eof())
        {
        double pt ; file >> pt ;
  
        for(int j = pip; j <= elm; j++)
        {
          double val,sta,sym,syp;
          file >> val >> sta >> sym >> syp;
  
          if(i == j)
          if( (j != pip && j != pim) || pt > 0.4)
          if(val != 0.) 
          {
            double sys = (sym + syp)/2;
  
             x.push_back(pt);
            ex.push_back(0.);
  
            // 1/mt dn/dmt 
             y.push_back(val / pt);
            ey.push_back(sqrt(sta*sta + sys*sym) / pt); 
  
             pid.push_back(j);
          }
        }
        }
      }
  
      file.close();
      }
    }

    // Fit
    double arglist[100];

    arglist[0] = -1;  // verbosity
    func->ExecuteCommand("SET PRI",   arglist, 1);
    func->ExecuteCommand("CLEAR",     arglist, 0);
    func->ExecuteCommand("SET NOWAR", arglist, 0);

    func->SetFCN(function);

    int mult_;
    if(mult == 0) mult_ = maxMult/3;
             else mult_ = mult;

//    cerr << " mult_ = " << mult_ << endl;

    double api  = 20*mult_;
    double aka  = 15*mult_;
    double aprp = 80*mult_;
    double aprm = 80*mult_;

    func->SetParameter(0, "api",    api , 10.0, 1e+0, 1e+4);
    func->SetParameter(1, "aka",    aka , 10.0, 1e+0, 1e+4);
    func->SetParameter(2, "aprp",   aprp, 10.0, 1e+0, 1e+4);
    func->SetParameter(3, "aprm",   aprm, 10.0, 1e+0, 1e+4);

    func->SetParameter(4, "T",       0.15, 0.05, 0.05, 1.0);
    func->SetParameter(5, "beta_s",  0.70, 0.10, 0.05, 1.0);

    func->ExecuteCommand("SIMPLEX", arglist, 0);
    func->ExecuteCommand("HESSE"  , arglist, 0);

    fileRes << " " << mult
            << " " << func->GetParameter(4)
            << " " << func->GetParError(4)
            << " " << func->GetParameter(5)
            << " " << func->GetParError(5)
            << " " << getChi2 
            << " " << getNdf
            << " " << func->GetParameter(0)
            << " " << func->GetParameter(1)
            << " " << func->GetParameter(2)
            << " " << func->GetParameter(3)
            << endl;

    cerr    << " " << func->GetParameter(4)
            << " " << func->GetParError(4)
            << " " << func->GetParameter(5)
            << " " << func->GetParError(5);

    plot = true;
    func->ExecuteCommand("CALL FCN" , arglist, 0);

    cerr << " chi2/ndf = " << getChi2 << "/" << getNdf << endl;

    fileOut << endl << endl;
  }

  fileOut.close(); 
  fileRes.close(); 

  
}
