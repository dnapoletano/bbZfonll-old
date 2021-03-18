////
// Calculate the 4F at LO and NLO both in the matched and in the standard scheme
// NB: need a file with a list of numbers (to be changed accordingly with the LHCHXSWG requests
////

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdexcept>

#include "4F.H"
#include "Massless_lim.H"


using namespace std;

LUMI * lumi_instance;
//double * temp;

void ReadFile(const char * filename, const size_t length, double vec[])
{
  ifstream infile;
  int k = 0;
  infile.open(filename);
  if (infile.fail())
    {
      cout << "error: couldn't open " << filename << endl;
    }
  while (!infile.eof())
    {
      infile >> vec[k];
      ++k;
    }
   infile.close();
}

const int len = 11;

double CalcggLO(const int pos, const int ss)
{
  double xsecggLO[len];
  if (ss == 1) // central value
    ReadFile("XS4F/4F_gg_central.dat",len, xsecggLO);
  if (ss == 2) // min scale
    ReadFile("XS4F/4F_gg_min.dat",len, xsecggLO);
  if (ss == 3) // max scale
    ReadFile("XS4F/4F_gg_max.dat",len, xsecggLO);
  return xsecggLO[pos];
}

double CalcLO(const int pos,const int ss)
{
  double xsecLO[len];
  if (ss == 1) // central value
    ReadFile("XS4F/4F_LO_central.dat",len,xsecLO);
  if (ss == 2) // min scale
    ReadFile("XS4F/4F_LO_min.dat",len,xsecLO);
  if (ss == 3) // max scale
    ReadFile("XS4F/4F_LO_max.dat",len,xsecLO);
  return xsecLO[pos];
}

double CalcNLO(const int pos,const int ss)
{
  double xsecNLO[len];
  if (ss == 1) // central value
    ReadFile("XS4F/4F_NLO_central.dat",len,xsecNLO);
  if (ss == 2) // min scale
    ReadFile("XS4F/4F_NLO_min.dat",len,xsecNLO);
  if (ss == 3) // max scale
    ReadFile("XS4F/4F_NLO_max.dat",len,xsecNLO);
  return xsecNLO[pos];
}

double Calc4F_mb(const int pos)
{
  double xsecmb[21];
  ReadFile("XS4F/XSmb.dat",21,xsecmb);
  
  return xsecmb[pos];
}

double Calc4F0_mb(const double mb, const int order)
{
  Scales sd(-1,order,mb,125.,125.,125.,13000.,1000,0,1.e-7);
  double a = xsecDiff("NNPDF30_nnlo_as_0118",sd); 
  return a;
}

double Calc4FLO(const int pos)
{
  double xsecLO[len];
  ReadFile("XS4F/4F_Only_LO.dat",len,xsecLO);
  return xsecLO[pos];
}

double Calc4FNLO(const int pos)
{
  double xsecNLO[len];
  ReadFile("XS4F/4F_Only_NLO.dat",len,xsecNLO);
  return xsecNLO[pos];
}

double kappa(const double x)
{
  double Tr(0.5);
  return 4.*Tr*x/(6.*M_PI);
}

int findpos(const double muR, const double muRmin, 
	    const double muRmax, const int Npoints)
{
  int pos = 0;
  pos = (int)ceil((muR-muRmin)*Npoints/(muRmax-muRmin));
  if(pos-1<0) return 0;
  else return pos-1;
}

double b0(const double nf)
{
  return (33. - 2.*nf)/(12.*M_PI);
}

double match()
{
  return 0.;
}

double Calc_4F_matched(const char* argv, const Scales scale, 
		       const int o_fonll, const int ss)
{
  double xs(0.), muR(scale.mu_R), 
    muF(scale.mu_F),m2b, mh;
  if(scale.i_or == -1)
    lumi_instance = new LUMI(argv,0,muF,scale.mbmb, true);
  else
    lumi_instance = new LUMI(argv,0,muF,4.92, false);
  m2b = lumi_instance->m2bpdf();
  double alphasR(2.*M_PI*lumi_instance->ass(muR*muR));
  double asF(lumi_instance->ass(muF*muF)); //asF is alphaS/2pi
  mh = scale.m_H;
  int pos = findpos(mh, 25., 525., len);

  if(o_fonll == 1) { // FONLL-A
    xs = CalcLO(pos, ss);
  }

  if(o_fonll == 2) { // FONLL-B
    double LR,LF;
    LR = log(muR*muR/m2b); LF = log(muF*muF/m2b);
    xs = CalcNLO(pos, ss);
    xs+= (alphasR*kappa(LR)
	  +2.*alphasR*match()*LR )*CalcLO(pos, ss)
      - asF*kappa(LF)*CalcggLO(pos, ss); 
  }
  return xs;
}
