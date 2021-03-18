///
// Calculate the massless limit at LO and NLO
///

#include "Massless_lim.H"
#include "Massless_Kernels.H"

#include <iostream>
#include <vector>
#include <cmath>
using namespace LHAPDF;
using namespace std;

double Massless_Limit::dsigmagg_lo(double y)
{
  double ds(0.), L(log(m2h/m2b)/(2.*M_PI));
  double asr(2.*M_PI*p_lumi->ass(m2r));
  ds = asr*asr*(Agg_2_2(y)*L*L + Agg_2_1(y)*L + Agg_2_0(y)); // Eq (29) of FONLL-A paper
  return ds;                                                        // A20 is not convoluted with anything, z variable of convolution
}

double Massless_Limit::dsigmaqqb_lo(double y)
{
  double ds(0.);
  double asr(2.*M_PI*p_lumi->ass(m2r));
  double mylfh=0., mylfr=0.;
  ds = asr*asr*deltaqqba_(&y,&mylfh,&mylfr); 
  return ds;
}

double Massless_Limit::Agg_2_0(double y)
{
  // constant in z!
  double mylfh=0., mylfr=0.;
  return deltagga_(&y,&mylfh,&mylfr);
}

double Massless_Limit::Agg_2_1(double y)
{

  return (y*pow(M_PI,-1)*(-3*(-1 + y)*(7 + 67*y) + 48*(-1 + y)*(1 + 3*y)*log(1 - y) + 
			  24*dilog(y)*pow(1 + 2*y,2) + 
			  6*log(y)*(1 - 4*(-2 + y)*y + log(y)*pow(1 + 2*y,2)) -
			  4*pow(M_PI + 2*M_PI*y,2)))/24.;
}

double Massless_Limit::Agg_2_2(double y)
{
  return (y*(-2 - 4*y + 6*pow(y,2) - log(y)*pow(1 + 2*y,2)))/2.;
}

double Massless_Limit::L_qq(const double x, const double y)
{
  double lqq(0.);
  lqq +=  p_lumi->Lumi(4,6,x,y,1.,1.,10);//d
  lqq +=  p_lumi->Lumi(6,4,x,y,1.,1.,10);
  lqq +=  p_lumi->Lumi(3,7,x,y,1.,1.,10);//u
  lqq +=  p_lumi->Lumi(7,3,x,y,1.,1.,10);
  lqq +=  p_lumi->Lumi(2,8,x,y,1.,1.,10);//s
  lqq +=  p_lumi->Lumi(8,2,x,y,1.,1.,10);
  lqq +=  p_lumi->Lumi(1,9,x,y,1.,1.,10);//c
  lqq +=  p_lumi->Lumi(9,1,x,y,1.,1.,10);
  return lqq;
  
}

double Massless_Limit::L_gq(const double x, const double y)
{
  double lgq(0.);
  lgq +=  p_lumi->Lumi(5,6,x,y,1.,1.,10);
  lgq +=  p_lumi->Lumi(5,4,x,y,1.,1.,10);
  lgq +=  p_lumi->Lumi(5,7,x,y,1.,1.,10);
  lgq +=  p_lumi->Lumi(5,3,x,y,1.,1.,10);
  lgq +=  p_lumi->Lumi(5,8,x,y,1.,1.,10);
  lgq +=  p_lumi->Lumi(5,2,x,y,1.,1.,10);
  lgq +=  p_lumi->Lumi(5,9,x,y,1.,1.,10);
  lgq +=  p_lumi->Lumi(5,1,x,y,1.,1.,10);
  return lgq;
  
}

double Massless_Limit::dsigma0LO(double z[], size_t dim, void* p)
{
  double dsg(0.), x,y;
  double tau, jac(1.), xi;
  double t, tmin;

  double lgg,lqq;
  double sgg(0.), sqq(0.);

  tau  = tauh + (1. - tauh)*z[0];
  jac *= (1. - tauh);
  double ximax = - 0.5*log(tau);
  xi   = -ximax + 2.*ximax*z[1];
  jac *=  2.*ximax;

  //standard variables  x, y ;
  x = sqrt(tau)*exp(xi);
  y = sqrt(tau)*exp(-xi);
  
  // gg
  sgg        =       dsigmagg_lo(y);
  // q q~
  sqq        =       dsigmaqqb_lo(y);
  
  lgg     = p_lumi->Lumi(5,5,x,tauh/tau,1.,1.,10);
  lqq     = L_qq(x,tauh/tau);

  dsg = jac*s0*( lqq * sqq + lgg * sgg )/tau;
  return dsg;
}



