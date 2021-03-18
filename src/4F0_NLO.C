#include "Massless_lim.H"
#include "Massless_Kernels.H"
#include "PDF.H"

#include <cmath>
#include <iostream>


using namespace std;


double Massless_Limit::dsigma0NLO(double z[], size_t dim, void *p)
{
  double res(0.);
  // first round of variables
  double tau,xi,ximax,jac(1.),x,y;
  // lumis and MEs
  double lgg,lqg,lgq,lqq,
    sgg,sgq,sqg, sqq;
  //second variables
  double zz, tt;
  // double asr(2.*M_PI*Massless_Limit::p_lumi->ass(m2r));
  double acut = m_amin; //0.7e-2;
  double MAX  = 1. - acut;

  //  tau  = tauh + (1. - tauh)*z[0];
  //  jac *= (1. - tauh);
  tau = tauh*pow(MAX/tauh,z[0]);
  jac *= log(MAX/tauh)*tau;

  ximax = -1./2.*log(tau);
  xi    = -ximax + 2.*ximax*z[1];
  jac  *= 2.*ximax;

  x = sqrt(tau)*exp(xi);
  y = sqrt(tau)*exp(-xi);

  /// 5 <-> 0, as I somehow put an offset of 5
  lgg = Massless_Limit::p_lumi->Lumi(5,5,x,tauh/tau,1.,1.,10);
  lgq = Massless_Limit::L_gq(x, tauh/tau);
  lqg = Massless_Limit::L_gq(tauh/tau, x);
  lqq = Massless_Limit::L_qq(x,tauh/tau);

  //  zz = y + (1. - y)*z[2];
  //  jac *= (1. - y);
  zz  = y*pow(MAX/y,z[2]);
  double jyx = jac;
  double jz  = log(MAX/y)*zz;
  jac *= jz;


  // tt = y/zz + (MAX - y/zz)*z[3];
  // jac *= (MAX - y/zz);
  tt = y/zz*pow(MAX*zz/y,z[3]);
  jac *= log(MAX*zz/y)*tt;

  sgg  = dsigmagg_nlo(y,zz,tt);

  sgq  = dsigmagq_nlo(y,zz)/log(MAX*zz/y)/tt;

  sqg  = dsigmaqg_nlo(y,zz)/log(MAX*zz/y)/tt;

  sqq  = 0.;

  res  = jac*s0*(lgg*sgg + lqq*sqq
		 + lqg*sqg + lgq*sgq)/tau;

  return res;
}


double Massless_Limit::dsigmaqqb_nlo(double y, double z)
{
  return 0.;
}

double Massless_Limit::dsigmagg_nlo(double y, double z, double t)
{
  double res(0.);
  res = 2.*bb(y,z,t)+ 4.*bg(y,z)/log((1.-m_amin)*z/y)/t;
  double asr(masr(m2r));
  return pow(asr,3)*res;
}

double Massless_Limit::bb(double y, double z, double t){
  return bb_LO(y,z)/log((1.-m_amin)*z/y)/t + bb_NLO(y,z,t);
}

double Massless_Limit::bb_LO(double y, double z)
{
  double res(0.);
  double L(log(m2h/m2b)), agb, pqg;
  pqg = L*Massless_Limit::p_lumi->Pqg(y/z)/(2.*M_PI);
  agb = (L*L*Massless_Limit::p_lumi->agb(z,2,2)
         + L*Massless_Limit::p_lumi->agb(z,2,1)
         + Massless_Limit::p_lumi->agb(z,2,0))/pow(2.*M_PI,2);

  res = 2.*y*pqg*agb/z;

  return res;
}

double Massless_Limit::bb_NLO(double y, double z, double t)
{
  double res(0.);
  double sc(0.);
  double L(log(m2h/m2b)), pqg1, pqg2, pqg20;
  double k=y/z;
  double sbb_d, sbb_p, sbb_sub, sbb_rest;

  sbb_d = delta1_(&sc); sbb_p = dterms1_(&t,&sc);
  sbb_sub = -dtsub1_(&sc,&k); sbb_rest = sbbq1_(&t,&sc);

  pqg1  = Massless_Limit::p_lumi->Pqg(z)/(2.*M_PI)*L;
  pqg20 = Massless_Limit::p_lumi->Pqg(y/z)/(2.*M_PI)*L;

  pqg2  = Massless_Limit::p_lumi->Pqg(y/z/t)/(2.*M_PI)*L;

  // delta bit:
  res = (sbb_d - sbb_sub)*pqg20/log((1.-m_amin)*z/y)/t;

  res += y*sbb_p*(pqg2/t/t - pqg20);

  // rest:
  res += pqg2*sbb_rest/t/t;

  res *= y*pqg1/z;

  return res;
}

double Massless_Limit::bg(double y, double z)
{
  return (bg_LO(y,z) + bg_NLO(y,z));
}

double Massless_Limit::bg_LO(double y, double z)
{
  double res(0.);
  double L(log(m2h/m2b)), agb;
  double sbg, sc(0.);

  agb = (Massless_Limit::p_lumi->agb(y/z,2,2)*L*L
         + Massless_Limit::p_lumi->agb(y/z,2,1)*L
         + Massless_Limit::p_lumi->agb(y/z,2,0))/pow(2.*M_PI,2);

  sbg = sbg1_(&z,&sc);

  res = y*sbg*agb/z/z;

  return res;
}

double Massless_Limit::bg_NLO(double y, double z)
{
  double res(0.);
  double L(log(m2h/m2b)), pqg;
  double sbg, sc(0.),nf(5.);

  pqg = Massless_Limit::p_lumi->Pqg(y/z)*L/(2.*M_PI);
  sbg = (deltabga_(&z,&sc,&sc) +
         nf*deltabgf_(&z,&sc,&sc));

  res = y*pqg*sbg/z/z;

  return res;
}

double Massless_Limit::dsigmagq_nlo(double y, double z)
{
  double res(0.);
  res = bbs(y,z) + 2.*bgs(y,z) + 2.*bqs(y,z);
  double asr(masr(m2r));
  return pow(asr,3)*res;
}

double Massless_Limit::dsigmaqg_nlo(double y, double z)
{
  return Massless_Limit::dsigmagq_nlo(y,z);
}

double Massless_Limit::bbs(double y, double z)
{
  double res(0.);
  double L(log(m2h/m2b)), asb, pqg;
  pqg = Massless_Limit::p_lumi->Pqg(y/z)/(2.*M_PI)*L;
  asb = (Massless_Limit::p_lumi->asb(z,2,2)*L*L
         + Massless_Limit::p_lumi->asb(z,2,1)*L
         + Massless_Limit::p_lumi->asb(z,2,0))/pow(2.*M_PI,2);

  res = 2.*y*pqg*asb/z;

  return res;
}

double Massless_Limit::bgs(double y, double z)
{
  double res(0.);
  double L(log(m2h/m2b)), asb;
  double sbg, sc(0.);
  asb = (Massless_Limit::p_lumi->asb(y/z,2,2)*L*L
         + Massless_Limit::p_lumi->asb(y/z,2,1)*L
         + Massless_Limit::p_lumi->asb(y/z,2,0))/pow(2.*M_PI,2);
  sbg = sbg1_(&z,&sc);

  res = y*sbg*asb/z/z;

  return res;
}

double Massless_Limit::bqs(double y, double z)
{
  double res(0.);
  double L(log(m2h/m2b)), pqg;
  double sbq, sc(0.), nf(5);
  pqg = Massless_Limit::p_lumi->Pqg(y/z)/(2.*M_PI)*L;
  sbq = (deltabqa_(&z,&sc,&sc) +
         nf*deltabqf_(&z,&sc,&sc));

  res = y*pqg*sbq/z/z;

  return res;

}

double Massless_Limit::muFct(double z[], size_t dim, void *p)
{
  double res(0.);
  // first round of variables
  double tau,xi,ximax,jac(1.),x,y;
  // lumis and MEs
  double lgg,lqg,lgq,lqq,
    sgg,sgq,sqg, sqq;
  //second variables
  double zz, tt;
  // double asr(2.*M_PI*Massless_Limit::p_lumi->ass(m2r));
  double acut = m_amin; //0.7e-2;
  double MAX  = 1. - acut;

  //  tau  = tauh + (1. - tauh)*z[0];
  //  jac *= (1. - tauh);
  tau = tauh*pow(MAX/tauh,z[0]);
  jac *= log(MAX/tauh)*tau;

  ximax = -1./2.*log(tau);
  xi    = -ximax + 2.*ximax*z[1];
  jac  *= 2.*ximax;

  x = sqrt(tau)*exp(xi);
  y = sqrt(tau)*exp(-xi);

  p_lumi->SetmuF(m2h);
  /// 5 <-> 0, as I somehow put an offset of 5
  lgg = Massless_Limit::p_lumi->Lumi(5,5,x,tauh/tau,1.,1.,10);
  lgq = Massless_Limit::L_gq(x, tauh/tau);
  lqg = Massless_Limit::L_gq(tauh/tau, x);
  lqq = Massless_Limit::L_qq(x,tauh/tau);
  p_lumi->SetmuF(m2f);
  //  zz = y + (1. - y)*z[2];
  //  jac *= (1. - y);
  zz  = y*pow(MAX/y,z[2]);
  double jz  = log(MAX/y)*zz;

  sgg = slogg(y,zz,jz);
  sqg = slogq(y,zz,jz);
  sgq = sqg;
  sqq = sloqq(y,zz,jz);
  res  = jac*///
    (lgg*sgg + lqq*sqq + lqg*sqg + lgq*sgq)/tau;
  return res;
}


double Massless_Limit::slogg(double y, double z, double jac)
{
  double res(0.);
  double pgg(p_lumi->Pgg_pt(z)),
    pgg_plus(1./(1.-z)), pgg_0(p_lumi->Pgg_pt(1.)),
    pgg_delta(p_lumi->Pgg_d()), pgg_inte(-log(1.-y)),
    pgg_rest(p_lumi->Pgg(z));

  double delta(0.), plus(0.), rest(0.);
  delta = (pgg_delta - pgg_inte*pgg_0)*dsigmagg_lo(y);
  if(z<1.-1.e-5)
    plus  = jac*pgg_plus*(pgg*dsigmagg_lo(y/z) - pgg_0*dsigmagg_lo(y));
  rest = jac*pgg_rest*dsigmagg_lo(y/z);
  res = delta+plus+rest;
  return 2.*res;
}


double Massless_Limit::slogq(double y, double z, double jac)
{
  double res(0.);
  double pgq(p_lumi->Pgq(z)),pqg(p_lumi->Pqg(z)),nf(5.);
  res = jac*(nf*dsigmaqqb_lo(y/z)*pqg + dsigmagg_lo(y/z)*pgq );
  return 2.*res;
}

double Massless_Limit::sloqq(double y, double z,double jac)
{
  double res(0.);
  double cf=4./3.;
  double pqq(p_lumi->Pqq_pt(z)),
    pqq_plus(1./(1.-z)), pqq_0(p_lumi->Pqq_pt(1.)),
    pqq_delta(cf*3./2.), pqq_inte(-log(1.-y));

  double delta(0.), plus(0.), rest(0.);
  delta = (pqq_delta - pqq_inte*pqq_0)*dsigmaqqb_lo(y);
  if(z<1.-1.e-5)
    plus  = jac*pqq_plus*(pqq*dsigmaqqb_lo(y/z) - pqq_0*dsigmaqqb_lo(y));
  res = delta+plus;
  return 2.*res;
}

