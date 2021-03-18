#ifndef MASSLESS_LIM_H
#define MASSLESS_LIM_H

#include <vector>
#include <cmath>
#include <string>

#include "PDF.H"

class Scales {
public:
  int i_or, i_fonll; 
  double mbmb, mu_R, mu_F, m_H, ss;
  double m_vev, m_amin;
  int i_evs,  i_member;
public:
  Scales(const int order, const int ofonll,
	 const double mf, const double mur, const double muF, 
	 const double mH, const double S, const int evs,
	 const int PDFmem, const double amin)://, char *proc_name):
    i_or(order),i_fonll(ofonll),mbmb(mf),
    mu_R(mur),mu_F(muF),m_H(mH),ss(S),i_evs(evs),
    i_member(PDFmem)
    //    m_proc_name(proc_name)
  {
    m_amin = amin;
  }
};


class Massless_Limit {
private:
  LUMI * p_lumi;
  double lfr, lfh, tauh,m_inte,m_err;
  double m2r,m2f,m2h,m2b,S;
  double s0,m_amin,m_cu,m_cd;
  int i_nev, m_mode,m_proc;
  //m_mode = 0 ---> 5F or Massless_lim
  //m_mode = 1 ---> Difference Term
public:
  Massless_Limit(const char *pdfname, const Scales scale); 
  //const int mode);
  ~Massless_Limit();
  
  double dsigmaLO(double x[], size_t dim, void* p);
  double dsigmaNLO(double x[], size_t dim, 
		   void* p, const size_t plus);
  double dsigmaNNLO(double x[], size_t dim, 
		    void* p, const size_t plus);

  double dsigma0LO(double x[], size_t dim, void* p);
  double dsigma0NLO(double x[], size_t dim, void* p);

  double dsigmagg_lo(double y);
  
  double dsigmaqqb_lo(double y);
  double Agg_2_2(double y);
  double Agg_2_1(double y);
  double Agg_2_0(double y);
  
  double dsigmagg_nlo(double y,double z, double t);
  double dsigmaqqb_nlo(double y, double z);
  double dsigmagq_nlo(double y, double z);
  double dsigmaqg_nlo(double y, double z);
 

  double muFct(double z[], size_t dim, void* p);
  double slogg(double y,double z, double jz);
  double sloqq(double y, double z, double jz);
  double slogq(double y, double z, double jz);

  double bb(double y, double z, double t);
  double bg(double y, double z);
  double bb_LO(double y, double z);
  double bb_NLO(double y, double z, double t);
  double bg_LO(double y, double z);
  double bg_NLO(double y, double z);
  
  double bbs(double y, double z);
  double bgs(double y, double z);
  double bqs(double y, double z);

  double L_qq(const double x, const double y);
  double L_gq(const double x, const double y);
  inline double masr(const double Q2) { return 2.*M_PI*p_lumi->ass(Q2); }
};

double wrap_0LO(double x[], size_t dim, void* p);

double wrap_0NLO(double x[], size_t dim, void* p);

double wrap_muFTerms(double x[], size_t dim, void* p);

double wrap_LO(double x[], size_t dim, void* p);

double wrap_NLO(double x[], size_t dim, void* p);

double wrap_NNLO(double x[], size_t dim, void* p);

double xsec5F_LO(const Scales scal);

double xsec5F_NLO(const Scales scal);

double xsec5F_NNLO(const Scales scal);

double xsec4F0_LO(const Scales scal);

double xsec4F0_NLO(const Scales scal);

double muFTerm(const Scales scal);

double xsec5F(const char *argv, const Scales scal);

double xsec4F0(const char *argv, const Scales scal);

double xsecDiff(const char *argv, const Scales scal);


#endif
