#ifndef PDF_H
#define PDF_H

#include "LHAPDF/LHAPDF.h"
#include <cmath>
#include <vector>


class LUMI {
private :

  LHAPDF::PDF *m_pdf;
  double m_as,Q,m2b;
  double mm2b;
  bool m_flag, m_lim;

public :  
  LUMI(const char* setname, const int setmem, const double muF,
       const double mbext, const bool ismlim);
  ~LUMI();
  double Lumi(const int i, const int j, 
	      const double x, const double y,
	      const double z1, const double z2, 
	      const int order);
  inline double alphaS()             { return m_as;  }
  inline double ass(const double Q2) { return m_pdf->alphasQ(sqrt(Q2))/(2.*M_PI); }
  double mb(const double mf, const double Q,const int nloop);
  inline double m2bpdf() { return m2b; }
  inline double m2bext() {return mm2b; }
  
  double Pqg(const double x);
  
  double Pgq(const double x);

  double Pgg(const double x); // finite bits
  double Pgg_d(); // delta term P_gg
  double Pgg_plus(const double x); // plus distribution
  double Pgg_pt(const double x); // plus coefficient
  double Pgg_psub(const double x); // plus reminder
  
  double Pqq_plus(const double x); // same as above
  double Pqq_pt(const double x);
  double Pqq_psub(const double x);
  double Pqq_d();
  
  double phi(const double x, const double z, 
	     const double Q2, const int o);
  double agb(const double z, const int fo, const int lo);
  double asb(const double z, const int fo, const int lo);
  void SetmuF(const double muf)  {Q=muf;}
    
  bool Init_Lumi() { return m_flag;}
};

#endif  
double dilog(double x);
double trilog(double x);
