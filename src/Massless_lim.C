#include "Massless_lim.H"
#include "Massless_Kernels.H"

#include <iostream>
#include <vector>
#include <cmath>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
using namespace LHAPDF;
using namespace std;

Massless_Limit::Massless_Limit(const char *pdfname, const Scales scale) : 
  m2r(pow(scale.mu_R,2)),m2f(pow(scale.mu_F,2)),
  m2h(pow(scale.m_H,2)),m2b(pow(scale.mbmb,2)),S(pow(scale.ss,2)),i_nev(scale.i_evs)
{
  m_mode = scale.i_fonll;
  m_amin = scale.m_amin;
  p_lumi=new LUMI(pdfname,scale.i_member,scale.mu_F,scale.mbmb,false);
  m_cu = m_cd = 1.;
  m_proc = 1; // bbh
  // if(scale.m_proc_name == "bbh" || scale.m_proc_name != "bbz"){
  //   if(p_lumi->Init_Lumi()){
  //     m_inte=m_err=0;  
  //     lfr = log(m2f/m2r); 
  //     lfh = log(m2f/m2h);
  //     tauh = m2h/S;
  //     double mbr, cv, qmb, vev;
  //     qmb = 4.18;
			
  //     mbr   = p_lumi->mb(qmb,sqrt(m2r),4);
			
  //     cv   = 3.8937966e8;
  //     vev  = 246.221;
  //     cout << "bbH ==== \n Parameters: \n MB("<< sqrt(m2r)<< " ): " << mbr << "\n"
  // 	   << "tauh("<< sqrt(m2h) <<  "): " << tauh << "\n"
  // 	   << endl;
  //     s0   = M_PI/6. * cv/m2h * pow(mbr/vev,2);
  //     m_cu = m_cd = 1.;
			
  //   }
  // }
  // else if(scale.m_proc_name == "bbz"){
    if(p_lumi->Init_Lumi()){
      m_inte=m_err=0;  
      lfr = log(m2f/m2r); 
      lfh = log(m2f/m2h);
      tauh = m2h/S;
      double mbr, cv, qmb, vev;
      qmb = 4.18;
      vev  = 246.221;		
      mbr   = p_lumi->mb(qmb,sqrt(m2r),4);
			
      cv   = 3.8937966e8;
      vev  = 246.221;
      cout << "bbZ ==== \n Parameters: \n MB("<< sqrt(m2r)<< " ): " << mbr << "\n"
	   << "tauh("<< sqrt(m2h) <<  "): " << tauh/vev << "\n"
	   << endl;
      double GF   = 1.1663787e-5;
 
      double mw = 80.385, mz = sqrt(m2h);
      double sw = 1. - sqr(mw/mz);

      m_cu = 1 + sqr(-8./3.*sw + 1.);
      m_cd = 1 + sqr(4./3.*sw  - 1.);
      m_proc = 1;
      
      s0   = cv*M_PI*sqrt(2.)*GF/12.*m_cd;
      //      cout <<scientific << setprecision(15)<<sw << endl;
    }
    //  }

}

Massless_Limit::~Massless_Limit()
{
  if(p_lumi) delete p_lumi;
}


