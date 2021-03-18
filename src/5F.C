////
// Calculate the massless scheme at LO, NLO and NNLO
///

#include "Massless_lim.H"
#include "Massless_Kernels.H"

#include <iostream>
#include <vector>
#include <cmath>
using namespace LHAPDF;
using namespace std;

double Massless_Limit::dsigmaLO(double x[], size_t dim, void* p)
{
  double dsg(0.);
  double tau, jac(1.);
  tau  = tauh + (1. - tauh)*x[0];
  jac *= (1. - tauh);
  dsg = 2.*jac*p_lumi->Lumi(10,0,tau,tauh/tau,1.,1.,10)/tau;
  return s0*dsg;
}

double Massless_Limit::dsigmaNLO(double z[], size_t dim, 
				 void* p, const size_t plus) 
// here plus means should I do the subtraction? yes=1,no=0
{
  double dsg(0.);
  // delta contr
  if(dim==1){
    double tau, jac(1.);
    double sbb(0.);
    tau  = tauh + (1. - tauh)*z[0];
    jac *= (1. - tauh);
    if(m_proc==0)
      sbb    = delta1_(&lfr);
    else if(m_proc==1) // bb->Z
      sbb    = delta1_(&lfh);;
    dsg = 2.*
      sbb*jac*p_lumi->Lumi(10,0,tau,tauh/tau,1.,1.,10)/tau;
    return s0*dsg;
  }
  // main thing
  if(dim==2 && plus == 0){
    double tau, jac(1.), xi, x, y;
    double sbbp(0.), sbb(0.), sbg(0.);
    tau  = tauh + (1. - tauh)*z[0];
    jac *= (1. - tauh);
    double ximax = - 0.5*log(tau);
    xi   = -ximax + 2.*ximax*z[1];
    jac *=  2.*ximax;

    // standard variables  x, y ;
    x = sqrt(tau)*exp(xi);
    y = sqrt(tau)*exp(-xi);
    sbbp    = dterms1_(&y, &lfh);
    sbb     = sbbq1_(&y, &lfh);
    sbg     = sbg1_(&y,&lfh); /////////
    dsg     = 
      jac*(
	      2.*p_lumi->Lumi(10,0,x,tauh/tau,1.,1.,10)*
	      (sbbp*y + sbb)
	      - 2.*p_lumi->Lumi(10,0,x,tauh/x,1.,1.,10)*(y*sbbp)
	      + (p_lumi->Lumi(10,5,x,tauh/tau,1.,1.,10) + 
		 p_lumi->Lumi(5,10,x,tauh/tau,1.,1.,10) + 
		 p_lumi->Lumi(0,5,x,tauh/tau,1.,1.,10)  + 
		 p_lumi->Lumi(5,0,x,tauh/tau,1.,1.,10))*(sbg)
	      )/tau;
    return s0*dsg;
  }

  if(dim==2 && plus == 1){
    double tau, jac(1.), xi, x, y;
    double sbbp(0.);
    tau  = tauh + (1. - tauh)*z[0];
    jac *= (1. - tauh);
    xi   = (tauh/tau)*z[1];
    jac *= (tauh/tau);
    x = tau; y = xi;
    sbbp    = dterms1_(&y, &lfh);

    dsg = 2.*jac*sbbp*p_lumi->Lumi(10,0,x,tauh/x,1.,1.,10)/x;
    return s0*dsg;
  }
  else if (dim > 2){
    cout << " No more than 2 dim for dsigmaNLO! " << std::endl;
    return 0.;
  }
  return s0*dsg;
}

double Massless_Limit::dsigmaNNLO(double z[], size_t dim, 
				  void* p, const size_t plus) 
// here plus means should I do the subtraction? yes=1,no=0
{
  //cout << " MODE:: "<< m_mode << std::endl;
  double dsg(0.), nf(5.);
  // delta contr
  if(dim==1){
    double tau, jac(1.);
    double sbba(0.);
    tau  = tauh + (1. - tauh)*z[0];
    jac *= (1. - tauh);
    sbba    = delta2_(&lfr, &lfh);
    dsg = 2.*sbba *
      jac*p_lumi->Lumi(10,0,tau,tauh/tau,1.,1.,10)/tau;
    return s0*dsg;
  }
  // main thing
  if(dim==2){ 
    if (plus == 0) {
      double tau, jac(1.), xi, x, y;
      double lbb, lbb0, lbg, lbq, lgg,lqq, lbs;
      double sbbp(0.), sbba(0.), sbbf(0.),
	sbg(0.), sbq(0.), sgg(0.), sbs(0.), sqq(0.);
      
      tau  = tauh + (1. - tauh)*z[0];
      jac *= (1. - tauh);
      double ximax = - 0.5*log(tau);
      xi   = -ximax + 2.*ximax*z[1];
      jac *=  2.*ximax;

      // standard variables  x, y ;
      x = sqrt(tau)*exp(xi);
      y = sqrt(tau)*exp(-xi);

      //soft-terms

      sbbp      =       dterms2_(&y,&lfh,&lfr);
      //      finite-terms
      //b b~
      sbba       =       deltabbqa_(&y,&lfh,&lfr);
      sbbf       =       deltabbqf_(&y,&lfh,&lfr);
      
      // // b g
      sbg        =       deltabga_(&y,&lfh,&lfr) +
       	nf* deltabgf_(&y,&lfh,&lfr);
      
      // // g g
      sgg        =       deltagga_(&y,&lfh,&lfr);
      
      // // b b
      sbs        =       deltabba_(&y,&lfh,&lfr);

      // // b q      
      sbq        =       deltabqa_(&y,&lfh,&lfr);

      // // q q~
      sqq        =       deltaqqba_(&y,&lfh,&lfr);

      // Luminosities 

      lbb     = 2.*p_lumi->Lumi(10,0,x,tauh/tau,1.,1.,10);
    
      lbb0    = 2.*p_lumi->Lumi(10,0,x,tauh/x,1.,1.,10);

      lbg     = p_lumi->Lumi(10,5,x,tauh/tau,1.,1.,10) + 
	p_lumi->Lumi(5,10,x,tauh/tau,1.,1.,10) +
	p_lumi->Lumi(0,5,x,tauh/tau,1.,1.,10) + 
	p_lumi->Lumi(5,0,x,tauh/tau,1.,1.,10);

      lbq     = 0.;
   
      lbq += 2.*( p_lumi->Lumi(10,4,x,tauh/tau,1.,1.,10)  // d
		  + p_lumi->Lumi(10,6,x,tauh/tau,1.,1.,10));
      lbq += 2.*( p_lumi->Lumi(4,10,x,tauh/tau,1.,1.,10)			 
		  + p_lumi->Lumi(6,10,x,tauh/tau,1.,1.,10));

      lbq += 2.*( p_lumi->Lumi(10,3,x,tauh/tau,1.,1.,10) // u
		  + p_lumi->Lumi(10,7,x,tauh/tau,1.,1.,10));
      lbq += 2.*( p_lumi->Lumi(3,10,x,tauh/tau,1.,1.,10) 
		  + p_lumi->Lumi(7,10,x,tauh/tau,1.,1.,10));

      lbq += 2.*( p_lumi->Lumi(10,2,x,tauh/tau,1.,1.,10) // s
		  + p_lumi->Lumi(10,8,x,tauh/tau,1.,1.,10));
      lbq += 2.*( p_lumi->Lumi(3,10,x,tauh/tau,1.,1.,10) 
		  + p_lumi->Lumi(8,10,x,tauh/tau,1.,1.,10));

      // (c = c~)
      lbq += 4.*( p_lumi->Lumi(10,9,x,tauh/tau,1.,1.,10) );
      lbq += 4.*( p_lumi->Lumi(9,10,x,tauh/tau,1.,1.,10) );
     

      lgg     = p_lumi->Lumi(5,5,x,tauh/tau,1.,1.,10);

      lqq     = L_qq(x,tauh/tau);
      
      lbs     = 2.*p_lumi->Lumi(10,10,x,tauh/tau,1.,1.,10);

      // cout << sbbp << " , " << sbba << " , " << sbbf << " , "
      // 	   << sbg << " , " << sgg << " , " << sbs << " , " 
      // 	   << sbq << " , " << sqq << endl;
      
      dsg     = 
	jac*( lbb * ( sbbp*y + sbba + nf * sbbf)
		  - lbb0 * sbbp *y
		  + lbg  * sbg
		  + lgg  * sgg
		  + lbs  * sbs
		  + lbq  * sbq
		  + lqq  * sqq
		  )/tau; 
      return s0*dsg;
    }
    
    if(plus == 1){
      double tau, jac(1.), xi, x, y;
      double sbbp(0.);
      tau  = tauh + (1. - tauh)*z[0];
      jac *= (1. - tauh);
      xi   = (tauh/tau)*z[1];
      jac *= (tauh/tau);
      x = tau; y = xi;
      sbbp    = dterms2_(&y,&lfh,&lfr);

      dsg = 2.*jac*sbbp*p_lumi->Lumi(10,0,x,tauh/x,1.,1.,10)/x;
      return s0*dsg;
    }
  }
  else if (dim > 2){
    cout << " No more than 2 dim for dsigmaNLO! " << std::endl;
    return 0.;
  } 
  return s0*dsg;
}


