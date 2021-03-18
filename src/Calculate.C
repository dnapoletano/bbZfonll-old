#include "Massless_lim.H"
#include "Massless_Kernels.H"
#include "PDF.H"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;

Massless_Limit* mlim_instance;

double wrap_LO(double x[], size_t dim, void* p){
  return mlim_instance->dsigmaLO(x, dim, p);
}

double wrap_NLO(double x[], size_t dim, void* p){
  if(dim == 1)  return  mlim_instance->dsigmaNLO(x, dim, p, 0);
  if(dim == 2)
    return
        mlim_instance->dsigmaNLO(x, dim, p, 0)
        - mlim_instance->dsigmaNLO(x, dim, p, 1);
  else return 0.;
}

double wrap_NNLO(double x[], size_t dim, void* p){
  if(dim == 1)   return  mlim_instance->dsigmaNNLO(x, dim, p, 0);
  if(dim == 2)
    return
        mlim_instance->dsigmaNNLO(x, dim, p, 0)
        - mlim_instance->dsigmaNNLO(x, dim, p, 1);
  else return 0.;
}

double xsec5F_LO(const Scales scal)
{
  double inte(0.), w2(0.);
  size_t l=1;
  gsl_monte_function F;
  double xmin[l], xmax[l];
  for(size_t i(0); i<l; ++i){
      xmin[i] = 0.; xmax[i] =1.;
    }

  F.f = &(wrap_LO);
  F.dim = l;
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(l);
  gsl_monte_vegas_init(state);
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&F,xmin,xmax,l,scal.i_evs,rng,state,&inte,&w2);
  gsl_monte_vegas_integrate(&F,xmin,xmax,l,scal.i_evs,rng,state,&inte,&w2);
  gsl_rng_free(rng);
  gsl_monte_vegas_free(state);
  // cout <<" Delta-LO ... : " << inte << " ± " << w2 << std::endl;
  return inte;
}

double xsec5F_NLO(const Scales scal)
{
  double inte_1(0.), w2_1(0.);
  double inte_2(0.), w2_2(0.);
  size_t l1(1), l2(2);
  gsl_monte_function Fd, Fp;
  double xmin1[l1], xmax1[l1];
  double xmin2[l2], xmax2[l2];
  xmin1[0]=0.; xmax1[0]=1.;
  for(size_t i(0); i<l2; ++i){
      xmin2[i] = 0.; xmax2[i] =1.;
    }

  Fd.dim = l1;  // delta-terms
  Fp.dim = l2;  // other terms
  Fd.f = &(wrap_NLO);
  Fp.f = &(wrap_NLO);

  gsl_monte_vegas_state* state1 = gsl_monte_vegas_alloc(l1);
  gsl_monte_vegas_init(state1);
  gsl_rng * rng1 = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&Fd,xmin1,xmax1,l1,
                            scal.i_evs,rng1,state1,
                            &inte_1,&w2_1);
  gsl_monte_vegas_integrate(&Fd,xmin1,xmax1,l1,
                            scal.i_evs,rng1,state1,
                            &inte_1,&w2_1);

  gsl_rng_free(rng1);
  gsl_monte_vegas_free(state1);
  // cout <<" Delta-NLO ... : " << inte_1 << " ± " << w2_1 << std::endl;


  gsl_monte_vegas_state* state2 = gsl_monte_vegas_alloc(l2);
  gsl_monte_vegas_init(state2);
  gsl_rng * rng2 = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&Fp,xmin2,xmax2,l2,
                            scal.i_evs,rng2,state2,
                            &inte_2,&w2_2);
  gsl_monte_vegas_integrate(&Fp,xmin2,xmax2,l2,
                            scal.i_evs,rng2,state2,
                            &inte_2,&w2_2);
  gsl_rng_free(rng2);
  gsl_monte_vegas_free(state2);
  // cout <<" Rest ........ : " << inte_2 << " ± " << w2_2 << std::endl;
  // cout << " ----------------- " << endl;
  // cout << inte_1+inte_2 << " ± " << w2_1 + w2_2 << endl;

  return inte_1+inte_2;
}

double xsec5F_NNLO(const Scales scal)
{
  double inte_1(0.), w2_1(0.);
  double inte_2(0.), w2_2(0.);
  size_t l1(1), l2(2);
  gsl_monte_function Fd, Fp;
  double xmin1[l1], xmax1[l1];
  double xmin2[l2], xmax2[l2];
  xmin1[0]=0.; xmax1[0]=1.;
  for(size_t i(0); i<l2; ++i){
      xmin2[i] = 0.; xmax2[i] =1.;
    }

  Fd.dim = l1;  // delta-terms
  Fp.dim = l2;  // other terms
  Fd.f = &(wrap_NNLO);
  Fp.f = &(wrap_NNLO);

  gsl_monte_vegas_state* state1 = gsl_monte_vegas_alloc(l1);
  gsl_monte_vegas_init(state1);
  gsl_rng * rng1 = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&Fd,xmin1,xmax1,l1,
                            scal.i_evs,rng1,state1,
                            &inte_1,&w2_1);
  gsl_monte_vegas_integrate(&Fd,xmin1,xmax1,l1,
                            scal.i_evs,rng1,state1,
                            &inte_1,&w2_1);

  gsl_rng_free(rng1);
  gsl_monte_vegas_free(state1);
  // cout <<" Delta-NNLO ... : " << inte_1 << " ± " << w2_1 << std::endl;


  gsl_monte_vegas_state* state2 = gsl_monte_vegas_alloc(l2);
  gsl_monte_vegas_init(state2);
  gsl_rng * rng2 = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&Fp,xmin2,xmax2,l2,
                            scal.i_evs,rng2,state2,
                            &inte_2,&w2_2);
  gsl_monte_vegas_integrate(&Fp,xmin2,xmax2,l2,
                            scal.i_evs,rng2,state2,
                            &inte_2,&w2_2);

  gsl_rng_free(rng2);
  gsl_monte_vegas_free(state2);
  // cout <<" Rest ........ : " << inte_2 << " ± " << w2_2 << std::endl;
  // cout << " ----------------- " << endl;

  return inte_1 + inte_2;
}


double xsec5F(const char *argv, const Scales scal)
{
  mlim_instance = new Massless_Limit(argv, scal);
  double asr = mlim_instance->masr(pow(scal.mu_R,2));


  double xs(0.);
  if ( scal.i_or == 0 ){
    // cout << " -----------  Computing LO 5F XS ---------  " << endl;
      xs  = xsec5F_LO(scal);
    }
  if ( scal.i_or == 1 ){
    // cout << " -----------  Computing NLO 5F XS ---------  " << endl;
      xs  = xsec5F_LO(scal) + asr*xsec5F_NLO(scal);
      // the nlo-res is slowly convergent
      //  cout << "NLO xs = " << xs << endl;
    }
  if ( scal.i_or == 2 ){
    // cout << " -----------  Computing NNLO 5F XS ---------  " << endl;
    double LO = xsec5F_LO(scal);
    double NLO =  xsec5F_NLO(scal);
    double NNLO = xsec5F_NNLO(scal);
    cout <<  LO << " , " << NLO << " , " << NNLO  << endl;
    xs = LO+asr*NLO+asr*asr*NNLO;
    }
  if (scal.i_or > 2 || scal.i_or < 0 ) xs = 0;
  delete mlim_instance;
  return xs;
}


double wrap_0LO(double x[], size_t dim, void* p){
  if(dim == 2) return mlim_instance->dsigma0LO(x, dim, p);
  else return 0.;
}

double wrap_0NLO(double x[], size_t dim, void* p){
  if(dim == 4) return mlim_instance->dsigma0NLO(x, dim, p);
  else return 0.;
}

double xsec4F0_LO(const Scales scal)
{
  double inte(0.), w2(0.);
  size_t l=2;
  gsl_monte_function F;
  double xmin[l], xmax[l];
  for(size_t i(0); i<l; ++i){
      xmin[i] = 0.; xmax[i] =1.;
    }

  F.f = &(wrap_0LO);
  F.dim = l;
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(l);
  gsl_monte_vegas_init(state);
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&F,xmin,xmax,l,scal.i_evs,rng,state,&inte,&w2);
  gsl_monte_vegas_integrate(&F,xmin,xmax,l,scal.i_evs,rng,state,&inte,&w2);
  gsl_rng_free(rng);
  gsl_monte_vegas_free(state);

  return inte;
}

double xsec4F0_NLO(const Scales scal)
{
  double inte(0.), w2(0.);

  size_t l(4);
  gsl_monte_function F;
  double xmin[l], xmax[l];
 
  int Evs= scal.i_evs;

  for(size_t i(0); i<l; ++i){
      xmin[i] = 0.; xmax[i] =1.;
    }

  F.f = &(wrap_0NLO);
  F.dim=l;

  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(l);
  gsl_monte_vegas_init(state);
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_params vegas_par;
  gsl_monte_vegas_params_get(state,&vegas_par);
  vegas_par.iterations = 20;
  gsl_monte_vegas_params_set(state,&vegas_par);
  gsl_monte_vegas_integrate(&F,xmin,xmax,l,Evs,rng,state,&inte,&w2);
  
  size_t rounds(0);
  cout << "iteration : " << rounds << "\n"
	 << " % err : " << std::fabs(w2/inte)*100
	 <<  " , chi2 : " << gsl_monte_vegas_chisq(state)  << endl;
  
  while (fabs (gsl_monte_vegas_chisq(state) - 1.0) > 0.5){
    vegas_par.iterations += 10;
    // if(rounds <3) vegas_par.stage = rounds+1;
    gsl_monte_vegas_params_set(state,&vegas_par);
    gsl_monte_vegas_integrate(&F,xmin,xmax,l,
			      size_t(Evs*(1+rounds)),rng,state,
			      &inte,&w2);
    rounds++;
    cout << "iteration : " << rounds << "\n"
	 << "int : " << inte << " % err : " << std::fabs(w2/inte)*100
	 <<  " , chi2 : " << gsl_monte_vegas_chisq(state)  << endl;
  }
  cout << " success ! vegas converged after " << rounds << " iterations !" << endl;
  cout << "int : " << inte << " % err : " << std::fabs(w2/inte)*100
       <<  " , chi2 : " << gsl_monte_vegas_chisq(state)  << endl;
  gsl_rng_free(rng);
  gsl_monte_vegas_free(state);
  // cout << " XS 4F O(as^3) .... " << inte << " ± " <<w2 << std::endl;
  // cout << " ----------------- " << endl;

  return inte;

}

double wrap_muFTerms(double x[], size_t dim, void* p)
{
  if(dim == 3) return mlim_instance->muFct(x, dim, p);
  else return 0.;
}

double muFTerm(const Scales scal)
{
  double inte(0.), w2(0.);

  size_t l(3);
  gsl_monte_function F;
  double xmin[l], xmax[l];
  int Evs= scal.i_evs;

  for(size_t i(0); i<l; ++i){
      xmin[i] = 0.; xmax[i] =1.;
    }

  F.f = &(wrap_muFTerms);
  F.dim=l;
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(l);
  gsl_monte_vegas_init(state);
  gsl_monte_vegas_params vegas_par;
  gsl_monte_vegas_params_get(state,&vegas_par);
  vegas_par.iterations = 20;
  gsl_monte_vegas_params_set(state,&vegas_par);

  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&F,xmin,xmax,l,
                            Evs,rng,state,
                            &inte,&w2);
  size_t rounds(0);
  while (fabs (gsl_monte_vegas_chisq(state) - 1.0) > 0.5){
    vegas_par.iterations += 10;
    // if(rounds <3) vegas_par.stage = rounds+1;
    gsl_monte_vegas_params_set(state,&vegas_par);
    //apply this mod to the xs4F0 and include muF in that!
    gsl_monte_vegas_integrate(&F,xmin,xmax,l,
			      size_t(Evs*(1+rounds)),rng,state,
			      &inte,&w2);
    rounds++;
    cout << "iteration : " << rounds << "\n"
	 << "int : " << inte << " % err : " << std::fabs(w2/inte)*100
	 <<  " , chi2 : " << gsl_monte_vegas_chisq(state)  << endl;
  }
  cout << " success ! vegas converged after " << rounds << " iterations !" << endl;
  cout << "int : " << inte << " % err : " << std::fabs(w2/inte)*100
       <<  " , chi2 : " << gsl_monte_vegas_chisq(state)  << endl;
  gsl_rng_free(rng);
  gsl_monte_vegas_free(state);

  return inte;
}


double xsec4F0(const char *argv, const Scales scal)
{
  double xs(0.), xs40(0.);
  mlim_instance = new Massless_Limit(argv, scal);
  if( scal.i_fonll == 1) {
    // cout << " -----------  Computing 4F0 XS O(as^2) ---------  " << endl;
      xs40  = xsec4F0_LO(scal);
      xs = xs40;
      // cout << "Massless-lim O(as^2) xs = .... " << xs << endl;
    }
  if( scal.i_fonll == 2) {
    // cout << " -----------  Computing 4F0 XS O(as^3) ---------  " << endl;

      double m2r = pow(scal.mu_R,2);
      double m2f = pow(scal.mu_F,2);
      double m2h = pow(scal.m_H,2);
      double b0 = (33.-10.)/12./M_PI;
      double asr = mlim_instance->masr(m2r);
      double ash = mlim_instance->masr(m2h);
      double asf = mlim_instance->masr(m2f);
      double A = xsec4F0_LO(scal);
      double B = xsec4F0_NLO(scal);
      double C = 0.;
   
      if(m2f!=m2h)
        C = muFTerm(scal);

      xs40 = A +(B + 2.*ash*b0*A*log(m2r/m2h)
                 - ash*C*log(m2f/m2h)/2./M_PI);

      xs = xs40;
      // cout  << "  LO = " << A/pow(asr,2) << " , NLO = " << xs << endl;
      // cout << " b0 : " << b0 << " , log(m2r/m2h) : " << log(m2r/m2h)
      // 	   << " , asr : " << asr << " ===== > " << 2.*asr*b0*A*log(m2r/m2h)
      // 	   << "\n" << " asf : " << asf << " C : " << C << " ,log(m2f/m2h) : "
      // 	   << log(m2f/m2h) << " ====== > " << -asf*C*log(m2f/m2h)/2./M_PI <<
      // 	endl;
      // cout << " -------------------------------------------\n" ;
      // cout << " ------ Massless lim O(as^3) xs = ... " << xs<< " ------"<<endl;
      // cout << " -------------------------------------------\n" ;
    }

  if (scal.i_fonll > 2 || scal.i_fonll < 0) xs40=0;
  
  xs = xs40;
  delete mlim_instance;
  return xs;
}


double xsecDiff(const char *argv, const Scales scal)
{
  double xs(0.), xs5(0.), xs40(0.);
  xs5 = xsec5F(argv,scal); xs40 = xsec4F0(argv,scal);
  // cout << " ############################################### " << endl
  //      << " XS 5F : \t " <<  xs5 << endl
  //      << " XS 4F : \t " <<  xs40 << endl
  //      << " ############################################### " << endl;
  xs = xs5 - xs40;

  return xs;
}

// int main()
// {
//   double mu0 = 91.1876,mu1=(91.1876+ 2.*4.92)/3.;
//   //order 5f, order fonll, mb_pole, muR, muF, mV, eCM, evs, member, amin
//   Scales sc(-1,2,4.92, 91.1876, 91.1876, 91.1876, 13000., 10000,0,1.e-12);
//   // cout << xsec5F("NNPDF31_nnlo_as_0118",sc) << endl;
//   // sc.mu_F=0.167*mu1; sc.mu_R=mu1;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   // sc.mu_F=mu1; sc.mu_R=0.167*mu1;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   // sc.mu_F=mu1; sc.mu_R=mu1;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   // sc.mu_F=mu1; sc.mu_R=2.*mu1;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   // sc.mu_F=2.*mu1; sc.mu_R=mu1;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   // sc.mu_F=0.167*mu0; sc.mu_R=mu0;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   sc.mu_F=mu0; sc.mu_R=mu0;
//   cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   sc.mu_F=0.167*mu0; sc.mu_R=mu0;
//   cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   // sc.mu_F=mu0; sc.mu_R=2.*mu0;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   // sc.mu_F=2.*mu0; sc.mu_R=mu0;
//   // cout << xsec4F0("NNPDF31_nnlo_as_0118",sc) << endl;

//   return 0;
// }
