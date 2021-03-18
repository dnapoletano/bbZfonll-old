#ifndef MASSLESS_KERNELS_H
#define MASSLESS_KERNELS_H

extern "C" {
  // delta coefficients
  double delta1_(double *lfr);  // delta NLO
  double delta2_(double *lfr, double *lfh); // delta NNLO

  // plus coefficients
  double dterms1_(double *xt, double *lfh); // plus NLO
  double dterms2_(double *xt, double *lfh, double *lfr); //plus NNLO
  
  // plus subtraction terms
  double dtsub1_(double *lfh, double *tauh); //sub NLO
  double dtsub2_(double *lfh, double *lfr, double *tauh); // sub NNLO

  // NLO terms
  double sbbq1_(double *xt, double *lfh); // b b~ NLO minus pure soft
  double sbg1_(double *xt, double *lfh); // b g NLO exact

  // NNLO terms
  // b b~
  double deltabbqa_(double *xx, double *lfh, double *lfr);
  double deltabbqf_(double *xx, double *lfh, double *lfr);

  // b g
  double deltabga_(double *xx, double *lfh, double *lfr);
  double deltabgf_(double *xx, double *lfh, double *lfr);

  // g g
  double deltagga_(double *xx, double *lfh, double *lfr);
  double deltaggf_(double *xx, double *lfh, double *lfr);

  // b b
  double deltabba_(double *xx, double *lfh, double *lfr);
  double deltabbf_(double *xx, double *lfh, double *lfr);

  // b q
  double deltabqa_(double *xx, double *lfh, double *lfr);
  double deltabqf_(double *xx, double *lfh, double *lfr);

  // q q~
  double deltaqqba_(double *xx, double *lfh, double *lfr);
  double deltaqqbf_(double *xx, double *lfh, double *lfr);

  // generalized polylogs
  double wgplg_(int *N, int*P, double *x);
  double trilog_(double *xx);

  // runmass @ 4loops
  double runmass_(const double *mass0, const double *api0,
		  const double *apif, const double *nf, const int *nloop);
};

#endif
