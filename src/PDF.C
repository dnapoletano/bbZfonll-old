#include "LHAPDF/LHAPDF.h"
#include "PDF.H"
#include "Massless_Kernels.H"
#include <cmath>

#define ZETA2 1.6449340668482264365
#define ZETA3 1.2020569031595942854

using namespace LHAPDF;
using namespace std;

LUMI::LUMI(const char* setname, const int setmem,
	   const double muF, const double mbext,
	   const bool ismlim) :
  m_pdf(mkPDF(setname,setmem)),m_as(0.),Q(muF),m_flag(false), m_lim(ismlim)
{
  if(m_pdf) m_flag=true;
  m_as = m_pdf->alphasQ(Q)/(2.*M_PI);
  mm2b = mbext*mbext;
  m2b = m_pdf->alphaS().quarkMass(5)*m_pdf->alphaS().quarkMass(5);

}
LUMI::~LUMI()
{
  if(m_pdf) delete m_pdf;
}


double LUMI::mb(const double mf, const double mu, const int nloop)
{
  double as,asmf,nf(5.);
  double mfrun(0.);
  double mu0=mf;
  as     = m_pdf->alphasQ(mu)/M_PI;
  asmf   = m_pdf->alphasQ(mu0)/M_PI;
  mfrun = runmass_(&mf,&asmf,&as,&nf,&nloop);
 return mfrun;
}



double LUMI::Lumi(const int i, const int j,
		  const double x, const double y,
		  const double z1, const double z2,
		  const int order)
{
  double Q2 = Q*Q;
  // cout << " muF = " << Q << "\n";

  if(m_lim == true) m2b = mm2b;

  if(order==10){
      return m_pdf->xfxQ2(i-5,x,Q2)*m_pdf->xfxQ2(j-5,y,Q2);
    }
  if(order==1){
      if( (i == 10 || i == 0) && (j!=10 && j!=0)){
          // double btildex = m_as*log(Q2/m2b)*Pqg(z1)*m_pdf->xfxQ2(0,x/z1,Q2);
          double btildex = log(Q2/m2b)*Pqg(z1)*m_pdf->xfxQ2(0,x/z1,Q2);
          return btildex*m_pdf->xfxQ2(j-5,y,Q2);
        }
      if( (i!=10 && i!=0) && (j == 10 || j == 0)){
          // double btildey = m_as*log(Q2/m2b)*Pqg(z2)*m_pdf->xfxQ2(0,y/z2,Q2);
          double btildey = log(Q2/m2b)*Pqg(z2)*m_pdf->xfxQ2(0,y/z2,Q2);
          return m_pdf->xfxQ2(i-5,x,Q2)*btildey;
        }

      if( (i == 10 || i == 0) && (j == 10 || j == 0) ){
          // double btildex = m_as*log(Q2/m2b)*Pqg(z1)*m_pdf->xfxQ2(0,x/z1,Q2);
          // double btildey = m_as*log(Q2/m2b)*Pqg(z2)*m_pdf->xfxQ2(0,y/z2,Q2);
          double btildex = log(Q2/m2b)*Pqg(z1)*m_pdf->xfxQ2(0,x/z1,Q2);
          double btildey = log(Q2/m2b)*Pqg(z2)*m_pdf->xfxQ2(0,y/z2,Q2);
          return btildex*btildey;
        }
    }
  if(order==2){
      double L = log(Q2/m2b);
      if( (i == 10 || i == 0) && (j!=10 && j!=0 )){
          // double bt2x = m_as*m_as*(L*L*phi(x,z1,Q2,2)+L*phi(x,z1,Q2,1));
          double bt2x = L*L*phi(x,z1,Q2,2) + L*phi(x,z1,Q2,1);
          if(sqrt(Q2) != 125.){
              bt2x += phi(x,z1,Q2,0);
            }
          return bt2x*m_pdf->xfxQ2(j-5,y,Q2);
        }
      if( (i!=10 && i!=0) && (j == 10 || j == 0)){
          // double bt2y = m_as*m_as*(L*L*phi(y,z2,Q2,2)+L*phi(y,z2,Q2,1));
          double bt2y = L*L*phi(y,z2,Q2,2)+L*phi(y,z2,Q2,1);
          if(sqrt(Q2) != 125.){
              bt2y += phi(y,z2,Q2,0);
            }
          return m_pdf->xfxQ2(i-5,x,Q2)*bt2y;
        }
      if( (i == 10 || i == 0) && (j == 10 || j == 0) ){
          // double bt1x = m_as*L*Pqg(z1)*m_pdf->xfxQ2(0,x/z1,Q2);
          // double bt1y = m_as*L*Pqg(z2)*m_pdf->xfxQ2(0,y/z2,Q2);
          // double bt2x = m_as*m_as*(L*L*phi(x,z1,Q2,2)+L*phi(x,z1,Q2,1));
          // double bt2y = m_as*m_as*(L*L*phi(y,z2,Q2,2)+L*phi(y,z2,Q2,1));

          double bt1x = L*Pqg(z1)*m_pdf->xfxQ2(0,x/z1,Q2);
          double bt1y = L*Pqg(z2)*m_pdf->xfxQ2(0,y/z2,Q2);
          double bt2x = L*L*phi(x,z1,Q2,2) + L*phi(x,z1,Q2,1);
          double bt2y = L*L*phi(y,z2,Q2,2) + L*phi(y,z2,Q2,1);

          if(sqrt(Q2) != 125.){
              bt2x += phi(x,z1,Q2,0);
              bt2y += phi(y,z2,Q2,0);
            }

          return bt1x*bt2y + bt2x*bt1y;
        }
    }
  return 0.;
}

double polevl(double x, double* coef, int N )
{
  double ans;
  int i;
  double *p;
  
  p = coef;
  ans = *p++;
  i = N;
  
  do
    ans = ans * x  +  *p++;
  while( --i );
  
  return ans;
}

double dilog(double x)
{
  static double cof_A[8] = {
    4.65128586073990045278E-5,
    7.31589045238094711071E-3,
    1.33847639578309018650E-1,
    8.79691311754530315341E-1,
    2.71149851196553469920E0,
    4.25697156008121755724E0,
    3.29771340985225106936E0,
    1.00000000000000000126E0,
  };
  static double cof_B[8] = {
    6.90990488912553276999E-4,
    2.54043763932544379113E-2,
    2.82974860602568089943E-1,
    1.41172597751831069617E0,
    3.63800533345137075418E0,
    5.03278880143316990390E0,
    3.54771340985225096217E0,
    9.99999999999999998740E-1,
  };
  if( x >1. ) {
      return -dilog(1./x)+M_PI*M_PI/3.-0.5*sqr(log(x));
    }
  x = 1.-x;
  double w, y, z;
  int flag;
  if( x == 1.0 )
    return( 0.0 );
  if( x == 0.0 )
    return( M_PI*M_PI/6.0 );
  
  flag = 0;
  
  if( x > 2.0 )
    {
      x = 1.0/x;
      flag |= 2;
    }
  
  if( x > 1.5 )
    {
      w = (1.0/x) - 1.0;
      flag |= 2;
    }
  
  else if( x < 0.5 )
    {
      w = -x;
      flag |= 1;
    }
  
  else
    w = x - 1.0;
  
  
  y = -w * polevl( w, cof_A, 7) / polevl( w, cof_B, 7 );
  
  if( flag & 1 )
    y = (M_PI * M_PI)/6.0  - log(x) * log(1.0-x) - y;
  
  if( flag & 2 )
    {
      z = log(x);
      y = -0.5 * z * z  -  y;
    }

  return y;

}

double trilog(double x)
{
  double z = x;
  return trilog_(&z);
}


double LUMI::Pqg(const double x)
{
  return 0.5*(x*x+pow((1.0-x),2));
}

double LUMI::Pqq_pt(const double x)
{
  return 4./3.*(1.+x*x);
}
double LUMI::Pqq_plus(const double x)
{
  return 1./(1.-x);
}

double LUMI::Pqq_psub(const double x)
{
  return -log(1. - x);
}

double LUMI::Pqq_d()
{
  return 2.;
}

double LUMI::Pgq(const double x)
{
  return 4./3.*(1.+(1.-x)*(1.-x))/x;
}

double LUMI::Pgg_pt(const double x)
{
  return 6.*x;
}

double LUMI::Pgg_plus(const double x)
{
  return Pqq_plus(x);
}

double LUMI::Pgg_psub(const double x)
{
  return -log(1. - x);
}

double LUMI::Pgg_d()
{
  return (33.-10.)/6.;
}

double LUMI::Pgg(const double x)
{
  return 6.*((1.-x)/x + x*(1.-x));
}

double LUMI::phi(const double x, const double z, 
                 const double Q2, const int o)
{
  double sigma(0.);
  for(size_t i(0);i<4;++i){
      sigma+=m_pdf->xfxQ2(i+1,x/z,Q2)+m_pdf->xfxQ2(-(i+1),x/z,Q2);
    }
  return agb(z,2,o)*m_pdf->xfxQ2(0,x/z,Q2) + asb(z,2,o)*sigma;
}

double LUMI::agb(const double z, const int fo, const int lo)
{
  double CF(4./3.), CA(3.), TR(0.5), lz(log(z)),lIz(log(1.-z)),
      z2(z*z);
  double res(0.), cftr(0.), catr(0.), tr2(0.);
  if(fo==1 && lo==1) res=Pqg(z);
  if(fo==2 && lo==0) {
      double Iz(1. - z);
      int n(1),p(2);
      cftr = 2.*Pqg(z)*(8.*ZETA3 + 4./3.*pow(lIz,3) - 8.*lIz*dilog(Iz)
                        +8.*ZETA2*lz - 4.*lz*pow(lIz,2) + 2./3.*pow(lz,3)
                        - 8.*lz*dilog(1.-z) + 8.*trilog(1.-z) - 24.*wgplg_(&n,&p,&Iz))
          + z2*(-16.*ZETA2*lz + 4./3.*pow(lz,3) + 16.*lz*dilog(Iz)+32*wgplg_(&n,&p,&Iz))
          - (4. + 96.*z - 64.*z2)*dilog(Iz) - (4.-48.*z + 40.*z2)*ZETA2 - (8.+48.*z - 24.*z2)*lz*lIz
          + (4.+8.*z- 12.*z2)*pow(lIz,2) - (1.+12.*z-20.*z2)*pow(lz,2)
          - (52.*z - 48.*z2)*lIz - (16.+18.*z+48.*z2)*lz
          + 26. - 82.*z + 80.*z2;

      double mz = -z;
      catr = 2.*Pqg(z)*(-4./3.*pow(lIz,3) + 8.*lIz*dilog(Iz) - 8.*trilog(Iz))
          +(1. + 2.*z + 2.*z2)*(-8.*ZETA2*log(1.+z) - 16.*log(1.+z)*dilog(-z)
                                -8.*lz*pow(log(1.+z),2) + 4.*pow(lz,2)*log(1.+z) + 8.*lz*dilog(-z)
                                -8.*trilog(-z) - 16.*wgplg_(&n,&p,&mz))
          + (16. + 64.*z)*(2.*wgplg_(&n,&p,&Iz)+lz*dilog(Iz))
          -(4./3. + 8./3.*z)*pow(lz,3)
          + (8.-32.*z + 16.*z2)*ZETA3
          -(16.+64.*z)*ZETA2*lz
          + (16.*z + 16.*z2)*(dilog(-z) +
                              lz*log(1.+z))
          + (32./3./z + 12. + 64.*z - 272./3.*z2)*dilog(Iz)
          - (12. + 48.*z - 260./3.*z2 + 32./3./z)*ZETA2
          - 4.*z2*lz*lIz
          - (2.+8.*z- 10.*z2)*pow(lIz,2)
          + (2. + 8.*z + 46./3.*z2)*pow(lz,2)
          + (4.+16.*z - 16.*z2)*lIz
          - (56./3. + 172./3.*z + 1600./9.*z2)*lz
          - 448./27./z - 4./3. - 628./3.*z + 6352./27.*z2;

      res = (CF*TR*cftr + CA*TR*catr)/8.;  // 8 is = 4 wrt Maria & 2 wrt to Buza...

    }
  if(fo==2 && lo==1) {
      cftr = 8.*Pqg(z)*(2.*lz*lIz - lIz*lIz + 2.*ZETA2)
          -(2. - 4.*z + 8.*z2)*lz*lz - 16.*z*(1.-z)*lIz
          -(6. - 8.*z + 16.*z2)*lz - 28. + 58.*z - 40.*z2;

      catr = (8. + 16.*z + 16.*z2)*(dilog(-z)+log(z)*log(1+z))
          + 8.*Pqg(z)*lIz*lIz
          + (4. + 8.*z)*lz*lz + 16.*z*ZETA2 + 16.*z*(1.-z)*lIz
          - (4. + 32.*z + 176.*z2/3.)*lz
          - 80./(9.*z) + 8. - 100.*z + 872.*z2/9.;

      res = - (CF*TR*cftr + CA*TR*catr)/4.;
    }
  if(fo==2 && lo==2){
      cftr = 8.*Pqg(z)*lIz - (2. - 4.*z + 8.*z2)*lz - (1. - 4.*z);

      catr = 8.*Pqg(z)*lIz + (4. + 16.*z)*lz
          + 8./(3.*z) + 2. + 16.*z - 62.*z2/3.;

      tr2 = -8.*(z2+(1.-z)*(1.-z))/3.;
      double tr(0.);
      tr = 8.*Pqg(z)/3.; // 8./3.*Pqg(z);

      res = (CF*TR*cftr - CA*TR*catr + TR*TR*tr2 + TR*tr)/4.;
    }
  return res;
}

double LUMI::asb(const double z, const int fo, const int lo)
{
  double CF(4./3.), TR(0.5), lz(log(z)), z2(z*z);
  double res(0.), cftr(0.);
  if(fo==2 && lo==0) {
      double a1(0.), a2(0.), a3(0.), a4(0.), a5(0.),
          a6(0.);
      double Iz(1. - z);
      int n(1),p(2);

      a1 = 32.*wgplg_(&n,&p,&Iz) + 16.*lz*dilog(Iz) - 16.*ZETA2*lz
          - 4./3.*pow(lz,3);
      a2 = 32./3./z + 8. -8.*z - 32./3.*z2;
      a3 = -a2;
      a4 = 2. + 10.*z+ 16./3.*z2;
      a5 = -(56./3.+88./3.*z + 448./9.*z2);
      a6 = -448./27./z - 4./3. - 124./3.*z + 1600./27.*z2;

      res =CF*TR*((1.+z)*a1 + dilog(Iz)*a2 + ZETA2*a3 + pow(lz,2)*a4
                  + lz*a5 + a6);
      res /= 8.;
//      res = 0.;

    }
  if(fo==2 && lo==2) {
      cftr = - 4.*(1. + z)*lz - 8./(3.*z) - 2. + 2.*z + 8.*z2/3. ;
      res = (CF*TR*cftr)/4.;
    }
  if(fo==2 && lo==1){
      cftr =- 4.*(1. + z)*lz*lz + (4. + 20.*z + 32.*z2/3.)*lz
          + 80./(9.*z) - 8. + 24.*z - 224.*z2/9.;

      res = (CF*TR*cftr)/4.;
    }
  return res;
}
