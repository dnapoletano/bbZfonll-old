#include "Massless_lim.H"
using namespace std;

double CalcggLO(const int pos, const int ss);
double CalcLO(const int pos, const int ss);
double CalcNLO(const int pos, const int ss);
double Calc4FLO(const int pos);
double Calc4FNLO(const int pos);
double Calc4F_mb(const int pos);
double Calc4F0_mb(const double mb,const int order);
double kappa(const double x);
int findpos(const double muR, const double muRmin, 
	    const double muRmax, const int Npoints);
double Calc_4F_matched(const char* argv, 
		       const Scales scale, 
		       const int o_fonll, const int ss);
