/*******************************************************************************
 File:   limiter.c
 Author: Nicola         
 Date:   Mon Mar 16 10:23:53 CET 1998
 *******************************************************************************/
#include <float.h>
#include <math.h>
#include "Common.h"
#include "math_own.h"
#include "enumerator.h"
#include "limiter.h" 


double smoothabs(double x, double xabs, double Xabs);


double NoSlope(const double a, const double b, const double k) {
	return 0.0;
}

double None(const double a, const double b, const double k) {
	return 0.5 * (a + b);
}

double MinMod(const double a, const double b, const double k) {
	return (a * b > 0.0) ? ((fabs(a) < fabs(b)) ? a : b) : 0.0;
}


double VanLeer(const double a, const double b, const double k) {
	const double c = 2.0 * a * b;
	return (c > 0.0) ? c / (a + b) : 0.0;
}

double VanLeerS(const double a, const double b, const double k) {
	const double absa  = fabs(a);
	const double absb  = fabs(b);
    const double abmin = MIN_own(absa,absb);
    const double abmax = MAX_own(absa,absb);
    const double x     = (absa < absb ? a : b);
    const double xabs  = abmin;
    const double X     = (absa < absb ? b : a);
    const double Xabs  = abmax;
    const double xabss = smoothabs(x, xabs, Xabs);
	return ((x*Xabs + X*xabss) / (Xabs + xabss + DBL_EPSILON));
}

double smoothabs(double x, double xabs, double Xabs) {
    double abseps = 0.1;
    double Xeps = abseps*Xabs;
    double sabs0  = sqrt(x*x + Xeps*Xeps);
    return (sabs0 + (xabs/(Xabs+DBL_EPSILON)) * (xabs - sabs0));
}


double Superbee(const double a, const double b, const double k) {
	return MAX_own(MinMod(a, 2.0 * b, 42.), MinMod(2.0 * a, b, 42.));
}


double MonotonizedCentral(const double a, const double b, const double k) {
	if(a * b > 0.0) {
		const double s = SIGN(a);
		const double d = 0.5 * (a + b);
		double c = a;
		if(s * b < s * a) c = b;
		c *= 2.0;
		if(s * c < s * d) return c;
		return d;
	}
	return 0.0;
}


double SwebyMunz(const double a, const double b, const double k) {
	const double sab = SIGN(a * b);
	const double sa = SIGN(a);
	const double as = fabs(a);
	const double bs = fabs(b);
	return MAX_own(0,sab) * sa * MAX_own(MIN_own(k * as, bs), MIN_own(as, k * bs));
}

#define PSI(r) (1.0 + r * (1-r) * (1-r*r)) 
/* #define PSI(r) (1.0) */
/* #define PSI(r) (1.0 + r * (1-r) * (1-r)) */
/* #define PSI(r) (1.0 + r * (1-r) * (1-r*r)) */
/* #define PSI(r) (1.0 + r * (1-r) * (1-r*r*r)) */
/* #define PSI(r) (1.0 + r * (1-r) * (1-r*r*r*r))  */

#define PSIZ(r) (1.0 + r * (1-r) * (1-r*r))
/* #define PSIZ(r) (1.0 - r * (1-r) * (1-r)) */
/* #define PSIZ(r) (1.0 - r * (1-r) * (1-r*r))  */
/* #define PSIZ(r) (1.0 - r * (1-r) * (1-r*r*r*r)) */
/* #define PSIZ(r) (1.0 - 0.5 * (1-r) * pow((1-r),0.25)) */

#ifdef EXTREMA_TRICK
double Rupe(const double a, const double b, const double k) {
	const double fa = fabs(a);
	const double fb = fabs(b);
	const double c = 2.0 * fa * fb;
	
	return ( c > DBL_EPSILON 
			? (fa/fb < 1 ? PSI((fa/fb)) : PSI((fb/fa))) * c / (fa + fb)
	       : 0.0);
}
#else
double Rupe(const double a, const double b, const double k) {
	const double c = 2.0 * a * b;
	return ((c > 0.0) 
            ? (a/b < 1 ? PSI((a/b)) : PSI((b/a))) * c / (a + b) 
			: 0.0);
}
#endif

#ifdef EXTREMA_TRICK
double RupeZ(const double a, const double b, const double k) {
	const double fa = fabs(a);
	const double fb = fabs(b);
	const double c = 2.0 * fa * fb;
	
	return ( c > DBL_EPSILON 
			? (fa/fb < 1 ? PSI((fa/fb)) : PSI((fb/fa))) * c / (fa + fb)
			: 0.0);
}
#else
double RupeZ(const double a, const double b, const double k) {
	const double c = 2.0 * a * b;
	return ((c > 0.0) 
            ? (a/b < 1 ? PSIZ((a/b)) : PSIZ((b/a))) * c / (a + b) 
			: 0.0);
}
#endif

/*
double RupeZ(const double a, const double b, const double k) {
	const double c = 2.0 * a * b;
	return ((c > 0.0) 
            ? (a/b < 1 ? PSI((a/b)) : PSI((b/a))) * c / (a + b) 
			: 0.0);
}



 A smooth non-standard limiter with MINMOD-type limiting near extrema:
 
 double Rupe(const double a, const double b, const double k) {
 const double c = 2.0 * a * b;
 double fa, fb;
 return ((c > 0.0) 
 ? c * (a/b < 1 ? PSI((a/b)) : PSI((b/a))) / fabs(a + b) 
 : ((fa = fabs(a)) < (fb = fabs(b))) ? fa : fb);
 }
 
 */


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
