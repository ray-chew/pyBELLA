/*******************************************************************************
 File:   limiter.h
 Author: Nicola         
 Date:   Mon Mar 16 10:23:53 CET 1998
 *******************************************************************************/
#ifndef LIMITER_H
#define LIMITER_H

#include "enumerator.h"


/*-----------------------------------------------------------------------------
 Limiter function
 double (*limiter[])(
 const double a, 
 const double b, 
 const double k);
 -----------------------------------------------------------------------------*/

double NoSlope(const double a, const double b, const double k);
double None(const double a, const double b, const double k);
double MinMod(const double a, const double b, const double k);
double MonotonizedCentral(const double a, const double b, const double k);
double Superbee(const double a, const double b, const double k);
double VanLeer(const double a, const double b, const double k);
double VanLeerS(const double a, const double b, const double k);
double SwebyMunz(const double a, const double b, const double k);
double Rupe(const double a, const double b, const double k);
double Rupe_velo(const double a, const double b, const double k);
double RupeZ(const double a, const double b, const double k);

#define NONE_LIM(a,b,k)      (0.5 * (a + b))
#define MINMOD_LIM(a,b,k)    ((a * b > 0.0) ? ((fabs(a) < fabs(b)) ? a : b) : 0.0)
#define MONCENT_LIM(a,b,k)   (a*b <= 0 ? 0.0 : (a*b < a*a ? (a*2.0*b < a*0.5*(a+b) ? 2.0*b : 0.5*(a+b)) : (a*2.0*a < a*0.5*(a+b) ? 2.0*a : 0.5*(a+b)))) 
#define SUPERBEE_LIM(a,b,k)  MAX_own(MINMOD_LIM(a, 2.0 * b, 42.), MINMOD_LIM(2.0 * a, b, 42.))
#define VANLEER_LIM(a,b,k)   ((2.0 * a * b > 0.0) ? 2.0 * a * b / (a + b) : 0.0)
#define SWEBYMUNZ_LIM(a,b,k) ( MAX_own(0,SIGN(a * b)) * SIGN(a) * MAX_own(MIN_own(k * fabs(a), fabs(b)), MIN_own(fabs(a), k * fabs(b))) )

#endif /* LIMITER_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
