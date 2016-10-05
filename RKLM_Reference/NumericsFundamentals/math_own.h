/*******************************************************************************
 File:   math_own.h
 Author: Nicola						Rupert
 Date:   Thu Feb 12 01:22:28 CET 1998  Feb. 2004
 *******************************************************************************/
#ifndef MATH_H
#define MATH_H

#include <math.h>
/* Contants ---------------------------------------------------------------- */
#define  PI          3.141592653589793100 
#define  ONE_THIRD   0.333333333333333333 
#define  ONE_SIXTH   0.166666666666666667 
#define  FIVE_SIXTHS 0.833333333333333333 

#define SIGNint(x)  ( (x) > 0  ?   1  : ( (x) < 0 ?   -1   :  0  ) )
#define SIGN(x)     ( (x) >= 0  ?  1.0 : ( (x) < 0 ? (-1.0) : 0.0 ) )
#define SIGNnull(x)     ( (x) > 0   ?  1.0 : ( (x) < 0 ? (-1.0) : 0.0 ) )
#define SIGNplus(x)     ( (x) >= 0  ?  1.0 :  (-1.0) )
#define SIGNminus(x)    ( (x) >  0  ?  1.0 :  (-1.0) )
#define MIN_own(a,b)    ( a < b   ?  (a) :   (b) )
#define MAX_own(a,b)    ( a >= b  ?  (a) :   (b) )
#define ABS(x)      (SIGN(x) * (x))
#define HEAVISIDE(x) ( (x) > 0.0 ? 1 : 0 )
#define SQR(x) ((x) * (x))
#define LIMIT(as,bs,sa,sb,sab,k)  MAX_own(0,sab)*sa*MAX_own(MIN_own(k*as,bs),MIN_own(as,k*bs))

/* #define SMOOTHSIGN(x,xref)     ( atan(x/xref)/PI ) */
#define SMOOTHSIGN(x,xref)     SIGN(x)
#define UPWIND_EPS 0.01

#endif /* MATH_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: math.h,v $
 Revision 1.1  1998/03/01 18:43:34  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/


