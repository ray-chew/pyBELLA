/*******************************************************************************
 File:   recovery.h
 Author: Nicola
 Date:   Wed Feb 25 07:58:40 CET 1998
 *******************************************************************************/
#ifndef RECOVERY_H
#define RECOVERY_H

#include "variable.h"


/*------------------------------------------------------------------------------
 Input:  Sol    = conservative and primitive Variables
 lambda = velocity of signal
 nmax   = value for cache controlled calculation
 
 Output: Lefts  =
 Rights =
 ------------------------------------------------------------------------------*/

void recovery_gravity(
					  States* Lefts, 
					  States* Rights,
					  double* gravity_source, 
                      double* buoy,
                      double* Yinv_ave,
                      double* Yinvbg,
					  const double gravity_strength,
					  States* Sol, 
					  double* S2,
					  double* p2,
                      const double* dp2,
					  const double lambda, 
					  const int nmax,
                      const int stage);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void pressure_gradient_and_gravity(
                                   double* gravity_source,
                                   const double strength,
                                   States* Sol,
                                   double* S2,
                                   double* p2,
                                   const double lambda, 
                                   const int nmax);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void recovery_malloc(const int size);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void recovery_free();

#endif /* RECOVERY_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: recovery.h,v $
 Revision 1.1  1998/03/01 18:43:35  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
