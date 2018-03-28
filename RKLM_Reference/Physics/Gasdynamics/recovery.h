/*******************************************************************************
 File:   recovery.h
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

void recovery(States* Lefts,
              States* Rights,
              States* Sol,
              ConsVars* Fluxes,
              const double* force,
              const double lambda_input, 
              const double dt_force, 
              const int nmax,
              const enum FluxesFrom adv_fluxes_from, 
              const enum MUSCL_ON_OFF muscl_on_off) ;

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void recovery_malloc(const int size);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void recovery_free( void );

#endif /* RECOVERY_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: recovery.h,v $
 Revision 1.1  1998/03/01 18:43:35  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
