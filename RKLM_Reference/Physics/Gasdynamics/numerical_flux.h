/*******************************************************************************
 File:   numerical_flux.h
 Author: Nicola
 Date:   Wed Feb 25 13:21:46 WET 1998
 *******************************************************************************/
#ifndef NUMERICAL_FLUX_H
#define NUMERICAL_FLUX_H

#include "variable.h"


/*------------------------------------------------------------------------------
 HLLE (Harten, Lax, van Leer and Einfeldt) numerical flux for the System I*
 
 Notice: this is Thomas Schneider's extension of the HLLE (Harten, Lax, van 
 Leer and Einfeldt) numerical flux for the System I*. For 
 ------------------------------------------------------------------------------*/
void hllestar(
			  ConsVars* Fluxes, 
			  States* Lefts, 
			  States* Rights,
              States* Sol,
              const double lambda,
			  const int n);

void Advective_Fluxes(VectorField* adv_flux, 
                      const ConsVars* Sol, 
                      const ElemSpaceDiscr* elem);

void Advective_Fluxes_x(double* rhoYu, 
                        const ConsVars* Sol, 
                        const ElemSpaceDiscr* elem, 
                        const int SplitStep);

#endif /* NUMERICAL_FLUX_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: numerical_flux.h,v $
 Revision 1.1  1998/03/01 18:43:35  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
