/*******************************************************************************
 File:   flux_correction.h
 Author: Nicola
 Date:   Fri Mar 13 14:42:38 WET 1998
 *******************************************************************************/
#ifndef FLUX_CORRECTION_H
#define FLUX_CORRECTION_H

#include "variable.h"
#include "mpv.h"


/*------------------------------------------------------------------------------
 Remark: this is a tentative implementation of the flux correction algorithm
 described in the hand notes 11-15 March '98. The following simplifications 
 havebeen done:
 
 - interface states are not decoded from the correspondent fluxes. Instead
 the average value of the last computed left and right states is used.
 
 - in the construction of the right hand side the total energy is assumed to
 be constant in time and the terms (p_{nl} - p) are neglected. 
 
 - in the correction all the blue terms are neglected.
 
 ------------------------------------------------------------------------------*/

void flux_correction(
					 ConsVars* flux[3],
					 VectorField* buoy,
					 const ElemSpaceDiscr* elem,
					 ConsVars* Sol, 
					 ConsVars* Sol0, 
					 const double t,
					 const double dt,
					 const double theta,
                     const int step);

void operator_coefficients(
                           double* hplus[3], 
                           double* wcenter, 
                           double* hS, 
                           const ElemSpaceDiscr* elem,
                           const ConsVars* Sol,
                           const ConsVars* Sol0,
                           const MPV* mpv,
                           const double dt); 


#endif /* FLUX_CORRECTION_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
