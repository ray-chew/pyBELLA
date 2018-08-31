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

void flux_correction(ConsVars* flux[3],
					 ConsVars* Sol, 
					 const ConsVars* Sol0, 
                     const ElemSpaceDiscr* elem,
                     const NodeSpaceDiscr* node,
					 const double t,
					 const double dt_in,
                     const int step);

double controlled_variable_flux_divergence(double* rhs, 
                                           const ConsVars* flux[3],
                                           const ElemSpaceDiscr* elem);

void update_SI_MIDPT_buoyancy(ConsVars* Sol, 
                              const ConsVars* flux[3], 
                              const MPV* mpv,
                              const ElemSpaceDiscr* elem,
                              const double dt);

#endif /* FLUX_CORRECTION_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
