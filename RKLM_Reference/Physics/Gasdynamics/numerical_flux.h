/*******************************************************************************
 File:   numerical_flux.h
 *******************************************************************************/
#ifndef NUMERICAL_FLUX_H
#define NUMERICAL_FLUX_H

#include "variable.h"
#include "boundary.h"

/*------------------------------------------------------------------------------
 plain upwind flux
 ------------------------------------------------------------------------------*/
void hllestar(
			  ConsVars* Fluxes, 
			  States* Lefts, 
			  States* Rights,
              States* Sol,
              const double lambda,
              const int n,
              const enum FluxesFrom adv_fluxes_from);

void recompute_advective_fluxes(ConsVars* flux[3], 
                                const ConsVars* Sol, 
                                const BDRY* bdry,
                                const ElemSpaceDiscr* elem,
                                const double dt);

#endif /* NUMERICAL_FLUX_H */
