/*******************************************************************************
 File:   numerical_flux.h
 *******************************************************************************/
#ifndef NUMERICAL_FLUX_H
#define NUMERICAL_FLUX_H

#include "variable.h"


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

void store_advective_fluxes(VectorField* adv_flux, 
                            const ConsVars* flux[3], 
                            const ElemSpaceDiscr* elem);

void add_advective_fluxes(VectorField* fd,
                          const ConsVars* flux[3], 
                          const int sign, 
                          const VectorField* ff,
                          const ElemSpaceDiscr* elem);

void update_advective_fluxes(ConsVars* flux[3], 
                             const VectorField* adv_flux,  
                             const ElemSpaceDiscr* elem,  
                             const NodeSpaceDiscr* node,
                             const double dt);

void recompute_advective_fluxes(ConsVars* flux[3], 
                                const ConsVars* Sol, 
                                const ElemSpaceDiscr* elem);

void Advective_Fluxes(VectorField* adv_flux, 
                      const ConsVars* Sol, 
                      const ElemSpaceDiscr* elem);

void Advective_Fluxes_x(double* rhoYu, 
                        const ConsVars* Sol, 
                        const ElemSpaceDiscr* elem, 
                        const int SplitStep);

#endif /* NUMERICAL_FLUX_H */
