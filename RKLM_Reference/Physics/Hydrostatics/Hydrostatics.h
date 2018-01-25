//
//  Hydrostatics.h
//  RKLM_Reference
//
//  Created by Klein, Rupert on 19/10/16.
//  Copyright © 2016 Klein, Rupert. All rights reserved.
//

#ifndef Hydrostatics_h
#define Hydrostatics_h

#include <stdio.h>
#include "space_discretization.h"
#include "mpv.h"

/* ================================================================================== */

void Hydrostatics_Column(States* HydroState, 
                         States* HydroState_n, 
                         double* Y, 
                         double* Yn, 
                         const ElemSpaceDiscr* elem, 
                         const NodeSpaceDiscr* node);

void Hydrostatics_State(MPV* mpv, double *Sbg, const ElemSpaceDiscr* elem);

void Hydrostatic_Exner_pressure(
                                double *pi, 
                                const double pi0, 
                                const double *S, 
                                const double S0,
                                const double dh,
                                const int n, 
                                const int ig);
#endif /* Hydrostatics_h */
