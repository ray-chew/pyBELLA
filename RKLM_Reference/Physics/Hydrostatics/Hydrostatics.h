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
#include "mpv.h"

/* ================================================================================== */

void Hydrostatics_Column(States* HydroState, 
                         States* HydroState_n, 
                         double* Y, 
                         double* Yn, 
                         const ElemSpaceDiscr* elem, 
                         const NodeSpaceDiscr* node);

void Hydrostatics_State(MPV* mpv, 
                        const ElemSpaceDiscr* elem, 
                        const NodeSpaceDiscr* node);

void Hydrostatic_Exner_pressure(
                                double *pi, 
                                const double pi0, 
                                const double *S, 
                                const double S0,
                                const double dh,
                                const int n, 
                                const int ig);

void Hydrostatic_Initial_Pressure(ConsVars* Sol, 
                                  MPV* mpv,
                                  const ElemSpaceDiscr *elem,
                                  const NodeSpaceDiscr *node);

#endif /* Hydrostatics_h */
