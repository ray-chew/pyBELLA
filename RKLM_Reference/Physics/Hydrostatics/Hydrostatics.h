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

void Hydrostatics_State(MPV* mpv, double *Yinvbg, const ElemSpaceDiscr* elem);

void Hydrostatic_Exner_pressure(
                                double *pi, 
                                const double pi0, 
                                const double *Yinv, 
                                const double Yinv0,
                                const double dh,
                                const int n, 
                                const int ig);
#endif /* Hydrostatics_h */
