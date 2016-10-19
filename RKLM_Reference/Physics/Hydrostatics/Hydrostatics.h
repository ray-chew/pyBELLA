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

void Hydrostatics_State(MPV* mpv, const ElemSpaceDiscr* elem);

#endif /* Hydrostatics_h */
