//
//  FastSubsystem.h
//  RKLM
//
//  Created by Klein, Rupert on 25/08/16.
//  Copyright © 2016 Klein, Rupert. All rights reserved.
//

#ifndef FastSubsystem_h
#define FastSubsystem_h

#include "Common.h"
#include "variable.h"
#include "mpv.h"

void Predictor_for_fast_Subsystem(ConsVars *Sol, 
                                  MPV *mpv, 
                                  const ElemSpaceDiscr *elem, 
                                  const NodeSpaceDiscr *node, 
                                  const double dt);

#endif /* FastSubsystem_h */
