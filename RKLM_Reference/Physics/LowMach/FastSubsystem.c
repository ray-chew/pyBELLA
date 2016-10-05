//
//  FastSubsystem.c
//  RKLM
//
//  Created by Klein, Rupert on 25/08/16.
//  Copyright © 2016 Klein, Rupert. All rights reserved.
//

#include "FastSubsystem.h"

void Predictor_for_fast_Subsystem(ConsVars *Sol, 
                                  MPV *mpv, 
                                  const ElemSpaceDiscr *elem, 
                                  const NodeSpaceDiscr *node, 
                                  const double dt)
{
    /*
     Following Piotr Smolarkiewicz' time integration strategy, 
     this routine updates all variables by adding 
     $dt/2 R^n$, where $R^n$ denotes all non-advective terms 
     in the governing equations evaluated at the starting point
     of a time step.
     */
    
    
    
}
