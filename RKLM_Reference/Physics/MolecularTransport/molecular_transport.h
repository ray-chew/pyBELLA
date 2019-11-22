//
//  molecular_transport.h
//  RKLM_Reference
//
//  Created by Klein, Rupert on 16.02.19.
//  Copyright © 2019 Klein, Rupert. All rights reserved.
//

#ifndef molecular_transport_h
#define molecular_transport_h

#include <stdio.h>
#include "variable.h"
#include "kgrid.h"

void molecular_transport(ConsVars* Sol, 
                         double* diss,
                         const ElemSpaceDiscr* elem, 
                         const double dt);

void diss_to_rhoY(ConsVars* Sol,
                  const double* diss,
                  const ElemSpaceDiscr* elem,
                  const NodeSpaceDiscr* node,
                  const double dt);

#endif /* molecular_transport_h */
