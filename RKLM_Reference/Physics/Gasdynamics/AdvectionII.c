//
//  AdvectionII.c
//  
//
//  Created by Klein, Rupert on 05.11.18.
//
#include <assert.h>
#include <string.h>
#include <math.h>



#include "AdvectionII.h"

void advect_mpdata(
                   ConsVars *Sol, 
                   ConsVars* flux[3],
                   double* force[3],
                   const double dt, 
                   const ElemSpaceDiscr* elem,
                   const enum FluxesFrom adv_fluxes_from, 
                   const enum MUSCL_ON_OFF muscl_on_off, 
                   const enum No_of_Strang_Sweeps no_of_sweeps,
                   const int odd)
{
    
    /* reconstructing Piotr's FORTRAN routine here in C */
    assert(elem->ndim == 2);
    
    int iord  = 2;
    int nonos = 1;
    
    double *U1 = flux[0]->rhoY; 
    double *U2 = flux[1]->rhoY; 
    
    
    
}
