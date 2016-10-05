/*******************************************************************************
 File:   thermodynamic.c
 Author: Nicola
 Date:   Thu Feb 26 16:17:40 WET 1998 
 *******************************************************************************/

#include "Common.h"
#include "userdata.h"
#include "thermodynamic.h"

/* ========================================================================= */

void Thermodynamic_init(Thermodynamic* th, const User_Data* ud) {
	
	const double g = ud->gamm;
	
	/* ideal gas stuff */
    th->gamm = g;
	th->gamminv = 1.0 / g;
	th->gm1 = g - 1.0;
	th->gm1inv = 1.0 / (g - 1.0);
	th->Gamma = (g-1.0)/g;
	th->Gammainv = g / (g-1.0);
    
    /* moisture stuff */
    th->Rg_over_Rv = ud->Rg_over_Rv;
    th->Q          = ud->Q;
}

/* ========================================================================= */

void set_Y_to_unity(ConsVars* Sol, ElemSpaceDiscr* elem) {
    
    for (int j=0; j<elem->icy; j++) {int mc = j*elem->icx;
        for (int i=0; i<elem->icx; i++) {
            int nc = mc + i;
            Sol->rho[nc] = Sol->rhoY[nc];
        }
    }
    
}


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: thermodynamic.c,v $
 Revision 1.1  1998/03/01 18:43:36  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
