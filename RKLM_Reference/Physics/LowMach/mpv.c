/*******************************************************************************
 File:   mpv.c
 Author: Nicola
 Date:   Fri Mar  6 16:11:51 WET 1998
 *******************************************************************************/
#include <math.h>
#include "Common.h"
#include "thermodynamic.h"
#include "mpv.h"
/* #include "Hydrostatics.h"*/
#include "kgrid.h"
#include "error.h"
#include "SimpleUtilities.h"

/* ========================================================================== */

void HydroState_init(
                     MPV* mpv,
                     const int nc,
                     const int nn) {
    
    extern User_Data ud;
    int nsp;
    
    mpv->HydroState        = (States*)malloc(sizeof(States));
    mpv->HydroState_n      = (States*)malloc(sizeof(States));
    
    mpv->HydroState->p0     = (double*)malloc(nc*sizeof(double));
    mpv->HydroState->pi0    = (double*)malloc(nc*sizeof(double));
    mpv->HydroState->p20    = (double*)malloc(nc*sizeof(double));
    mpv->HydroState->rho0   = (double*)malloc(nc*sizeof(double));
    mpv->HydroState->rhoY0  = (double*)malloc(nc*sizeof(double));
    mpv->HydroState->S0     = (double*)malloc(nc*sizeof(double));
    mpv->HydroState->S10    = (double*)malloc(nc*sizeof(double));
    mpv->HydroState->Y0     = (double*)malloc(nc*sizeof(double));
    for (nsp = 0; nsp < ud.nspec; nsp++) {
        mpv->HydroState->X[nsp]    = (double*)malloc(nc*sizeof(double));
    }
    
    mpv->HydroState_n->p0     = (double*)malloc(nn*sizeof(double));
    mpv->HydroState_n->pi0    = (double*)malloc(nn*sizeof(double));
    mpv->HydroState_n->p20    = (double*)malloc(nn*sizeof(double));
    mpv->HydroState_n->rho0   = (double*)malloc(nn*sizeof(double));
    mpv->HydroState_n->rhoY0  = (double*)malloc(nn*sizeof(double));
    mpv->HydroState_n->S0     = (double*)malloc(nn*sizeof(double));
    mpv->HydroState_n->S10    = (double*)malloc(nn*sizeof(double));
    mpv->HydroState_n->Y0     = (double*)malloc(nn*sizeof(double));
    for (nsp = 0; nsp < ud.nspec; nsp++) {
        mpv->HydroState_n->X[nsp]    = (double*)malloc(nn*sizeof(double));
    }
    
}

/* ========================================================================== */

void HydroState_free(
                     MPV* mpv) 
{
    
    extern User_Data ud;
        
    free(mpv->HydroState->p0   );
    free(mpv->HydroState->pi0  );
    free(mpv->HydroState->p20  );
    free(mpv->HydroState->rho0 );
    free(mpv->HydroState->rhoY0);
    free(mpv->HydroState->S0   );
    free(mpv->HydroState->S10  );
    free(mpv->HydroState->Y0   );
    for (int nsp = 0; nsp < ud.nspec; nsp++) {
        free(mpv->HydroState->X[nsp]);
    }
    
    free(mpv->HydroState_n->p0   );
    free(mpv->HydroState_n->pi0  );
    free(mpv->HydroState_n->p20  );
    free(mpv->HydroState_n->rho0 );
    free(mpv->HydroState_n->rhoY0);
    free(mpv->HydroState_n->S0   );
    free(mpv->HydroState_n->S10  );
    free(mpv->HydroState_n->Y0   );
    for (int nsp = 0; nsp < ud.nspec; nsp++) {
        free(mpv->HydroState_n->X[nsp]);
    }
    
    free(mpv->HydroState  );
    free(mpv->HydroState_n);
}
/* ========================================================================== */

MPV* MPV_new(
             const ElemSpaceDiscr *elem,
             const NodeSpaceDiscr *node)
{
    extern User_Data ud;
    
	/* allocate infrastructure for the semi-implicit pressure solves */
	int i, nc, nn;
    	
    MPV* mpv=(MPV*)malloc(sizeof(MPV));	
	
    nc = elem->nc;	
	nn = node->nc;	
	mpv->p0  = 1.0;
	mpv->p00 = 1.0;
	
    mpv->p2_cells  = (double*)malloc(nc*sizeof(double));
    mpv->dp2_cells = (double*)malloc(nc*sizeof(double));
    mpv->p2_nodes  = (double*)malloc(nn*sizeof(double));
    mpv->dp2_nodes = (double*)malloc(nn*sizeof(double));
	
    for(i=0; i<nc; i++){
		mpv->p2_cells[i]  = 0.0;
		mpv->dp2_cells[i] = 0.0;
	}

    for(i=0; i<nn; i++){
		mpv->p2_nodes[i]  = 0.0;
		mpv->dp2_nodes[i] = 0.0;
	}

    mpv->rhs     = (double*)malloc(nc*sizeof(double));
    mpv->diaginv = (double*)malloc(nc*sizeof(double));
    mpv->wcenter = (double*)malloc(nn*sizeof(double));
    mpv->wgrav   = (double*)malloc(nn*sizeof(double));
    for (int idim = 0; idim < elem->ndim; idim++) {
        mpv->wplus[idim] = (double*)malloc(nn*sizeof(double));
    }

    /* if (ud.g_ref > ud.eps_Machine) */
        HydroState_init(mpv, elem->ic[1], node->ic[1]);

    return mpv;
}

/* ========================================================================== */

void MPV_free( 
              MPV* mpv,
              const ElemSpaceDiscr* elem)
{
    extern User_Data ud;
                
    free(mpv->p2_cells );
    free(mpv->dp2_cells);
    free(mpv->p2_nodes );
    free(mpv->dp2_nodes);
    
    free(mpv->rhs    );
    free(mpv->diaginv);
    free(mpv->wcenter);
    free(mpv->wgrav  );
    for (int idim = 0; idim < elem->ndim; idim++) {
        free(mpv->wplus[idim]);
    }
    
    if (ud.g_ref > ud.eps_Machine) HydroState_free(mpv);
    
    free(mpv);    
}


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: mpv.c,v $
 Revision 1.1  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
