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

MPV* mpv;                

/* ========================================================================== */

void initialize_MG_projection(
							  const int inx,
							  const int iny,
							  const int inz,
							  const double x0, 
							  const double x1, 
							  const double y0, 
							  const double y1, 
							  const double z0, 
							  const double z1,
							  const enum BdryType left, 
							  const enum BdryType right, 
							  const enum BdryType bottom, 
							  const enum BdryType top, 
							  const enum BdryType back, 
							  const enum BdryType front)
{
	/* allocate multigrid structure for first and second projections */
    Grid* grid;
	int i, nc, nn;
    	
    grid = Grid_new(inx, iny, inz, 
                    x0, x1, y0, y1, z0, z1, 
                    left, right, bottom, top, back, front); 
	
    mpv = MPV_new(grid);
    Grid_free(grid); 
	
	nc = mpv->Level[0]->elem->nc;	
	nn = mpv->Level[0]->node->nc;	
	mpv->p0  = 1.0;
	mpv->p00 = 1.0;
	
    mpv->p2_cells  = (double*)malloc(nc*sizeof(double));
    mpv->dp2_cells = (double*)malloc(nc*sizeof(double));
    mpv->p2_nodes  = (double*)malloc(nn*sizeof(double));
    mpv->dp2_nodes = (double*)malloc(nn*sizeof(double));
    mpv->eta_hyp0  = (double*)malloc(nn*sizeof(double));
    mpv->eta_hyp   = (double*)malloc(nn*sizeof(double));
    mpv->eta       = (double*)malloc(nn*sizeof(double));
    mpv->eta0      = (double*)malloc(nn*sizeof(double));
	
    for(i=0; i<nc; i++){
		mpv->p2_cells[i]  = 0.0;
		mpv->dp2_cells[i] = 0.0;
	}

    for(i=0; i<nn; i++){
		mpv->p2_nodes[i]  = 0.0;
		mpv->dp2_nodes[i] = 0.0;
        mpv->eta_hyp0[i]  = 0.0;
        mpv->eta_hyp[i]   = 0.0;
        mpv->eta0[i]      = 0.0;
        mpv->eta[i]       = 0.0;
	}
}

/* ========================================================================== */

void terminate_MG_projection(void)
{
    MPV_free(mpv);
}

/* ========================================================================== */

void initialize_HydroState(
						   int nc,
                           int nn) {
    
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

MPV* MPV_new(const Grid* base_grid)
{	
    MPV* mpv=(MPV*)malloc(sizeof(MPV));
        
    mpv->Level    = (LowMachMGLevel**)malloc(sizeof(LowMachMGLevel*));
    mpv->Level[0] = LowMachMGLevel_new(0, base_grid);
	
    return(mpv);
}

/* ========================================================================== */

LowMachMGLevel* LowMachMGLevel_new(const int level, const Grid* base_grid)
{
    int dir, ndim, inx0, iny0, inz0, inxi, inyi, inzi, nn, nf, nfs[3];
    Grid* grid;
	
    LowMachMGLevel* mg_level=(LowMachMGLevel*)malloc(sizeof(LowMachMGLevel));
    
    mg_level->level_number = level;
	
    mg_level->elem  = (ElemSpaceDiscr*)malloc(sizeof(ElemSpaceDiscr));
    mg_level->node  = (NodeSpaceDiscr*)malloc(sizeof(NodeSpaceDiscr));
	
    inx0 = base_grid->inx;
    iny0 = base_grid->iny;
    inz0 = base_grid->inz;
	
    inxi = 1+(inx0-1)/integer_power(2,level);
    inyi = 1+(iny0-1)/integer_power(2,level); 
    inzi = 1+(inz0-1)/integer_power(2,level);
	
    grid = Grid_new(
					inxi, 
					inyi, 
					inzi, 
					base_grid->x0, 
					base_grid->x1, 
					base_grid->y0, 
					base_grid->y1, 
					base_grid->z0, 
					base_grid->z1, 
					base_grid->left,
					base_grid->right,
					base_grid->bottom,
					base_grid->top,
					base_grid->back,
					base_grid->front); 
	
    mg_level->elem = ElemSpaceDiscr_new(grid);
	mg_level->node = NodeSpaceDiscr_new(grid);
	
	/* #define LOW_MACH_SCALINGS  */
#ifdef LOW_MACH_SCALINGS
	mg_level->elem->scale_factor = pow(2.0,level);
	mg_level->node->scale_factor = pow(2.0,level);
#else
	mg_level->elem->scale_factor = 1.0;
	mg_level->node->scale_factor = 1.0;
#endif
	
    /* cell-centered pressure fields */
    nn = mg_level->node->nc;
    
    mg_level->p_coarse = (double*)malloc(nn*sizeof(double));
    mg_level->p        = (double*)malloc(nn*sizeof(double));
    mg_level->rhs      = (double*)malloc(nn*sizeof(double));
    mg_level->diaginv  = (double*)malloc(nn*sizeof(double)); 
	
    /* cell-face-centered coefficients */
    ndim = mg_level->elem->ndim;
    
    nfs[0] = mg_level->node->nfx;       
    if(ndim > 1) nfs[1] = mg_level->node->nfy;
    if(ndim > 2) nfs[2] = mg_level->node->nfz;
	
    for(dir = 0; dir < ndim; dir++) {
		nf = nfs[dir];
		
		mg_level->wplus[dir]   = (double*)malloc( nf * sizeof(double));
        /*
		mg_level->wminus[dir]  = (double*)malloc(nf*sizeof(double));
		mg_level->lwplus[dir]  = (double*)malloc(nf*sizeof(double));
		mg_level->lwminus[dir] = (double*)malloc(nf*sizeof(double));
         */
    }
	mg_level->wcenter = (double*)malloc(nn*sizeof(double));
	mg_level->wgrav   = (double*)malloc(nn*sizeof(double));
	
    Grid_free(grid); 
	
    return(mg_level);
}


/* ========================================================================== */

void MPV_free(MPV* var)
{    
    LowMachMGLevel_free(var->Level[0]);
	
    free(var);
}

/* ========================================================================== */

void LowMachMGLevel_free(LowMachMGLevel* var)
{
    int dir, ndim, nfs[3];
	
    ndim = var->elem->ndim;
    nfs[0] = var->elem->nfx;
    if(ndim > 1) nfs[1] = var->elem->nfy;
    if(ndim > 2) nfs[2] = var->elem->nfz;
	
    for(dir = 0; dir < ndim; dir++) {		
		free(var->wplus[dir]);
        /*
		free(var->wminus[dir]);
		free(var->lwplus[dir]);
		free(var->lwminus[dir]);
         */
	}
	
    free(var->p_coarse);
    free(var->p);
    free(var->rhs);
	
    ElemSpaceDiscr_free(var->elem);
	
    free(var);
}


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: mpv.c,v $
 Revision 1.1  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
