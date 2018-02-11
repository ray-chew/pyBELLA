/*******************************************************************************
 File:   numerical_flux.c
 *******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdarg.h>
#include <math.h>
#include "Common.h"
#include "math_own.h"
#include "Eos.h"   
#include "thermodynamic.h"
#include "variable.h"
#include "userdata.h"
#include "mpv.h"
#include "error.h"
#include "numerical_flux.h"
#include "boundary.h"
#include "boundary.h"
#include "flux_correction.h"
#include "io.h"

/*------------------------------------------------------------------------------
 plain upwind flux
 ------------------------------------------------------------------------------*/
void hllestar(
              ConsVars* Fluxes, 
              States* Lefts, 
              States* Rights,
              States* Sol,
              const double lambda,
              const int n,
              const enum FluxesFrom adv_fluxes_from) {
    
    extern User_Data ud;
    
    double rhol, ul, vl, wl, pl, rhoul, Hl, Yl;
    double rhor, ur, vr, wr, pr, rhour, Hr, Yr;
    double Xl[NSPEC], Xr[NSPEC];
    double upwind, upl, upr;
    int i, nsp;
    
    const double given_flux = (adv_fluxes_from == FLUX_EXTERNAL ? 1.0 : 0.0);
    
    primitives(Lefts,  1, n - 1);
    primitives(Rights, 1, n - 1);
    
    /* sanity values for fluxes at "0", "n-1" and "n-2" */
    Fluxes->rho[0]  = Fluxes->rho[n-2]  = Fluxes->rho[n-1]  = 0.0;
    Fluxes->rhou[0] = Fluxes->rhou[n-2] = Fluxes->rhou[n-1] = 0.0;
    Fluxes->rhov[0] = Fluxes->rhov[n-2] = Fluxes->rhov[n-1] = 0.0;
    Fluxes->rhow[0] = Fluxes->rhow[n-2] = Fluxes->rhow[n-1] = 0.0;
    Fluxes->rhoe[0] = Fluxes->rhoe[n-2] = Fluxes->rhoe[n-1] = 0.0;
    Fluxes->rhoZ[PRES][0] = Fluxes->rhoZ[PRES][n-2] = Fluxes->rhoZ[PRES][n-1] = 0.0;

    for (nsp = 0; nsp < ud.nspec; nsp++) {
        Fluxes->rhoX[nsp][0] = Fluxes->rhoX[nsp][n-2] = Fluxes->rhoX[nsp][n-1] = 0.0;
    }

    for(i = 1; i < n - 2; i++) {	
        rhol  = Lefts->rho[i];
        ul    = Lefts->u[i];
        vl    = Lefts->v[i];
        wl    = Lefts->w[i];
        pl    = Lefts->p[i];
        Yl    = Lefts->Y[i];
        
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Xl[nsp]    = Lefts->X[nsp][i];
        }
        
        rhoul = Lefts->rhou[i];
        Hl = Lefts->rhoe[i] + pl;
        
        rhor  = Rights->rho[i+1];
        ur    = Rights->u[i+1];
        vr    = Rights->v[i+1];
        wr    = Rights->w[i+1];
        pr    = Rights->p[i+1];
        Yr    = Rights->Y[i+1];
        
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Xr[nsp]    = Rights->X[nsp][i+1];
        }
        
        rhour = Rights->rhou[i+1];
        Hr    = Rights->rhoe[i+1] + pr;
        
        Fluxes->rhoY[i] = given_flux * Fluxes->rhoY[i] + (1.0-given_flux) * 0.25 * (rhol*Yl+rhor*Yr)*(ul + ur);
        
        upwind = 0.5 * ( 1.0 + SIGN(Fluxes->rhoY[i]));
        
        upl    = upwind / Yl;
        upr    = (1.0 - upwind) / Yr;
        Fluxes->rhou[i] = Fluxes->rhoY[i] * (upl * ul  + upr * ur) ;
                
        Fluxes->rho[i]  = Fluxes->rhoY[i] * (upl * 1.0 + upr * 1.0);
        Fluxes->rhoe[i] = Fluxes->rhoY[i] * (upl * Hl  + upr * Hr) ;
        Fluxes->rhov[i] = Fluxes->rhoY[i] * (upl * vl  + upr * vr) ;
        Fluxes->rhow[i] = Fluxes->rhoY[i] * (upl * wl  + upr * wr) ;
        
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Fluxes->rhoX[nsp][i] = Fluxes->rhoY[i] * (upl * Xl[nsp]  + upr * Xr[nsp]) ;
        }        
    }
}

/*------------------------------------------------------------------------------
 store advective flux
 ------------------------------------------------------------------------------*/
void store_advective_fluxes(VectorField* adv_flux_full, 
                            const ConsVars* flux[3], 
                            const ElemSpaceDiscr* elem)
{
    
    /* copy the advective fluxes */
    memcpy(adv_flux_full->x, flux[0]->rhoY, elem->nfx*sizeof(double));
    if (elem->ndim > 1) memcpy(adv_flux_full->y, flux[1]->rhoY, elem->nfy*sizeof(double));
    if (elem->ndim > 2) memcpy(adv_flux_full->z, flux[2]->rhoY, elem->nfz*sizeof(double));
}

/*------------------------------------------------------------------------------
 store advective flux difference
 ------------------------------------------------------------------------------*/
void add_advective_fluxes(VectorField* fd,
                          const ConsVars* flux[3], 
                          const int sign, 
                          const VectorField* ff,
                          const ElemSpaceDiscr* elem)
{
    for (int i=0; i<elem->nfx; i++) fd->x[i] = flux[0]->rhoY[i] + sign*ff->x[i]; 
    if (elem->ndim >1) for (int i=0; i<elem->nfy; i++) fd->y[i] = flux[1]->rhoY[i] + sign*ff->y[i]; 
    if (elem->ndim >2) for (int i=0; i<elem->nfz; i++) fd->z[i] = flux[2]->rhoY[i] + sign*ff->z[i]; 
}

/*------------------------------------------------------------------------------
 store advective flux
 ------------------------------------------------------------------------------*/
void recompute_advective_fluxes(ConsVars* flux[3], 
                      const ConsVars* Sol, 
                      const ElemSpaceDiscr* elem)
{
    /* recompute advective flux at fixed time level from cell averages */
    switch (elem->ndim) {
        case 1: {
            for(int i=1; i<elem->icx; i++) {
                double u_c    = Sol->rhou[i]/Sol->rho[i];
                double u_m    = Sol->rhou[i-1]/Sol->rho[i-1];
                double rhoY_c = Sol->rhoY[i];
                double rhoY_m = Sol->rhoY[i-1];
                flux[0]->rhoY[i] = 0.25*(u_c+u_m)*(rhoY_c+rhoY_m);
            }
            break;
        } 
            
        case 2: {
            int icx = elem->icx;
            int icy = elem->icy;
            int ifx = elem->ifx;
            int ify = elem->ify;
            for (int j=1; j<icy; j++) {
                int ncj  = j*icx;
                int nfxj = j*ifx;
                int nfyj = j;
                for(int i=1; i<icx; i++) {
                    int ncij  = ncj  + i;
                    int nfxij = nfxj + i;
                    int nfyij = nfyj + i*ify;
                    double u_c     = Sol->rhou[ncij]/Sol->rho[ncij];
                    double u_m     = Sol->rhou[ncij-1]/Sol->rho[ncij-1];
                    double rhoY_c  = Sol->rhoY[ncij];
                    double rhoY_mx = Sol->rhoY[ncij-1];
                    
                    double v_c     = Sol->rhov[ncij]/Sol->rho[ncij];
                    double v_m     = Sol->rhov[ncij-icx]/Sol->rho[ncij-icx];
                    double rhoY_my = Sol->rhoY[ncij-icx];
                    
                    flux[0]->rhoY[nfxij] = 0.25*(u_c+u_m)*(rhoY_c+rhoY_mx);
                    flux[1]->rhoY[nfyij] = 0.25*(v_c+v_m)*(rhoY_c+rhoY_my);
                }
            }
            break;
        }
            
        case 3: {
            ERROR("recompute_advective_fluxes() not implemented in 3D yet.");
            break;
        }
    }
}


/*------------------------------------------------------------------------------
 update advective flux
 ------------------------------------------------------------------------------*/
void update_advective_fluxes(ConsVars* flux[3], 
                             const VectorField* adv_flux, 
                             const ElemSpaceDiscr* elem, 
                             const NodeSpaceDiscr* node,
                             const double dt)
{
    extern User_Data ud;
    extern double *W0;
    double *rhs = W0;
    
    /* advective fluxes in adv_flux get updated by increments in flux and
       stored in the latter field
     */
    for (int nf=0; nf<elem->nfx; nf++) {
        flux[0]->rho[nf]  = 0.0;
        flux[0]->rhou[nf] = 0.0;
        flux[0]->rhov[nf] = 0.0;
        flux[0]->rhow[nf] = 0.0;
        flux[0]->rhoe[nf] = 0.0;
        for (int nsp=0; nsp<ud.nspec; nsp++) {
            flux[0]->rhoX[nsp][nf] = 0.0;
        }
        flux[0]->rhoY[nf] += adv_flux->x[nf];
    }            
    if (elem->ndim > 1) {
        for (int nf=0; nf<elem->nfy; nf++) {
            flux[1]->rho[nf]  = 0.0;
            flux[1]->rhou[nf] = 0.0;
            flux[1]->rhov[nf] = 0.0;
            flux[1]->rhow[nf] = 0.0;
            flux[1]->rhoe[nf] = 0.0;
            for (int nsp=0; nsp<ud.nspec; nsp++) {
                flux[1]->rhoX[nsp][nf] = 0.0;
            }
            flux[1]->rhoY[nf] += adv_flux->y[nf];
        }            
    }
    if (elem->ndim > 2) {
        for (int nf=0; nf<elem->nfz; nf++) {
            flux[2]->rho[nf]  = 0.0;
            flux[2]->rhou[nf] = 0.0;
            flux[2]->rhov[nf] = 0.0;
            flux[2]->rhow[nf] = 0.0;
            flux[2]->rhoe[nf] = 0.0;
            for (int nsp=0; nsp<ud.nspec; nsp++) {
                flux[2]->rhoX[nsp][nf] = 0.0;
            }
            flux[2]->rhoY[nf] += adv_flux->z[nf];
        }            
    }
#if 1
    extern User_Data ud;
    double rhsmax = controlled_variable_flux_divergence(rhs, (const ConsVars**)flux, dt, elem);
    printf("rhsmax = %e", rhsmax);
    FILE *prhs2file = NULL;
    char fn2[120], fieldname2[90];
    sprintf(fn2, "%s/rhs_cells/rhs_post.hdf", ud.file_name);
    sprintf(fieldname2, "rhs_post");
    WriteHDF(prhs2file, elem->icx, elem->icy, elem->icz, elem->ndim, rhs, fn2, fieldname2);
    sprintf(fn2, "%s/fluxes/flux_rhoY_x.hdf", ud.file_name);
    sprintf(fieldname2, "flux_rhoY_x");
    WriteHDF(prhs2file, elem->ifx, elem->icy, elem->icz, elem->ndim, flux[0]->rhoY, fn2, fieldname2);
    sprintf(fn2, "%s/fluxes/flux_rhoY_y.hdf", ud.file_name);
    sprintf(fieldname2, "flux_rhoY_y");
    WriteHDF(prhs2file, elem->ify, elem->icx, elem->icz, elem->ndim, flux[1]->rhoY, fn2, fieldname2);
#endif
    
    
}
