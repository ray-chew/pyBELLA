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
    
    double rhol, ul, vl, wl, pl, Hl, Yl;
    double rhor, ur, vr, wr, pr, Hr, Yr;
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
        
        Hr    = Rights->rhoe[i+1] + pr;
        
        Fluxes->rhoY[i] = given_flux * Fluxes->rhoY[i] + (1.0-given_flux) * 0.25 * (rhol*Yl+rhor*Yr)*(ul + ur);
        
        upwind = 0.5 * ( 1.0 + SIGN(Fluxes->rhoY[i]));
        //upwind = 0.5;
        
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
#define OUTPUT_ADVECTIVE_FLUXES 0
#if OUTPUT_ADVECTIVE_FLUXES
#include "io.h"
#include "memory.h"
static int flux_output_count = 0;
#endif



void recompute_advective_fluxes(ConsVars* flux[3], 
                                const ConsVars* Sol, 
                                const ElemSpaceDiscr* elem,
                                const double dt)
{
    extern User_Data ud;
    
    /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 
     UPWIND_RHOY should be off to revive the October 24 version. 
     Nov. 15, 2018; I tested the UPWIND_RHOY option, but it does 
     not seem to have any effect.
     */
#ifdef UPWIND_RHOY
    double (*limiter[])(const double a, 
                        const double b, 
                        const double k) = {
        None,
        MinMod,
        VanLeer,
        VanLeerS,
        Superbee,
        MonotonizedCentral,
        SwebyMunz,
        Rupe,
        NoSlope
    };

    const enum LimiterType limiter_type_scalars  = ud.limiter_type_scalars;
    const double kz = ud.kz;
    const double lambdax = dt/elem->dx;
    const double lambday = dt/elem->dy;
#endif
    
#ifndef FOURTH_ORDER_ADV_FLUXES
    
#ifdef NODAL_PROJECTION_ONLY
    /* recompute advective flux at fixed time level from cell averages */
    switch (elem->ndim) {
        case 1: {
            for(int i=1; i<elem->icx; i++) {
                double u_c    = Sol->rhou[i]/Sol->rho[i];
                double u_m    = Sol->rhou[i-1]/Sol->rho[i-1];
                double rhoY_c = Sol->rhoY[i];
                double rhoY_m = Sol->rhoY[i-1];
                flux[0]->rhoY[i] = 0.5*(u_c*rhoY_c+u_m*rhoY_m);
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
                    
                    double rhoYu_cm     = Sol->rhoY[ncij-icx]   * Sol->rhou[ncij-icx]/Sol->rho[ncij-icx];
                    double rhoYu_mm     = Sol->rhoY[ncij-1-icx] * Sol->rhou[ncij-1-icx]/Sol->rho[ncij-1-icx];
                    double rhoYu_cc     = Sol->rhoY[ncij]       * Sol->rhou[ncij]/Sol->rho[ncij];
                    double rhoYu_mc     = Sol->rhoY[ncij-1]     * Sol->rhou[ncij-1]/Sol->rho[ncij-1];
                    double rhoYu_cp     = Sol->rhoY[ncij+icx]   * Sol->rhou[ncij+icx]/Sol->rho[ncij+icx];
                    double rhoYu_mp     = Sol->rhoY[ncij-1+icx] * Sol->rhou[ncij-1+icx]/Sol->rho[ncij-1+icx];

                    double rhoYu_m = 0.5*(rhoYu_cm+rhoYu_mm);
                    double rhoYu_c = 0.5*(rhoYu_cc+rhoYu_mc);
                    double rhoYu_p = 0.5*(rhoYu_cp+rhoYu_mp);
                    
                    flux[0]->rhoY[nfxij] = 0.25*(rhoYu_m + 2.0*rhoYu_c + rhoYu_p);

                    double rhoYv_cm     = Sol->rhoY[ncij-1]     * Sol->rhov[ncij-1]/Sol->rho[ncij-1];
                    double rhoYv_mm     = Sol->rhoY[ncij-icx-1] * Sol->rhov[ncij-icx-1]/Sol->rho[ncij-icx-1];
                    double rhoYv_cc     = Sol->rhoY[ncij]       * Sol->rhov[ncij]/Sol->rho[ncij];
                    double rhoYv_mc     = Sol->rhoY[ncij-icx]   * Sol->rhov[ncij-icx]/Sol->rho[ncij-icx];
                    double rhoYv_cp     = Sol->rhoY[ncij+1]     * Sol->rhov[ncij+1]/Sol->rho[ncij+1];
                    double rhoYv_mp     = Sol->rhoY[ncij-icx+1] * Sol->rhov[ncij-icx+1]/Sol->rho[ncij-icx+1];
                    
                    double rhoYv_m = 0.5*(rhoYv_cm+rhoYv_mm);
                    double rhoYv_c = 0.5*(rhoYv_cc+rhoYv_mc);
                    double rhoYv_p = 0.5*(rhoYv_cp+rhoYv_mp);
                    
                    flux[1]->rhoY[nfyij] = 0.25*(rhoYv_m + 2.0*rhoYv_c + rhoYv_p);
                }
            }
            
            break;
        }
            
        case 3: {
            assert(0);  /* 3D variant of NODAL_PROJECTION_ONLY not yet implemented */
            int icx = elem->icx;
            int icy = elem->icy;
            int icz = elem->icz;
            int ifx = elem->ifx;
            int ify = elem->ify;
            int ifz = elem->ifz;
            for (int k=1; k<icz; k++) {
                int nck  = k*icx*icy;
                int nfxk = k*ifx*icy;
                int nfyk = k*ify;
                int nfzk = k;
                for (int j=1; j<icy; j++) {
                    int ncjk  = nck  + j*icx;
                    int nfxjk = nfxk + j*ifx;
                    int nfyjk = nfyk + j;
                    int nfzjk = nfzk + j*ifz*icx;
                    for(int i=1; i<icx; i++) {
                        int ncijk  = ncjk  + i;
                        int nfxijk = nfxjk + i;
                        int nfyijk = nfyjk + i*ify*icz;
                        int nfzijk = nfzjk + i*ifz;
                        
                        double u_c     = Sol->rhou[ncijk]/Sol->rho[ncijk];
                        double u_m     = Sol->rhou[ncijk-1]/Sol->rho[ncijk-1];
                        double rhoY_c  = Sol->rhoY[ncijk];
                        double rhoY_mx = Sol->rhoY[ncijk-1];
                        
                        double v_c     = Sol->rhov[ncijk]/Sol->rho[ncijk];
                        double v_m     = Sol->rhov[ncijk-icx]/Sol->rho[ncijk-icx];
                        double rhoY_my = Sol->rhoY[ncijk-icx];
                        
                        double w_c     = Sol->rhow[ncijk]/Sol->rho[ncijk];
                        double w_m     = Sol->rhow[ncijk-icx*icy]/Sol->rho[ncijk-icx*icy];
                        double rhoY_mz = Sol->rhoY[ncijk-icx*icy];
                        
                        flux[0]->rhoY[nfxijk] = 0.25*(u_c+u_m)*(rhoY_c+rhoY_mx);
                        flux[1]->rhoY[nfyijk] = 0.25*(v_c+v_m)*(rhoY_c+rhoY_my);
                        flux[2]->rhoY[nfzijk] = 0.25*(w_c+w_m)*(rhoY_c+rhoY_mz);
                    }
                }
            }
            break;
        }
    }
#else /* NODAL_PROJECTION_ONLY */
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
                    
#ifdef UPWIND_RHOY
                    double rhoY, rhoY_p, rhoY_pp, rhoY_m, rhoY_mm, s_m, s_p; 
                    double u_p      = Sol->rhou[ncij]/Sol->rho[ncij];
                    double u_m      = Sol->rhou[ncij-1]/Sol->rho[ncij-1];
                    double u        = 0.5 * (u_p+u_m);
                    rhoY_p  = Sol->rhoY[ncij];
                    rhoY_pp = Sol->rhoY[ncij+1];
                    rhoY_m  = Sol->rhoY[ncij-1];
                    rhoY_mm = Sol->rhoY[ncij-2];
                    s_m     = (*limiter[limiter_type_scalars])(rhoY_m - rhoY_mm, rhoY_p - rhoY_m, kz);
                    s_p     = (*limiter[limiter_type_scalars])(rhoY_p - rhoY_m, rhoY_pp - rhoY_p, kz);
                    rhoY    = (u >= 0.0? rhoY_m + 0.5 * (1.0-u*lambdax) * s_m : rhoY_p - 0.5 * (1.0 + u*lambdax) * s_p);

                    flux[0]->rhoY[nfxij] = u * rhoY;
                    
                    double v_p      = Sol->rhov[ncij]/Sol->rho[ncij];
                    double v_m      = Sol->rhov[ncij-icx]/Sol->rho[ncij-icx];
                    double v        = 0.5 * (v_p+v_m);
                    rhoY_p  = Sol->rhoY[ncij];
                    rhoY_pp = Sol->rhoY[ncij+icx];
                    rhoY_m  = Sol->rhoY[ncij-icx];
                    rhoY_mm = Sol->rhoY[ncij-2*icx];
                    s_m     = (*limiter[limiter_type_scalars])(rhoY_m - rhoY_mm, rhoY_p - rhoY_m, kz);
                    s_p     = (*limiter[limiter_type_scalars])(rhoY_p - rhoY_m, rhoY_pp - rhoY_p, kz);
                    rhoY    = (v >= 0.0? rhoY_m + 0.5 * (1.0-v*lambday) * s_m : rhoY_p - 0.5 * (1.0+v*lambday) * s_p);

                    flux[1]->rhoY[nfyij] = v * rhoY;
#else /* UPWIND_RHOY */
                    double u_c     = Sol->rhou[ncij]/Sol->rho[ncij];
                    double u_m     = Sol->rhou[ncij-1]/Sol->rho[ncij-1];
                    double rhoY_c  = Sol->rhoY[ncij];
                    double rhoY_mx = Sol->rhoY[ncij-1];
                    
                    double v_c     = Sol->rhov[ncij]/Sol->rho[ncij];
                    double v_m     = Sol->rhov[ncij-icx]/Sol->rho[ncij-icx];
                    double rhoY_my = Sol->rhoY[ncij-icx];
                    
#if 1
                    /* Nov. 15, 2018; Interesting observation: 
                     The following option gives very good results for the 
                     planetary IGW test, whereas the alternative (set 
                     compiler directive to  #if 0) does not work AT ALL. 
                     I am leaving the two options and the comment, 
                     because they observation is so surprising, that
                     one might need a warning when getting ideas later. 
                     */
                    flux[0]->rhoY[nfxij] = 0.25*(u_c+u_m)*(rhoY_c+rhoY_mx);
                    flux[1]->rhoY[nfyij] = 0.25*(v_c+v_m)*(rhoY_c+rhoY_my);
#else
                    flux[0]->rhoY[nfxij] = 0.5*(rhoY_c * u_c + rhoY_mx * u_m);
                    flux[1]->rhoY[nfyij] = 0.5*(rhoY_c * v_c + rhoY_my * v_m);
#endif
#endif /* UPWIND_RHOY */
                }
            }
            break;
        }
            
        case 3: {
            int icx = elem->icx;
            int icy = elem->icy;
            int icz = elem->icz;
            int ifx = elem->ifx;
            int ify = elem->ify;
            int ifz = elem->ifz;
            for (int k=1; k<icz; k++) {
                int nck  = k*icx*icy;
                int nfxk = k*ifx*icy;
                int nfyk = k*ify;
                int nfzk = k;
                for (int j=1; j<icy; j++) {
                    int ncjk  = nck  + j*icx;
                    int nfxjk = nfxk + j*ifx;
                    int nfyjk = nfyk + j;
                    int nfzjk = nfzk + j*ifz*icx;
                    for(int i=1; i<icx; i++) {
                        int ncijk  = ncjk  + i;
                        int nfxijk = nfxjk + i;
                        int nfyijk = nfyjk + i*ify*icz;
                        int nfzijk = nfzjk + i*ifz;
                        
                        double u_c     = Sol->rhou[ncijk]/Sol->rho[ncijk];
                        double u_m     = Sol->rhou[ncijk-1]/Sol->rho[ncijk-1];
                        double rhoY_c  = Sol->rhoY[ncijk];
                        double rhoY_mx = Sol->rhoY[ncijk-1];
                        
                        double v_c     = Sol->rhov[ncijk]/Sol->rho[ncijk];
                        double v_m     = Sol->rhov[ncijk-icx]/Sol->rho[ncijk-icx];
                        double rhoY_my = Sol->rhoY[ncijk-icx];
                        
                        double w_c     = Sol->rhow[ncijk]/Sol->rho[ncijk];
                        double w_m     = Sol->rhow[ncijk-icx*icy]/Sol->rho[ncijk-icx*icy];
                        double rhoY_mz = Sol->rhoY[ncijk-icx*icy];
                        
                        flux[0]->rhoY[nfxijk] = 0.25*(u_c+u_m)*(rhoY_c+rhoY_mx);
                        flux[1]->rhoY[nfyijk] = 0.25*(v_c+v_m)*(rhoY_c+rhoY_my);
                        flux[2]->rhoY[nfzijk] = 0.25*(w_c+w_m)*(rhoY_c+rhoY_mz);
                    }
                }
            }
            break;
        }
    }
#endif /* NODAL_PROJECTION_ONLY */
#else /* FOURTH_ORDER_ADV_FLUXES */ 
    
    /* TODO: 
     7) Make sure, the fourth order computation of the advective fluxes is implemented
     compatibly with rigid wall boundary conditions. 
     */
    int igx = elem->igx;
    int igy = elem->igy;
    int igz = elem->igz;
    
    int icx = elem->icx;
    int icy = elem->icy;
    int icz = elem->icz;
    
    int ifx = elem->ifx;
    int ify = elem->ify;
    int ifz = elem->ifz;

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
            for (int j=igy; j<icy-1; j++) {
                int ncj  = j*icx;
                int nfxj = j*ifx;
                int nfyj = j;
                
                for(int i=igx; i<icx-1; i++) {
                    int ncij  = ncj  + i;
                    int nfxij = nfxj + i;
                    int nfyij = nfyj + i*ify;
                    
                    double rhoYu_p  = (Sol->rhoY[ncij+1]/Sol->rho[ncij+1])*Sol->rhou[ncij+1];
                    double rhoYu_c  = (Sol->rhoY[ncij]/Sol->rho[ncij])*Sol->rhou[ncij];
                    double rhoYu_m  = (Sol->rhoY[ncij-1]/Sol->rho[ncij-1])*Sol->rhou[ncij-1];
                    double rhoYu_mm = (Sol->rhoY[ncij-2]/Sol->rho[ncij-2])*Sol->rhou[ncij-2];
                    
                    double rhoYv_p  = (Sol->rhoY[ncij+icx]/Sol->rho[ncij+icx])*Sol->rhov[ncij+icx];
                    double rhoYv_c  = (Sol->rhoY[ncij]/Sol->rho[ncij])*Sol->rhov[ncij];
                    double rhoYv_m  = (Sol->rhoY[ncij-icx]/Sol->rho[ncij-icx])*Sol->rhov[ncij-icx];
                    double rhoYv_mm = (Sol->rhoY[ncij-2*icx]/Sol->rho[ncij-2*icx])*Sol->rhov[ncij-2*icx];
                    
                    flux[0]->rhoY[nfxij] = (-rhoYu_p + 7.0*(rhoYu_c+rhoYu_m) - rhoYu_mm)/12.0;
                    flux[1]->rhoY[nfyij] = (-rhoYv_p + 7.0*(rhoYv_c+rhoYv_m) - rhoYv_mm)/12.0;
                }
            }
            break;
        }
            
        case 3: {            
            for (int k=igz; k<icz-1; k++) {
                int nck  = k*icx*icy;
                int nfxk = k*ifx*icy;
                int nfyk = k*ify;
                int nfzk = k;
                for (int j=igy; j<icy-1; j++) {
                    int ncjk  = nck  + j*icx;
                    int nfxjk = nfxk + j*ifx;
                    int nfyjk = nfyk + j;
                    int nfzjk = nfzk + j*ifz*icx;
                    for(int i=igx; i<icx-1; i++) {
                        int ncijk  = ncjk  + i;
                        int nfxijk = nfxjk + i;
                        int nfyijk = nfyjk + i*ify*icz;
                        int nfzijk = nfzjk + i*ifz;
                        
                        double rhoYu_p  = (Sol->rhoY[ncijk+1]/Sol->rho[ncijk+1])*Sol->rhou[ncijk+1];
                        double rhoYu_c  = (Sol->rhoY[ncijk]/Sol->rho[ncijk])*Sol->rhou[ncijk];
                        double rhoYu_m  = (Sol->rhoY[ncijk-1]/Sol->rho[ncijk-1])*Sol->rhou[ncijk-1];
                        double rhoYu_mm = (Sol->rhoY[ncijk-2]/Sol->rho[ncijk-2])*Sol->rhou[ncijk-2];
                        
                        double rhoYv_p  = (Sol->rhoY[ncijk+icx]/Sol->rho[ncijk+icx])*Sol->rhov[ncijk+icx];
                        double rhoYv_c  = (Sol->rhoY[ncijk]/Sol->rho[ncijk])*Sol->rhov[ncijk];
                        double rhoYv_m  = (Sol->rhoY[ncijk-icx]/Sol->rho[ncijk-icx])*Sol->rhov[ncijk-icx];
                        double rhoYv_mm = (Sol->rhoY[ncijk-2*icx]/Sol->rho[ncijk-2*icx])*Sol->rhov[ncijk-2*icx];
                        
                        int diz      = icx*icy;
                        double rhoYw_p  = (Sol->rhoY[ncijk+diz]/Sol->rho[ncijk+diz])*Sol->rhow[ncijk+diz];
                        double rhoYw_c  = (Sol->rhoY[ncijk]/Sol->rho[ncijk])*Sol->rhow[ncijk];
                        double rhoYw_m  = (Sol->rhoY[ncijk-diz]/Sol->rho[ncijk-diz])*Sol->rhow[ncijk-diz];
                        double rhoYw_mm = (Sol->rhoY[ncijk-2*diz]/Sol->rho[ncijk-2*diz])*Sol->rhow[ncijk-2*diz];
                        
                        flux[0]->rhoY[nfxijk] = (-rhoYu_p + 7.0*(rhoYu_c+rhoYu_m) - rhoYu_mm)/12.0;
                        flux[1]->rhoY[nfyijk] = (-rhoYv_p + 7.0*(rhoYv_c+rhoYv_m) - rhoYv_mm)/12.0;
                        flux[2]->rhoY[nfzijk] = (-rhoYw_p + 7.0*(rhoYw_c+rhoYw_m) - rhoYw_mm)/12.0;
                    }
                }
            }
            break;
        }
    }
    
    /* Check wall boundary conditions in the vertical */
    if (ud.bdrytype_min[1] == WALL) {
        for (int k=igz; k<icz-igz; k++) {
            int nfyk = k*ify;
            for(int i=igx; i<icx-1; i++) {
                int nfyik = nfyk + i*ify*icz;
                int j=igy;
                int nfyijk = nfyik + j;
                flux[1]->rhoY[nfyijk] = 0.0;
            }
        }
    }
    if (ud.bdrytype_max[1] == WALL) {
        for (int k=igz; k<icz-igz; k++) {
            int nfyk = k*ify;
            for(int i=igx; i<icx-1; i++) {
                int nfyik = nfyk + i*ify*icz;
                int j=ify-igx-1;
                int nfyijk = nfyik + j;
                flux[1]->rhoY[nfyijk] = 0.0;
            }
        }
    }
#endif /* FOURTH_ORDER_ADV_FLUXES */ 
    
#if OUTPUT_ADVECTIVE_FLUXES
    extern User_Data ud;
    FILE *pfluxfile = NULL;
    char fnx[120], fny[120], fieldname[90];
    if (flux_output_count < 10) {
        sprintf(fnx, "%s/flux_x/flux_x_00%d.hdf", ud.file_name, flux_output_count);
        sprintf(fny, "%s/flux_y/flux_y_00%d.hdf", ud.file_name, flux_output_count);
    } else if(flux_output_count < 100) {
        sprintf(fnx, "%s/flux_x/flux_x_0%d.hdf", ud.file_name, flux_output_count);
        sprintf(fny, "%s/flux_y/flux_y_0%d.hdf", ud.file_name, flux_output_count);
    } else {
        sprintf(fnx, "%s/flux_x/flux_x_%d.hdf", ud.file_name, flux_output_count);
        sprintf(fny, "%s/flux_y/flux_y_%d.hdf", ud.file_name, flux_output_count);
    }
    sprintf(fieldname, "flux_x");    
    WriteHDF(pfluxfile, elem->ifx, elem->icy, elem->icz, elem->ndim, flux[0]->rhoY, fnx, fieldname);
    sprintf(fieldname, "flux_y");   
    WriteHDF(pfluxfile, elem->ify, elem->icx, elem->icz, elem->ndim, flux[1]->rhoY, fny, fieldname);
    
    flux_output_count++;
#endif

}

