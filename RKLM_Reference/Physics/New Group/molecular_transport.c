//
//  molecular_transport.c
//  RKLM_Reference
//
//  Created by Klein, Rupert on 16.02.19.
//  Copyright © 2019 Klein, Rupert. All rights reserved.
//
#include <math.h>

#include "userdata.h"
#include "thermodynamic.h"
#include "molecular_transport.h"

void molecular_fluxes(ConsVars* flux[3],
                      double* diss,
                      const ConsVars* Sol, 
                      const ElemSpaceDiscr* elem);

/* =================================================================== */

void molecular_transport(ConsVars* Sol, 
                         double* diss,
                         const ElemSpaceDiscr* elem, 
                         const double dt)
{
    /* Here we compute viscous forces, heat conduction effects, and 
     diffusive transport of species. Note that we are not working
     with total energy as a conserved quantity, so that we will
     have to carefully account for the effect of the transport 
     processes on potential temperature to ensure total energy 
     conservation.
     */
    extern User_Data ud;
    extern Thermodynamic th;
    extern ConsVars* flux[3];

    if (ud.mol_trans == FULL_MOLECULAR_TRANSPORT) {
        assert(0); /* full transport branch of code not fully implemented yet */
        
        /* viscous momentum and conductive energy fluxes */
        molecular_fluxes(flux, diss, Sol, elem);
        
        /* Euler forward time stepping */
        switch (elem->ndim) {
            case 1:
                assert(0); 
                break;
                
            case 2: {
                const int icx = elem->icx;
                const int icy = elem->icy;
                const int ifx = elem->ifx;
                const int ify = elem->ify;
                const int igx = elem->igx;
                const int igy = elem->igy;
                
                const double dx = elem->dx;
                const double dy = elem->dy;
                
                for (int j=igy; j<icy-igy; j++) {
                    int ncj  = j*icx;
                    int nfxj = j*ifx;
                    int nfyj = j;
                    
                    for (int i=igx; i<icx-igx; i++) {
                        int ncij  = ncj + i;
                        int nfxij = nfxj + i;
                        int nfyij = nfyj + i*ify;
                        
                        double pold = pow(Sol->rhoY[ncij], th.gamm);
                        double dp   = -dt*th.gm1*( (flux[0]->rhoe[nfxij+1] - flux[0]->rhoe[nfxij])/dx + (flux[1]->rhoe[nfyij+icx] - flux[1]->rhoe[nfyij])/dy + diss[ncij]);
                        Sol->rhoY[ncij]  = pow(pold+dp, th.gamminv);
                        Sol->rhou[ncij] -= dt*( (flux[0]->rhou[nfxij+1] - flux[0]->rhou[nfxij])/dx + (flux[1]->rhou[nfyij+icx] - flux[1]->rhou[nfyij])/dy);
                        Sol->rhov[ncij] -= dt*( (flux[0]->rhov[nfxij+1] - flux[0]->rhov[nfxij])/dx + (flux[1]->rhov[nfyij+icx] - flux[1]->rhov[nfyij])/dy);                    
                    }
                }
            }
                break;
                
            case 3:
                assert(0);
                break;
                
            default:
                break;
        }
    } else if (ud.mol_trans == STRAKA_DIFFUSION_MODEL) {
        
        /* viscous diffusion of velocities and potential temperature */
        molecular_fluxes(flux, diss, Sol, elem);
        
        /* Euler forward time stepping */
        switch (elem->ndim) {
            case 1:
                assert(0); 
                break;
                
            case 2: {
                const int icx = elem->icx;
                const int icy = elem->icy;
                const int ifx = elem->ifx;
                const int ify = elem->ify;
                const int igx = elem->igx;
                const int igy = elem->igy;
                
                const double dx = elem->dx;
                const double dy = elem->dy;
                
                for (int j=igy; j<icy-igy; j++) {
                    int ncj  = j*icx;
                    int nfxj = j*ifx;
                    int nfyj = j;
                    
                    for (int i=igx; i<icx-igx; i++) {
                        int ncij  = ncj + i;
                        int nfxij = nfxj + i;
                        int nfyij = nfyj + i*ify;
                        
                        diss[ncij]      -= dt*Sol->rho[ncij]*( (flux[0]->rhoY[nfxij+1] - flux[0]->rhoY[nfxij])/dx + (flux[1]->rhoY[nfyij+icx] - flux[1]->rhoY[nfyij])/dy);
                        Sol->rhou[ncij] -= dt*Sol->rho[ncij]*( (flux[0]->rhou[nfxij+1] - flux[0]->rhou[nfxij])/dx + (flux[1]->rhou[nfyij+icx] - flux[1]->rhou[nfyij])/dy);
                        Sol->rhov[ncij] -= dt*Sol->rho[ncij]*( (flux[0]->rhov[nfxij+1] - flux[0]->rhov[nfxij])/dx + (flux[1]->rhov[nfyij+icx] - flux[1]->rhov[nfyij])/dy);                    
                    }
                }
            }
                break;
                
            case 3:
                assert(0);
                break;
                
            default:
                break;
        }

    }
}

/* =================================================================== */

void molecular_fluxes(ConsVars* flux[3],
                      double* diss,
                      const ConsVars* Sol, 
                      const ElemSpaceDiscr* elem)
{
    /* Here we compute viscous forces, heat conduction effects, and 
     diffusive transport of species. Note that we are not working
     with total energy as a conserved quantity, so that we will
     have to carefully account for the effect of the transport 
     processes on potential temperature to ensure total energy 
     conservation.
     */
    extern User_Data ud;
    extern Thermodynamic th;
    
    if (ud.mol_trans == FULL_MOLECULAR_TRANSPORT) {
        assert(0); /* full transport branch of code not fully implemented yet */

        const double twothird  = 2.0/3.0;
        
        for (int i=0; i<elem->nc; i++) {
            diss[i] = 0.0;
        }
        
        /* viscous forces */
        /* x-direction */
        switch (elem->ndim) {
            case 1:
                assert(0); /* not implemented for 1D yet */
                break;
                
            case 2: {
                const int icx = elem->icx;
                const int icy = elem->icy;
                const int ifx = elem->ifx;
                const int ify = elem->ify;
                const int igx = elem->igx;
                const int igy = elem->igy;
                
                const double dx = elem->dx;
                const double dy = elem->dy;
                
                /* fluxes in the x-direction */
                for (int j=igy; j<icy-igy; j++) {
                    int ncj = j*icx;
                    int nfj = j*ifx;
                    for (int i=igx; i<ifx-igx; i++) {
                        int ncij = ncj + i;
                        int nfij = nfj + i;
                        
                        double rhocm   = Sol->rho[ncij-1];
                        double pcm     = pow(Sol->rhoY[ncij-1],th.gamm);
                        double Tcm     = pcm/rhocm;
                        double ucm     = Sol->rhou[ncij-1]/rhocm;
                        double vcm     = Sol->rhov[ncij-1]/rhocm;
                        
                        double rhocp   = Sol->rho[ncij];
                        double pcp     = pow(Sol->rhoY[ncij],th.gamm);
                        double Tcp     = pcp/rhocp;
                        double ucp     = Sol->rhou[ncij]/rhocp;
                        double vcp     = Sol->rhov[ncij]/rhocp;
                        
                        double rhotm   = Sol->rho[ncij+icx-1];
                        double utm     = Sol->rhou[ncij+icx-1]/rhotm;
                        double vtm     = Sol->rhov[ncij+icx-1]/rhotm;
                        
                        double rhotp   = Sol->rho[ncij+icx];
                        double utp     = Sol->rhou[ncij+icx]/rhotp;
                        double vtp     = Sol->rhov[ncij+icx]/rhotp;
                        
                        double rhobm   = Sol->rho[ncij-icx-1];
                        double ubm     = Sol->rhou[ncij-icx-1]/rhobm;
                        double vbm     = Sol->rhov[ncij-icx-1]/rhobm;
                        
                        double rhobp   = Sol->rho[ncij-icx];
                        double ubp     = Sol->rhou[ncij-icx]/rhobp;
                        double vbp     = Sol->rhov[ncij-icx]/rhobp;
                        
                        double visc  = ud.viscm  + ud.visct;
                        double viscb = ud.viscbm + ud.viscbt; 
                        double cond  = ud.cond;
                        
                        double rhobar = 0.5*(rhocp+rhocm);
                        double dudx   = (ucp-ucm)/dx;
                        double dvdx   = (vcp-vcm)/dx;
                        double dTdx   = (Tcp-Tcm)/dx;                    
                        double dudy   = 0.25*(utp-ubp + utm-ubm)/dy;
                        double dvdy   = 0.25*(vtp-vbp + vtm-vbm)/dy;
                        
                        /* fluxes in the x-direction */
                        flux[0]->rhou[nfij] = - rhobar*(2.0*visc*dudx - twothird*viscb*(dudx+dvdy));
                        flux[0]->rhov[nfij] = - rhobar*visc*(dudy + dvdx);
                        flux[0]->rhoe[nfij] = - rhobar*cond*dTdx;
                    }
                }
                
                /* fluxes in the y-direction */
                for (int j=igy; j<icy-igy; j++) {
                    int ncj = j*icx;
                    int nfj = j;
                    for (int i=igx; i<ifx-igx; i++) {
                        int ncij = ncj + i;
                        int nfij = nfj + i*ify;
                        
                        double rhocm   = Sol->rho[ncij-icx];
                        double pcm     = pow(Sol->rhoY[ncij-icx],th.gamm);
                        double Tcm     = pcm/rhocm;
                        double ucm     = Sol->rhou[ncij-icx]/rhocm;
                        double vcm     = Sol->rhov[ncij-icx]/rhocm;
                        
                        double rhocp   = Sol->rho[ncij];
                        double pcp     = pow(Sol->rhoY[ncij],th.gamm);
                        double Tcp     = pcp/rhocp;
                        double ucp     = Sol->rhou[ncij]/rhocp;
                        double vcp     = Sol->rhov[ncij]/rhocp;
                        
                        double rhorm   = Sol->rho[ncij+1-icx];
                        double urm     = Sol->rhou[ncij+1-icx]/rhorm;
                        double vrm     = Sol->rhov[ncij+1-icx]/rhorm;
                        
                        double rhorp   = Sol->rho[ncij+1];
                        double urp     = Sol->rhou[ncij+1]/rhorp;
                        double vrp     = Sol->rhov[ncij+1]/rhorp;
                        
                        double rholm   = Sol->rho[ncij-1-icx];
                        double ulm     = Sol->rhou[ncij-1-icx]/rholm;
                        double vlm     = Sol->rhov[ncij-1-icx]/rholm;
                        
                        double rholp   = Sol->rho[ncij-1];
                        double ulp     = Sol->rhou[ncij-1]/rholp;
                        double vlp     = Sol->rhov[ncij-1]/rholp;
                        
                        double visc  = ud.viscm  + ud.visct;
                        double viscb = ud.viscbm + ud.viscbt; 
                        double cond  = ud.cond;
                        
                        double rhobar = 0.5*(rhocp+rhocm);
                        double dudy   = (ucp-ucm)/dy;
                        double dvdy   = (vcp-vcm)/dy;
                        double dTdy   = (Tcp-Tcm)/dy;                    
                        double dudx   = 0.25*(urp-ulp + urm-ulm)/dx;
                        double dvdx   = 0.25*(vrp-vlp + vrm-vlm)/dx;
                        
                        /* fluxes in the x-direction */
                        flux[1]->rhou[nfij] = - rhobar*visc*(dudy + dvdx);  
                        flux[1]->rhov[nfij] = - rhobar*(2.0*visc*dvdy - twothird*viscb*(dudx+dvdy));
                        flux[1]->rhoe[nfij] = - rhobar*cond*dTdy;
                    }
                }
            }
                break;
                
            case 3:
                assert(0); /* not implemented for 1D yet */
                break;
                
            default:
                break;
        }
    } else if (ud.mol_trans == STRAKA_DIFFUSION_MODEL) {
        /* this implements diffusive fluxes for the simplified diffusion model
         introduced by Straka et al., Int. J. Num. Meth. Fluids, 17, 1–22 (1993) 
         which amounts to constant coefficient diffusion of the velocity 
         components and the potential temperature.
         */
        
        double mu = ud.visct;
        
        for (int i=0; i<elem->nc; i++) {
            diss[i] = 0.0;
        }
        
        switch (elem->ndim) {
            case 1:
                assert(0); /* not implemented for 1D yet */
                break;
                
            case 2: {
                const int icx = elem->icx;
                const int icy = elem->icy;
                const int ifx = elem->ifx;
                const int ify = elem->ify;
                const int igx = elem->igx;
                const int igy = elem->igy;
                
                const double dx = elem->dx;
                const double dy = elem->dy;
                
                /* fluxes in the x-direction */
                for (int j=igy; j<icy-igy; j++) {
                    int ncj = j*icx;
                    int nfj = j*ifx;
                    for (int i=igx; i<ifx-igx; i++) {
                        int ncij = ncj + i;
                        int nfij = nfj + i;
                        
                        double rhocm  = Sol->rho[ncij-1];
                        double Ycm    = Sol->rhoY[ncij-1]/rhocm;
                        double ucm    = Sol->rhou[ncij-1]/rhocm;
                        double vcm    = Sol->rhov[ncij-1]/rhocm;
                        
                        double rhocp  = Sol->rho[ncij];
                        double Ycp    = Sol->rhoY[ncij]/rhocp;
                        double ucp    = Sol->rhou[ncij]/rhocp;
                        double vcp    = Sol->rhov[ncij]/rhocp;
                                                
                        /* fluxes in the x-direction */
                        flux[0]->rhoY[nfij] = - mu*(Ycp-Ycm)/dx;
                        flux[0]->rhou[nfij] = - mu*(ucp-ucm)/dx;
                        flux[0]->rhov[nfij] = - mu*(vcp-vcm)/dx;
                    }
                }
                
                /* fluxes in the y-direction */
                for (int j=igy; j<icy-igy; j++) {
                    int ncj = j*icx;
                    int nfj = j;
                    for (int i=igx; i<ifx-igx; i++) {
                        int ncij = ncj + i;
                        int nfij = nfj + i*ify;
                        
                        double rhocm   = Sol->rho[ncij-icx];
                        double Ycm     = Sol->rhoY[ncij-icx]/rhocm;
                        double ucm     = Sol->rhou[ncij-icx]/rhocm;
                        double vcm     = Sol->rhov[ncij-icx]/rhocm;
                        
                        double rhocp   = Sol->rho[ncij];
                        double Ycp     = Sol->rhoY[ncij]/rhocp;
                        double ucp     = Sol->rhou[ncij]/rhocp;
                        double vcp     = Sol->rhov[ncij]/rhocp;
                                                
                        /* fluxes in the x-direction */
                        flux[1]->rhoY[nfij] = - mu*(Ycp-Ycm)/dy;
                        flux[1]->rhou[nfij] = - mu*(ucp-ucm)/dy;  
                        flux[1]->rhov[nfij] = - mu*(vcp-vcm)/dy;
                    }
                }
            }
                break;
                
            case 3:
                assert(0); /* not implemented for 1D yet */
                break;
                
            default:
                break;
                
        } 
    }
}
