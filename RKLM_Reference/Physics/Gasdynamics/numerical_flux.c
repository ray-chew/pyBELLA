/*******************************************************************************
 File:   numerical_flux.c
 Author: Thomas (Nicola)
 Date:   Wed Feb 25 13:21:46 WET 1998
 *******************************************************************************/
#include <float.h>
#include <math.h>
#include "Common.h"
#include "ProjectionType.h"
#include "math_own.h"
#include "Eos.h"   
#include "thermodynamic.h"
#include "variable.h"
#include "userdata.h"
#include "mpv.h"
#include "error.h"
#include "numerical_flux.h"
#include "boundary.h"


void hllestar(
              ConsVars* Fluxes, 
              States* Lefts, 
              States* Rights,
              States* Sol,
              const double lambda,
              const int n) {
    
    extern User_Data ud;
    
    double rhol, ul, vl, wl, pl, rhoul, Hl, Yl, Zl;
    double rhor, ur, vr, wr, pr, rhour, Hr, Yr, Zr;
    double Xl[NSPEC], Xr[NSPEC];
    double upwind, upl, upr;
    int i, nsp;
    
    primitives(Lefts,  1, n - 1);
    primitives(Rights, 1, n - 1);
    
    /* sanity values for fluxes at "0", "n-1" and "n-2" */
    Fluxes->rho[0]  = Fluxes->rho[n-2]  = Fluxes->rho[n-1]  = 0.0;
    Fluxes->rhou[0] = Fluxes->rhou[n-2] = Fluxes->rhou[n-1] = 0.0;
    Fluxes->rhov[0] = Fluxes->rhov[n-2] = Fluxes->rhov[n-1] = 0.0;
    Fluxes->rhow[0] = Fluxes->rhow[n-2] = Fluxes->rhow[n-1] = 0.0;
    Fluxes->rhoe[0] = Fluxes->rhoe[n-2] = Fluxes->rhoe[n-1] = 0.0;
    Fluxes->rhoZ[0] = Fluxes->rhoZ[n-2] = Fluxes->rhoZ[n-1] = 0.0;

    for (nsp = 0; nsp < ud.nspec; nsp++) {
        Fluxes->rhoX[nsp][0] = Fluxes->rhoX[nsp][n-2] = Fluxes->rhoX[nsp][n-1] = 0.0;
    }

#ifndef NODAL_PROJECTION_ONLY
    Fluxes->rhoY[0] = Fluxes->rhoY[n-2] = Fluxes->rhoY[n-1] = 0.0;
#endif

    for(i = 1; i < n - 2; i++) {	
        rhol  = Lefts->rho[i];
        ul    = Lefts->u[i];
        vl    = Lefts->v[i];
        wl    = Lefts->w[i];
        pl    = Lefts->p[i];
        Yl    = Lefts->Y[i];
        Zl    = Lefts->Z[i];
        
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
        Zr    = Rights->Z[i+1];
        
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Xr[nsp]    = Rights->X[nsp][i+1];
        }
        
        rhour = Rights->rhou[i+1];
        Hr    = Rights->rhoe[i+1] + pr;
        
#ifdef NODAL_PROJECTION_ONLY
        {
            double rhoYu_l  = Sol->rhoY[i]*Sol->rhou[i]/Sol->rho[i];
            double rhoYu_lp = 0.5*(Lefts->rhoY[i]*Lefts->rhou[i]/Lefts->rho[i] + Rights->rhoY[i]*Rights->rhou[i]/Rights->rho[i]);
             
            double rhoYu_r  = Sol->rhoY[i+1]*Sol->rhou[i+1]/Sol->rho[i+1];
            double rhoYu_rp = 0.5*(Lefts->rhoY[i+1]*Lefts->rhou[i+1]/Lefts->rho[i+1] + Rights->rhoY[i+1]*Rights->rhou[i+1]/Rights->rho[i+1]);
            
            Fluxes->rhoY[i] += SECOND_ORDER_ADVECTING_FLUX * 0.5 * ((rhoYu_lp - rhoYu_l) + (rhoYu_rp - rhoYu_r));
        }
#else /* NODAL_PROJECTION_ONLY */
        Fluxes->rhoY[i] = 0.25 * (rhol*Yl+rhor*Yr)*(ul + ur);
#endif /* NODAL_PROJECTION_ONLY */
        
        upwind = 0.5 * ( 1.0 + SIGN(Fluxes->rhoY[i]));
        
        upl    = upwind / Yl;
        upr    = (1.0 - upwind) / Yr;
        Fluxes->rhou[i] = Fluxes->rhoY[i] * (upl * ul  + upr * ur) ;
                
        Fluxes->rho[i]  = Fluxes->rhoY[i] * (upl * 1.0 + upr * 1.0);
        Fluxes->rhoe[i] = Fluxes->rhoY[i] * (upl * Hl  + upr * Hr) ;
        Fluxes->rhov[i] = Fluxes->rhoY[i] * (upl * vl  + upr * vr) ;
        Fluxes->rhow[i] = Fluxes->rhoY[i] * (upl * wl  + upr * wr) ;
        Fluxes->rhoZ[i] = Fluxes->rhoY[i] * (upl * Zl  + upr * Zr) ;
        
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Fluxes->rhoX[nsp][i] = Fluxes->rhoY[i] * (upl * Xl[nsp]  + upr * Xr[nsp]) ;
        }        
    }
}


#ifdef NODAL_PROJECTION_ONLY
/* -------------------------------------------------------------------------- */

void Advective_Fluxes(VectorField* adv_flux, 
                      const ConsVars* Sol, 
                      const ElemSpaceDiscr* elem)
{
    
    /* 
     advective fluxes obtained from cell-centered data by the averaging
     procedure that makes the associated primary cell divergence
     compatible with the dual cell divergence used in the second projection
     */
    extern User_Data ud;
    
    int off_bcx_min[3] = {0,0,0}; 
    int off_bcx_max[3] = {0,0,0}; 
    int off_bcy_min[3] = {0,0,0}; 
    int off_bcy_max[3] = {0,0,0}; 
    
    off_bcx_min[0] = (ud.bdrytype_min[0] == PERIODIC ?   elem->icx-2*elem->igx  :  1);
    off_bcx_max[2] = (ud.bdrytype_max[0] == PERIODIC ? -(elem->icx-2*elem->igx) : -1);
    off_bcy_min[0] = (ud.bdrytype_min[1] == PERIODIC ?   elem->icy-2*elem->igy  :  1);
    off_bcy_max[2] = (ud.bdrytype_max[1] == PERIODIC ? -(elem->icy-2*elem->igy) : -1);
    
    switch (elem->ndim) {
        case 1:
            {                                
                const int ifx = elem->ifx;
                const int igx = elem->igx;
                
                double* prho  = Sol->rho;
                double* prhoY = Sol->rhoY;
                double* prhou = Sol->rhou;
                                
                /* x - advective fluxes */
                for (int nf=0; nf<elem->nfx; nf++) {
                    adv_flux->x[nf] = 0.0;
                }
                
                for (int i=igx-1; i<ifx-igx+1; i++) {
                    int noffc = i;
                    int noffm = i-1;
                    int ni_f  = i;
                    
                    double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                    double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                        
                    adv_flux->x[ni_f] += 0.5 * (mYc+mYm);
                }
            }
            
            Bound_adv_flux(adv_flux, elem);
            break;

        case 2:
            {
                const int icx = elem->icx;
                const int icy = elem->icy;
                                
                const int ifx = elem->ifx;
                const int ify = elem->ify;
                
                const int igx = elem->igx;
                const int igy = elem->igy;
                
                const int str[2] = {elem->stride[0],elem->stride[1]};

                double* prho  = Sol->rho;
                double* prhoY = Sol->rhoY;
                double* prhou = Sol->rhou;
                double* prhov = Sol->rhov;

                int weight[3] = {1,2,1};

                /* x - advective fluxes */
                for (int nf=0; nf<elem->nfx; nf++) {
                    adv_flux->x[nf] = 0.0;
                }
                
                /* bottom boundary row */
                for (int i=igx; i<ifx-igx; i++) {
                    int nji_c = igy*icx + i;
                    int nji_f = igy*ifx + i;
                    double wesum = 0.0;
                    
                    for (int jj=-1; jj<2; jj++) {
                        int jjj    = jj + off_bcy_min[jj+1];
                        int noffc  = nji_c + jjj*icx;
                        int noffm  = nji_c + jjj*icx - 1;
                        double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                        double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                        double we  = weight[jj+1];
                        
                        adv_flux->x[nji_f] += we * 0.5 * (mYc+mYm);
                        wesum += we;
                    }
                    adv_flux->x[nji_f] /= wesum;
                }

                /* core of the domain */
                for (int j=igy+1; j<icy-igy-1; j++) {
                    int nj_c = j*icx;
                    int nj_f = j*ifx;
                                        
                    for (int i=igx; i<ifx-igx; i++) {
                        int nji_c = nj_c + i;
                        int nji_f = nj_f + i;
                        double wesum = 0.0;
                        
                        for (int jj=-1; jj<2; jj++) {
                            int noffc  = nji_c + jj*str[1];
                            int noffm  = nji_c + jj*str[1] - str[0];
                            double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                            double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                            double we  = weight[jj+1];
                            
                            adv_flux->x[nji_f] += we * 0.5 * (mYc+mYm);
                            wesum += we;
                        }
                        adv_flux->x[nji_f] /= wesum;
                    }
                }
                
                /* top boundary row */
                for (int i=igx; i<ifx-igx; i++) {
                    int nji_c = (icy-igy-1)*icx + i;
                    int nji_f = (icy-igy-1)*ifx + i;
                    double wesum = 0.0;
                    
                    for (int jj=-1; jj<2; jj++) {
                        int jjj    = jj + off_bcy_max[jj+1];
                        int noffc  = nji_c + jjj*icx;
                        int noffm  = nji_c + jjj*icx - 1;
                        double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                        double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                        double we  = weight[jj+1];
                        
                        adv_flux->x[nji_f] += we * 0.5 * (mYc+mYm);
                        wesum += we;
                    }
                    adv_flux->x[nji_f] /= wesum;
                }

                
                /* y - advective fluxes */
                for (int nf=0; nf<elem->nfy; nf++) {
                    adv_flux->y[nf] = 0.0;
                }
                                                
                /* left boundary column */
                for (int j=igy; j<ify-igy; j++) {
                    int nj_c = j*icx;
                    int nj_f = j;
                    int i    = igx; 
                    int nji_c = nj_c + i;
                    int nji_f = nj_f + i*ify;
                    double wesum = 0.0;
                    
                    for (int ii=-1; ii<2; ii++) {
                        int iii    = ii + off_bcx_min[ii+1];
                        int noffc  = nji_c + iii;
                        int noffm  = nji_c + iii - icx;
                        double mYc = prhoY[noffc]*prhov[noffc]/prho[noffc];
                        double mYm = prhoY[noffm]*prhov[noffm]/prho[noffm];
                        double we  = weight[ii+1];
                        
                        adv_flux->y[nji_f] += we * 0.5 * (mYc+mYm);
                        wesum += we;
                    }
                    adv_flux->y[nji_f] /= wesum;
                }

                /* core of the domain */
                for (int j=igy; j<ify-igy; j++) {
                    int nj_c = j*icx;
                    int nj_f = j;
                    
                    for (int i=igx+1; i<icx-igx-1; i++) {
                        int nji_c = nj_c + i;
                        int nji_f = nj_f + i*ify;
                        double wesum = 0.0;
                                                
                        for (int ii=-1; ii<2; ii++) {
                            int noffc  = nji_c + ii;
                            int noffm  = nji_c + ii - icx;
                            double mYc = prhoY[noffc]*prhov[noffc]/prho[noffc];
                            double mYm = prhoY[noffm]*prhov[noffm]/prho[noffm];
                            double we  = weight[ii+1];
                            
                            adv_flux->y[nji_f] += we * 0.5 * (mYc+mYm);
                            wesum += we;
                        }
                        adv_flux->y[nji_f] /= wesum;
                    }
                }
                
                /* right boundary column */
                for (int j=igy; j<ify-igy; j++) {
                    int nj_c = j*icx;
                    int nj_f = j;
                    int i    = icx-igx-1; 
                    int nji_c = nj_c + i;
                    int nji_f = nj_f + i*ify;
                    double wesum = 0.0;
                    
                    for (int ii=-1; ii<2; ii++) {
                        int iii    = ii + off_bcx_max[ii+1];
                        int noffc  = nji_c + iii;
                        int noffm  = nji_c + iii - icx;
                        double mYc = prhoY[noffc]*prhov[noffc]/prho[noffc];
                        double mYm = prhoY[noffm]*prhov[noffm]/prho[noffm];
                        double we  = weight[ii+1];
                        
                        adv_flux->y[nji_f] += we * 0.5 * (mYc+mYm);
                        wesum += we;
                    }
                    adv_flux->y[nji_f] /= wesum;
                }

            }
            
            Bound_adv_flux(adv_flux, elem);
            break;

        case 3:
            {
                ERROR("Advective_Fluxes() implementation for 3D not checked yet");

                /*
                const int icx = elem->icx;
                const int icy = elem->icy;
                const int icz = elem->icz;
                
                const int str[3] = {elem->stride[0],elem->stride[1],elem->stride[2]};
                
                const int ifx = elem->ifx;
                const int ify = elem->ify;
                const int ifz = elem->ifz;
                
                const int igx = elem->igx;
                const int igy = elem->igy;
                const int igz = elem->igz;
                
                double* prho  = Sol->rho;
                double* prhoY = Sol->rhoY;
                double* prhou = Sol->rhou;
                double* prhov = Sol->rhov;
                double* prhow = Sol->rhow;
                
                int weight[3] = {1,2,1};
                  */
                
                /* x - advective fluxes 
                for (int nf=0; nf<elem->nfx; nf++) {
                    adv_flux->x[nf] = 0.0;
                }
                
                for (int k=igz-1; k<icz-igz+1; k++) {
                    int nk_c = k*icy*icx;
                    int nk_f = k*icy*ifx;
                    
                    for (int j=igy-1; j<icy-igy+1; j++) {
                        int njk_c = nk_c + j*icx;
                        int njk_f = nk_f + j*ifx;
                        
                        for (int i=igx-1; i<ifx-igx+1; i++) {
                            int nkji_c = njk_c + i;
                            int nkji_f = njk_f + i;
                            double wesum = 0.0;
                            
                            for (int kk=-1; kk<2; kk++) {
                                for (int jj=-1; jj<2; jj++) {
                                    int noffc  = nkji_c + kk*str[2] + jj*str[1];
                                    int noffm  = nkji_c + kk*str[2] + jj*str[1] - str[0];
                                    double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                                    double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                                    double we  = weight[kk+1] * weight[jj+1];
                                    
                                    adv_flux->x[nkji_f] += we * 0.5 * (mYc+mYm);
                                    wesum += we;
                                }
                            }
                            adv_flux->x[nkji_f] /= wesum;
                        }
                    }
                }
                */
                
                /* y - advective fluxes 
                for (int nf=0; nf<elem->nfy; nf++) {
                    adv_flux->y[nf] = 0.0;
                }
                
                for (int k=igz-1; k<icz-igz+1; k++) {
                    int nk_c = k*icy*icx;
                    int nk_f = k*ify;
                    
                    for (int j=igy-1; j<ify-igy+1; j++) {
                        int njk_c = nk_c + j*icx;
                        int njk_f = nk_f + j;
                        
                        for (int i=igx-1; i<icx-igx+1; i++) {
                            int nkji_c = njk_c + i;
                            int nkji_f = njk_f + i*icz*ify;
                            double wesum = 0.0;
                            
                            for (int ii=-1; ii<2; ii++) {
                                for (int kk=-1; kk<2; kk++) {
                                    int noffc  = nkji_c + ii*str[0] + kk*str[2];
                                    int noffm  = nkji_c + ii*str[0] + kk*str[2] - str[1];
                                    double mYc = prhoY[noffc]*prhov[noffc]/prho[noffc];
                                    double mYm = prhoY[noffm]*prhov[noffm]/prho[noffm];
                                    double we  = weight[ii+1] * weight[kk+1];

                                    adv_flux->y[nkji_f] += we * 0.5 * (mYc+mYm);
                                    wesum += we;
                                }
                            }
                            adv_flux->y[nkji_f] /= wesum;
                        }
                    }
                }
                */
                
                /* z - advective fluxes 
                for (int nf=0; nf<elem->nfz; nf++) {
                    adv_flux->z[nf] = 0.0;
                }
                
                for (int k=igz-1; k<icz-igz+1; k++) {
                    int nk_c = k*icy*icx;
                    int nk_f = k;
                    
                    for (int j=igy-1; j<icy-igy+1; j++) {
                        int njk_c = nk_c + j*icx;
                        int njk_f = nk_f + j*icx*ifz;
                        
                        for (int i=igx-1; i<icx-igx+1; i++) {
                            int nkji_c = njk_c + i;
                            int nkji_f = njk_f + i*ifz;
                            double wesum = 0.0;
                            
                            for (int jj=-1; jj<2; jj++) {
                                for (int ii=-1; ii<2; ii++) {
                                    int noffc  = nkji_c + jj*str[1] + ii*str[0];
                                    int noffm  = nkji_c + jj*str[1] + ii*str[0] - str[2];
                                    double mYc = prhoY[noffc]*prhow[noffc]/prho[noffc];
                                    double mYm = prhoY[noffm]*prhow[noffm]/prho[noffm];
                                    double we  = weight[jj+1] * weight[ii+1];

                                    adv_flux->z[nkji_f] += we * 0.5 * (mYc+mYm);
                                    wesum += we;
                                }
                            }
                            adv_flux->z[nkji_f] /= wesum;
                        }
                    }
                }
                
            }  
            
            Bound_adv_flux(adv_flux, elem);
                 */
            break;

        default:
            break;
        }
    }
}

/* -------------------------------------------------------------------------- */

void Advective_Fluxes_x(double* rhoYu, 
                        const ConsVars* Sol, 
                        const ElemSpaceDiscr* elem, 
                        const int SplitStep)
{
    /* 
     advective fluxes obtained from cell-centered data by the averaging
     procedure that makes the associated primary cell divergence
     compatible with the dual cell divergence used in the second projection

     Here: only for current x-direction in OpSplit mode.
     */
    
    extern User_Data ud;
    
    switch (elem->ndim) {
        case 1:
        {                                
            ERROR("Advective_Fluxes_x() implementation not checked for 1D yet");

            /*
            const int ifx = elem->ifx;
            const int igx = elem->igx;
            
            double* prho  = Sol->rho;
            double* prhoY = Sol->rhoY;
            double* prhou = Sol->rhou;
             */
            
            /* x - advective fluxes 
            for (int nf=0; nf<elem->nfx; nf++) {
                rhoYu[nf] = 0.0;
            }
            
            for (int i=igx-1; i<ifx-igx+1; i++) {
                int noffc = i;
                int noffm = i-1;
                int ni_f  = i;
                
                double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                
                rhoYu[ni_f] += 0.5 * (mYc+mYm);
            }
             */
        }
            break;
            
        case 2:
        {
            const int icx = elem->icx;
            const int icy = elem->icy;
            
            const int ifx = elem->ifx;
            
            const int igx = elem->igx;
            const int igy = elem->igy;
                        
            double* prho  = Sol->rho;
            double* prhoY = Sol->rhoY;
            double* prhou = Sol->rhou;

            int weight[3] = {1,2,1};

            int offset_bcy_min[3] = {0,0,0}; 
            int offset_bcy_max[3] = {0,0,0}; 
            
            offset_bcy_min[0] = (ud.bdrytype_min[1-SplitStep] == PERIODIC ?   elem->icy-2*elem->igy  :  1);
            offset_bcy_max[2] = (ud.bdrytype_max[1-SplitStep] == PERIODIC ? -(elem->icy-2*elem->igy) : -1);
            
            /* x - advective fluxes */
            for (int nf=0; nf<elem->nfx; nf++) {
                rhoYu[nf] = 0.0;
            }
            
            /* lower boundary row 
            for (int i=igx-1; i<ifx-igx+1; i++) {
                int nji_c = igy*icx + i;
                int nji_f = igy*ifx + i;
                double wesum = 0.0;
                
                for (int jj=-1; jj<2; jj++) {
                    int jjj    = jj + offset_bcy_min[jj+1];
                    int noffc  = nji_c + jjj*icx;
                    int noffm  = nji_c + jjj*icx - 1;
                    double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                    double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                    double we  = weight[jj+1];
                    
                    rhoYu[nji_f] += we * 0.5 * (mYc+mYm);
                    wesum += we;
                }
                rhoYu[nji_f] /= wesum;
            }
             */
            
            /* core of the domain */
            for (int j=igy; j<icy-igy; j++) {
                int nj_c = j*icx;
                int nj_f = j*ifx;

                /* symmetry test variables 
                int nj_c_mirr = (icy-1-j)*icx;
                int nj_f_mirr = (icy-1-j)*ifx;
                 */
                
                for (int i=igx-1; i<ifx-igx+1; i++) {
                    int nji_c = nj_c + i;
                    int nji_f = nj_f + i;
                    
                    /* symmetry test variables 
                    int nji_c_mirr = nj_c_mirr + i;
                    int nji_f_mirr = nj_f_mirr + i;
                     */
                    double wesum = 0.0;
                    double rhoYu_mirr = 0.0;
                    
                    for (int jj=-1; jj<2; jj++) {
                        double we  = weight[jj+1];
                        int noffc  = nji_c + jj*icx;
                        int noffm  = nji_c + jj*icx - 1;
                        double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                        double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                        
                        /* test variables 
                        int noffc_mirr  = nji_c_mirr - jj*icx;
                        int noffm_mirr  = nji_c_mirr - jj*icx - 1;
                        double mYc_mirr = prhoY[noffc_mirr]*prhou[noffc_mirr]/prho[noffc_mirr];
                        double mYm_mirr = prhoY[noffm_mirr]*prhou[noffm_mirr]/prho[noffm_mirr];
                         rhoYu_mirr   += we * 0.5 * (mYc_mirr+mYm_mirr);
                          */
                        rhoYu[nji_f] += we * 0.5 * (mYc+mYm);
                        wesum += we;
                    }
                    rhoYu[nji_f] /= wesum;
                    /*
                    rhoYu_mirr /= wesum;
                     */
                    wesum = 0.0;
                }
            }

            /* top boundary row 
            for (int i=igx-1; i<ifx-igx+1; i++) {
                int nji_c = (icy-igy-1)*icx + i;
                int nji_f = (icy-igy-1)*ifx + i;
                double wesum = 0.0;
                
                for (int jj=-1; jj<2; jj++) {
                    int jjj    = jj + offset_bcy_max[jj+1];
                    int noffc  = nji_c + jjj*icx;
                    int noffm  = nji_c + jjj*icx - 1;
                    double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                    double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                    double we  = weight[jj+1];
                    
                    rhoYu[nji_f] += we * 0.5 * (mYc+mYm);
                    wesum += we;
                }
                rhoYu[nji_f] /= wesum;
            }
            */
            
            /* call boundary routine */
            Bound_adv_flux_x(rhoYu, elem, SplitStep);
            /* printf("I am here to stop\n"); */
        }
            break;
            
        case 3:
        {
            ERROR("Advective_Fluxes_x() implementation not checked for 3D yet");

            /*
            const int icx = elem->icx;
            const int icy = elem->icy;
            const int icz = elem->icz;
            
            const int str[3] = {elem->stride[0],elem->stride[1],elem->stride[2]};
            
            const int ifx = elem->ifx;

            const int igx = elem->igx;
            const int igy = elem->igy;
            const int igz = elem->igz;
            
            double* prho  = Sol->rho;
            double* prhoY = Sol->rhoY;
            double* prhou = Sol->rhou;
            
            int weight[3] = {1,2,1};
             */
            
            /* x - advective fluxes 
            for (int nf=0; nf<elem->nfx; nf++) {
                rhoYu[nf] = 0.0;
            }
            
            for (int k=igz-1; k<icz-igz+1; k++) {
                int nk_c = k*icy*icx;
                int nk_f = k*icy*ifx;
                
                for (int j=igy-1; j<icy-igy+1; j++) {
                    int njk_c = nk_c + j*icx;
                    int njk_f = nk_f + j*ifx;
                    
                    for (int i=igx-1; i<ifx-igx+1; i++) {
                        int nkji_c = njk_c + i;
                        int nkji_f = njk_f + i;
                        double wesum = 0.0;
                        
                        for (int kk=-1; kk<2; kk++) {
                            for (int jj=-1; jj<2; jj++) {
                                int noffc  = nkji_c + kk*str[2] + jj*str[1];
                                int noffm  = nkji_c + kk*str[2] + jj*str[1] - str[0];
                                double mYc = prhoY[noffc]*prhou[noffc]/prho[noffc];
                                double mYm = prhoY[noffm]*prhou[noffm]/prho[noffm];
                                double we  = weight[kk+1] * weight[jj+1];
                                
                                rhoYu[nkji_f] += we * 0.5 * (mYc+mYm);
                                wesum += we;
                            }
                        }
                        rhoYu[nkji_f] /= wesum;
                    }
                }
            }
             */
        }            
            break;
            
        default:
            break;
    }
    
}
#endif /* NODAL_PROJECTION_ONLY */


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: numerical_flux.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:35  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
