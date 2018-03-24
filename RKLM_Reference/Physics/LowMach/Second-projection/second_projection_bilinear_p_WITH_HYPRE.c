
/*******************************************************************************
 File:   second_projection.c
 Author: Nicola
 Date:   Fri Mar 13 14:58:07 WET 1998
 *******************************************************************************/
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <tgmath.h>
#include "Common.h"
#include "BiCGSTAB.h"
#include "variable.h"
#include "mpv.h"
#include "error.h"
#include "warning.h"
#include "userdata.h"
#include "variable_coefficient_poisson_nodes.h"
#include "second_projection_bilinear_p.h" 
#include "io.h"
#include "userdata.h"
#include "enumerator.h"
#include "Eos.h"
#include "thermodynamic.h"
#include "boundary.h"
#include "limiter.h"
#include "math_own.h"
#include "memory.h"
#include "laplacian_nodes.h"
#include "numerical_flux.h"
#include "flux_correction.h"

static double divergence_nodes(
							 double* rhs,
							 const ElemSpaceDiscr* elem,
							 const NodeSpaceDiscr* node,
							 const ConsVars* Sol,
							 const MPV* mpv,
							 const BDRY* bdry,
							 const double dt,
							 const double weight);


static enum Constraint integral_condition_nodes(
												double* rhs,
												const NodeSpaceDiscr* node,
												const int x_periodic,
												const int y_periodic,
												const int z_periodic);

static void	catch_periodic_directions(
									  double* rhs,  
									  const NodeSpaceDiscr* node, 
									  const ElemSpaceDiscr* elem,
									  const int x_periodic,
									  const int y_periodic,
									  const int z_periodic);

static void operator_coefficients_nodes(
										double* Splus[3], 
										double* wcenter,
										const ElemSpaceDiscr* elem,
										const NodeSpaceDiscr* node,
										ConsVars* Sol,
										const ConsVars* Sol0,
										const MPV* mpv,
                                        const double dt);

void correction_nodes(
					  ConsVars* Sol,
					  const ElemSpaceDiscr* elem,
					  const NodeSpaceDiscr* node,
                      const double* hplus[3], 
					  double* p2, 
					  const double t,
					  const double dt);

/* ========================================================================== */

#define OUTPUT_RHS 0
#if OUTPUT_RHS
static int rhs_output_count = 0;
#endif

void second_projection(
                       ConsVars* Sol,
                       MPV* mpv,
                       const ConsVars* Sol0,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       const double t,
                       const double dt) 
{
	
	extern User_Data ud;
	extern BDRY* bdry;
        
	const int nc = node->nc;
	
	double** hplus   = mpv->Level[0]->wplus;
	double*  hcenter = mpv->Level[0]->wcenter;
	
	double* rhs      = mpv->Level[0]->rhs;
	double* p2       = mpv->Level[0]->p;
	
    double rhs_max;
    
	int x_periodic, y_periodic, z_periodic;
	int ii;
	
    
    printf("\n\n====================================================");
    printf("\nSecond Projection");
    printf("\n====================================================\n");

	/* switch for stopping loops in solver for periodic data before
	 last grid line in each direction */
	x_periodic = 0;
	y_periodic = 0;
	z_periodic = 0;
	if(ud.bdrytype_min[0] == PERIODIC) x_periodic = 1;
	if(ud.bdrytype_min[1] == PERIODIC) y_periodic = 1;
	if(ud.bdrytype_min[2] == PERIODIC) z_periodic = 1;
	
    /* KEEP_OLD_POISSON_SOLUTIONS */
	for(ii=0; ii<nc; ii++){
		p2[ii] = mpv->p2_nodes[ii];
		rhs[ii] = 0.0;
	}
    
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv, bdry, dt, 1.0);
    printf("\nrhsmax = %e\n", rhs_max);

    catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);

#if OUTPUT_RHS
    extern User_Data ud;
    FILE *prhsfile = NULL;
    char fn2[200], fieldname2[90], step_string[30];
    
    if(rhs_output_count<10) {
        sprintf(step_string, "00%d", rhs_output_count);
    }
    else if(rhs_output_count<100) {
        sprintf(step_string, "0%d", rhs_output_count);
    }
    else {
        sprintf(step_string, "%d", rhs_output_count);
    }
    

    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_%s.hdf", ud.file_name, step_string);
    sprintf(fieldname2, "rhs_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->node->icx,
             mpv->Level[0]->node->icy,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->node->ndim,
             rhs,
             fn2,
             fieldname2);
    
    rhs_output_count++;
    
#endif

    assert(integral_condition_nodes(rhs, node, x_periodic, y_periodic, z_periodic) != VIOLATED); 
    
    if (ud.is_compressible) {
        for (int nn=0; nn<node->nc; nn++) {
            rhs[nn] += hcenter[nn]*mpv->p2_nodes[nn];
        }
    }
     
         
#if 0  
    
    if(rhs_output_count<10) {
        sprintf(step_string, "00%d", rhs_output_count);
    }
    else if(rhs_output_count<100) {
        sprintf(step_string, "0%d", rhs_output_count);
    }
    else {
        sprintf(step_string, "%d", rhs_output_count);
    }

    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_%s.hdf", ud.file_name, step_string);
    sprintf(fieldname2, "rhs_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->node->icx,
             mpv->Level[0]->node->icy,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->node->ndim,
             rhs,
             fn2,
             fieldname2);
    
    rhs_output_count++;
    
#endif
    
	operator_coefficients_nodes(hplus, hcenter, elem, node, Sol, Sol0, mpv, dt);
	variable_coefficient_poisson_nodes(p2, (const double **)hplus, hcenter, rhs, x_periodic, y_periodic, z_periodic, dt);
    correction_nodes(Sol, elem, node, (const double**)hplus, p2, t, dt);
    
#if OUTPUT_RHS
    for(ii=0; ii<nc; ii++) rhs[ii] = 0.0;

    if(rhs_output_count<10) {
        sprintf(step_string, "00%d", rhs_output_count);
    }
    else if(rhs_output_count<100) {
        sprintf(step_string, "0%d", rhs_output_count);
    }
    else {
        sprintf(step_string, "%d", rhs_output_count);
    }

    
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv, bdry, dt, 1.0);
    printf("\nrhsmax = %e, test-output\n", rhs_max);
    catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);

    /**/
    for(int k = node->igz; k < node->icz - node->igz; k++) {
        const int ln = k * node->icy * node->icx; 
        for(int j = node->igy; j < node->icy - node->igy; j++) {
            const int mn = ln + j * node->icx; 
            for(int i = node->igx; i < node->icx - node->igx; i++) {
                const int nn = mn + i;
                if (fabs(rhs[nn]) > 1.0e-3) {
                    printf("\nHere we go! (i,j,k) = (%d, %d, %d); rhs = %e\n", i,j,k, rhs[nn]);
                }
            }
        }
    }
    
    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_%s.hdf", ud.file_name, step_string);
    sprintf(fieldname2, "rhs_nodes");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->node->icx,
             mpv->Level[0]->node->icy,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->node->ndim,
             rhs,
             fn2,
             fieldname2);
    
    /*
    EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, mpv->p2_nodes, (const double**)hplus, hcenter, x_periodic, y_periodic, z_periodic, rhs);
    
    sprintf(fn2, "%s/lap_nodes/lap_nodes_%s.hdf", ud.file_name, step_string);
    sprintf(fieldname2, "lap_nodes");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->node->icx,
             mpv->Level[0]->node->icy,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->node->ndim,
             rhs,
             fn2,
             fieldname2);
    
    assert(0);
     */
    /**/
    sprintf(fn2, "%s/p2_nodes/p2_n_%s.hdf", ud.file_name, step_string);
    sprintf(fieldname2, "p2_nodes");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->node->icx,
             mpv->Level[0]->node->icy,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->node->ndim,
             p2,
             fn2,
             fieldname2);
    
    
    rhs_output_count++;
    
#endif

    
    for(ii=0; ii<nc; ii++) {
        mpv->dp2_nodes[ii] = p2[ii] - mpv->p2_nodes[ii];
        mpv->p2_nodes[ii]  = p2[ii];
    }

    Set_Explicit_Boundary_Data(Sol, elem);
    set_ghostnodes_p2(mpv->p2_nodes, node, 2);    
}

/* ========================================================================== */

static double divergence_nodes(
							 double* rhs,
							 const ElemSpaceDiscr* elem,
							 const NodeSpaceDiscr* node,
							 const ConsVars* Sol,
							 const MPV* mpv,
							 const BDRY* bdry,
							 const double dt,
							 const double weight) {
	
    /* with weight = 1.0, this routine computes 
          rhs_out = rhs_in + (2/dt) * div(rhoY\vec{v}) 
       from the current Sol. 
     */

	extern User_Data ud;
	
	const int ndim = node->ndim;

    double div_max = 0.0;
    
    int i, j, k, mn;
    
#ifdef FOURTH_ORDER_ADV_FLUXES
    extern ConsVars* flux[3];
    extern double* W0;
    
    switch(ndim) {
        case 1: {
            ERROR("divergence_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
            
            const int icxn = node->icx;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
                        
            double Y, rhsmax;
            double *rhs_cell = W0;
            
            /* build nodal divergence from cell-centered divergence averaged to nodes */
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
            rhsmax = controlled_variable_flux_divergence(rhs_cell, (const ConsVars**)flux, dt, elem);

            /* predicted time level divergence via scattering */
            for(j = igye; j < icye - igye; j++) {
                const int me = j * icxe;
                const int mn = j * icxn; 
                for(i = igxe; i < icxe - igxe; i++) {
                    const int ne    = me + i; 
                    const int n     = mn + i;
                    const int nicx  = n  + icxn;
                    const int n1    = n  + 1;
                    const int n1icx = n1 + icxn;
                                        
                    rhs[n]     += 0.25*rhs_cell[ne];
                    rhs[n1]    += 0.25*rhs_cell[ne];
                    rhs[n1icx] += 0.25*rhs_cell[ne];
                    rhs[nicx]  += 0.25*rhs_cell[ne];
                }
            } 
            
            
            /* account for influx bottom boundary */
            j = igye;
            mn = j * icxn; 
            Y  = mpv->HydroState->Y0[j];
            for(i = igxe; i < icxe - igxe; i++) {
                const int n     = mn + i;
                const int n1    = n  + 1;
                
                double oowdydt   = (2.0/dt) * weight * 0.25 / elem->dy;
                double rhov_wall = bdry->wall_massflux[i]; 
                double tmpy      = oowdydt * Y * rhov_wall;
                
                rhs[n]  += - tmpy;
                rhs[n1] += - tmpy;
            }
            
            for(j = igye+1; j < icye - igye; j++) {
                const int mn = j * icxn; 
                for(i = igxe+1; i < icxe - igxe; i++) {
                    int nn    = mn + i;
                    div_max = MAX_own(div_max, fabs(rhs[nn]));
                }
            }
            break;
        }
            
        case 3: {
            
            const int icxn = node->icx;
            const int icyn = node->icy;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            const int igze = elem->igz;
            const int icze = elem->icz;
            
            const int dixn = 1;
            const int diyn = icxn;
            const int dizn = icxn*icyn;
            
            double Y, rhsmax;
            double *rhs_cell = W0;
            
            /* build nodal divergence from cell-centered divergence averaged to nodes */
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
            rhsmax = controlled_variable_flux_divergence(rhs_cell, (const ConsVars**)flux, dt, elem);
            
            /* predicted time level divergence via scattering */
            for(k = igze; k < icze - igze; k++) {
                const int le = k * icxe*icye;
                const int ln = k * icxn*icyn; 
                for(j = igye; j < icye - igye; j++) {
                    const int me = le + j * icxe;
                    const int mn = ln + j * icxn; 
                    for(i = igxe; i < icxe - igxe; i++) {
                        const int ne     = me + i; 
                        const int nn000  = mn + i;  /* foresee consistent interpretation of "abc" in nnabc between parts of code */
                        const int nn010  = nn000 + diyn;
                        const int nn011  = nn000 + diyn + dizn;
                        const int nn001  = nn000        + dizn;
                        const int nn100  = nn000 + dixn;
                        const int nn110  = nn000 + dixn + diyn;
                        const int nn111  = nn000 + dixn + diyn + dizn;
                        const int nn101  = nn000 + dixn        + dizn;
                        
                        rhs[nn000] += 0.125*rhs_cell[ne];
                        rhs[nn010] += 0.125*rhs_cell[ne];
                        rhs[nn011] += 0.125*rhs_cell[ne];
                        rhs[nn001] += 0.125*rhs_cell[ne];
                        rhs[nn100] += 0.125*rhs_cell[ne];
                        rhs[nn110] += 0.125*rhs_cell[ne];
                        rhs[nn111] += 0.125*rhs_cell[ne];
                        rhs[nn101] += 0.125*rhs_cell[ne];
                    }
                }
            }
            
            /* account for influx bottom boundary */
            j = igye;
            Y  = mpv->HydroState->Y0[j];
            for(k = igze; k < icze - igze; k++) {
                int ln = k * icyn*icxn;
                int mn = ln + j*icxn;
                for(i = igxe; i < icxe - igxe; i++) {
                    int nn   = mn + i;
                    int nn00 = nn;
                    int nn10 = nn + dixn;
                    int nn11 = nn + dixn + dizn;
                    int nn01 = nn +      + dizn;
                    
                    double oowdydt = (2.0/dt) * weight * 0.25 / elem->dy;
                    double rhov_wall = bdry->wall_massflux[i]; 
                    double tmpy = oowdydt * Y * rhov_wall;
                    
                    rhs[nn00] += - tmpy;
                    rhs[nn10] += - tmpy;
                    rhs[nn01] += - tmpy;
                    rhs[nn11] += - tmpy;
                }
            }        
            
            for(k = igze+1; k < icze - igze; k++) {
                int ln = k*icxn*icyn; 
                for(j = igye+1; j < icye - igye; j++) {
                    int mn = ln + j * icxn; 
                    for(i = igxe+1; i < icxe - igxe; i++) {
                        int nn  = mn + i;
                        div_max = MAX_own(div_max, fabs(rhs[nn]));
                    }
                }
            }
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
#else
    switch(ndim) {
        case 1: {
            ERROR("divergence_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
            
            const int icxn = node->icx;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double dx      = node->dx;
            const double dy      = node->dy;
            const double oodx    = 1.0 / dx;
            const double oody    = 1.0 / dy;
            const double oowdxdt = (2.0/dt) * weight * 0.5 * oodx;
            const double oowdydt = (2.0/dt) * weight * 0.5 * oody;
            
            double Y;
            
            /* predicted time level divergence via scattering */
            for(j = igye; j < icye - igye; j++) {
                const int me = j * icxe;
                const int mn = j * icxn; 
                for(i = igxe; i < icxe - igxe; i++) {
                    const int n     = mn + i;
                    const int nicx  = n  + icxn;
                    const int n1    = n  + 1;
                    const int n1icx = n1 + icxn;
                    
                    const int ne = me + i; 
                    double tmpfx, tmpfy;
                    
                    Y = Sol->rhoY[ne] / Sol->rho[ne];
                    tmpfx = oowdxdt * Y * Sol->rhou[ne];
                    tmpfy = oowdydt * Y * Sol->rhov[ne];
                    
                    rhs[n]           += + tmpfx + tmpfy;
                    rhs[n1]          += - tmpfx + tmpfy;
                    rhs[n1icx]       += - tmpfx - tmpfy;
                    rhs[nicx]        += + tmpfx - tmpfy;
                }
            } 
            
            
            /* account for influx bottom boundary */
            j = igye;
            mn = j * icxn; 
            Y  = mpv->HydroState->Y0[j];
            for(i = igxe; i < icxe - igxe; i++) {
                const int n     = mn + i;
                const int n1    = n  + 1;
                
                double rhov_wall = bdry->wall_massflux[i]; 
                double tmpy = oowdydt * Y * rhov_wall;
                
                rhs[n]  += - tmpy;
                rhs[n1] += - tmpy;
            }
            
            for(j = igye+1; j < icye - igye; j++) {
                const int mn = j * icxn; 
                for(i = igxe+1; i < icxe - igxe; i++) {
                    int nn    = mn + i;
                    div_max = MAX_own(div_max, fabs(rhs[nn]));
                }
            }
            break;
        }
            
        case 3: {
            
            const int icxn = node->icx;
            const int icyn = node->icy;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            const int igze = elem->igz;
            const int icze = elem->icz;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            
            const double oodx = 1.0 / dx;
            const double oody = 1.0 / dy;
            const double oodz = 1.0 / dz;
            
            const double oowdxdt = (2.0/dt) * weight * 0.25 * oodx;
            const double oowdydt = (2.0/dt) * weight * 0.25 * oody;
            const double oowdzdt = (2.0/dt) * weight * 0.25 * oodz;
            
            const int dixn = 1;
            const int diyn = icxn;
            const int dizn = icxn*icyn;
            
            double Y;
            
            /* predicted time level divergence via scattering */
            for(k = igze; k < icze - igze; k++) {
                const int le = k * icye * icxe;
                const int ln = k * icyn * icxn; 
                for(j = igye; j < icye - igye; j++) {
                    const int me = le + j * icxe;
                    const int mn = ln + j * icxn; 
                    for(i = igxe; i < icxe - igxe; i++) {
                        const int ne = me + i; 
                        const int nn000  = mn + i;  /* foresee consistent interpretation of "abc" in nnabc between parts of code */
                        const int nn010  = nn000 + diyn;
                        const int nn011  = nn000 + diyn + dizn;
                        const int nn001  = nn000        + dizn;
                        const int nn100  = nn000 + dixn;
                        const int nn110  = nn000 + dixn + diyn;
                        const int nn111  = nn000 + dixn + diyn + dizn;
                        const int nn101  = nn000 + dixn        + dizn;
                        
                        double tmpfx, tmpfy, tmpfz;
                        
                        Y = Sol->rhoY[ne] / Sol->rho[ne];
                        tmpfx = oowdxdt * Y * Sol->rhou[ne]; /* (rhou*Y) * 0.5 * / (dx*dt) */
                        tmpfy = oowdydt * Y * Sol->rhov[ne];
                        tmpfz = oowdzdt * Y * Sol->rhow[ne];
                        
                        rhs[nn000]       += + tmpfx + tmpfy + tmpfz;
                        rhs[nn010]       += + tmpfx - tmpfy + tmpfz;
                        rhs[nn011]       += + tmpfx - tmpfy - tmpfz;
                        rhs[nn001]       += + tmpfx + tmpfy - tmpfz;
                        rhs[nn100]       += - tmpfx + tmpfy + tmpfz;
                        rhs[nn110]       += - tmpfx - tmpfy + tmpfz;
                        rhs[nn111]       += - tmpfx - tmpfy - tmpfz;
                        rhs[nn101]       += - tmpfx + tmpfy - tmpfz;
                    }
                } 
            } 
            
            
            
            /* account for influx bottom boundary */
            j = igye;
            Y  = mpv->HydroState->Y0[j];
            for(k = igze; k < icze - igze; k++) {
                int ln = k * icyn*icxn;
                int mn = ln + j*icxn;
                for(i = igxe; i < icxe - igxe; i++) {
                    int nn   = mn + i;
                    int nn00 = nn;
                    int nn10 = nn + dixn;
                    int nn11 = nn + dixn + dizn;
                    int nn01 = nn +      + dizn;
                    
                    double rhov_wall = bdry->wall_massflux[i]; 
                    double tmpy = oowdydt * Y * rhov_wall;
                    
                    rhs[nn00] += - tmpy;
                    rhs[nn10] += - tmpy;
                    rhs[nn01] += - tmpy;
                    rhs[nn11] += - tmpy;
                }
            }        
            
            for(k = igze+1; k < icze - igze; k++) {
                int ln = k*icxn*icyn; 
                for(j = igye+1; j < icye - igye; j++) {
                    int mn = ln + j * icxn; 
                    for(i = igxe+1; i < icxe - igxe; i++) {
                        int nn  = mn + i;
                        div_max = MAX_own(div_max, fabs(rhs[nn]));
                    }
                }
            }
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
#endif
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_nodes.hdf");
    sprintf(fieldname, "rhs-nodes");    
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
    
    return div_max;
    
}

/* ========================================================================== */

static enum Constraint integral_condition_nodes(
												double* rhs,
												const NodeSpaceDiscr* node,
												const int x_periodic,
												const int y_periodic,
												const int z_periodic) {
	
	/* Rupert, this is valid only for periodic data so far!!! */
	
	int i, j, k, l, m, n;
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int igz = node->igz;
	const int icz = node->icz;
	const int icxicy = icx * icy;
    
    const int imax = MAX_own(1, icx - igx - x_periodic);
    const int jmax = MAX_own(1, icy - igy - y_periodic);
    const int kmax = MAX_own(1, icz - igz - z_periodic);
    
    const double del = DBL_EPSILON * node->nc;
    
	double tmp = 0.0, rhs_max=0.0;
	
	for(k = igz; k < kmax; k++) {l = k * icxicy;
		for(j = igy; j < jmax; j++) {m = l + j * icx; 
			for(i = igx; i < imax; i++) {n = m + i;
				tmp += rhs[n];
                if(fabs(rhs[n])>rhs_max) rhs_max = fabs(rhs[n]);
			}
		}
	}
    
    /* if(fabs(tmp) > 100*del) { */
    if(fabs(tmp) > 10000*del) {
        printf("CHEATING:  tmp = %e,  rhs_max = %e, (tmp/rhs_max)/del = %e\n", tmp/rhs_max, rhs_max, (tmp/rhs_max)/del);
        /* return VIOLATED; */
        return VIOLATED;
    }
    else if(fabs(tmp) > 100*del) {
        printf("integral_condition_nodes = barely OK; tmp = %e,  rhs_max = %e, (tmp/rhs_max)/del = %e\n", tmp/rhs_max, rhs_max, (tmp/rhs_max)/del);
        return SATISFIED;
    }
    else {
        /* printf("integral_condition_nodes = OK;    tmp = %e,  rhs_max = %e\n", tmp, rhs_max); */
        printf("integral_condition_nodes = OK; tmp = %e\n", tmp);
        return SATISFIED;
    }
}

/* ========================================================================== */

static void	catch_periodic_directions(
									  double* rhs,  
									  const NodeSpaceDiscr* node, 
									  const ElemSpaceDiscr* elem,
									  const int x_periodic,
									  const int y_periodic,
									  const int z_periodic) {
	
	/* in any of the periodicity directions, gather boundary values of
	 rhs on first nodes, set last nodes to zero
	 */
	
	int i, j, k, l, m, l0, l1, m0, m1, n0, n1;
    
	const int igx = node->igx;
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int igz = node->igz;
	const int icz = node->icz;
	const int icxicy = icx * icy;
	
	if(x_periodic == 1){
		for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
			for(j = igy; j < icy - igy; j++) {m = l + j * icx; 
				n0 = m + igx;
				n1 = m + icx-igx-1;
				rhs[n0] += rhs[n1];
				rhs[n1] = 0.0;
			}
		}
	}
	
	if(node->ndim > 1 && y_periodic == 1){
		for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
			m0 = l + igy * icx;
			m1 = l + (icy-igy-1) * icx;
			for(i = igx; i < icx - igx; i++) {
				n0 = m0 + i;
				n1 = m1 + i;
				rhs[n0] += rhs[n1];
				rhs[n1] = 0.0;
			}
		}
	}
	
	if(node->ndim > 2 && z_periodic == 1){
		l0 = igz * icxicy;
		l1 = (icz-igz-1) * icxicy; 
		for(j = igy; j < icy - igy; j++) {
			m0 = l0 + j * icx; 
			m1 = l1 + j * icx; 
			for(i = igx; i < icx - igx; i++) {
				n0 = m0 + i;
				n1 = m1 + i;
				rhs[n0] += rhs[n1];
				rhs[n1] = 0.0;
			}
		}
	}
}

/* ========================================================================== */

static void operator_coefficients_nodes(
										double* hplus[3], 
										double* hcenter,
										const ElemSpaceDiscr* elem,
										const NodeSpaceDiscr* node,
										ConsVars* Sol,
										const ConsVars* Sol0,
										const MPV* mpv,
                                        const double dt) {
	
    extern User_Data ud;
	extern Thermodynamic th;
	
    const double g   = ud.gravity_strength[1];
    const double Msq = ud.Msq;
    const double Gammainv = th.Gammainv;
    
    const int ndim = node->ndim;
                    
	switch(ndim) {
		case 1: {    
			ERROR("interface_enthalpy_nodes() not implemented for 1D\n");
			break;
		}
		case 2: {			
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int igy = elem->igy;
			const int icy = elem->icy;

            const double dy = elem->dy;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hc      = hcenter;

            /*  const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt); 
                const double cexp    = 1.0-th.gamm;
             */
			const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gm1inv/(dt*dt);
            const double cexp    = 2.0-th.gamm;
            
            int i, j, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = 0.0;
						                        
			for(j = igy; j < icy - igy; j++) {
                m = j * icx;
				
                double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j]) \
                                    /(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;

                for(i = igx; i < icx - igx; i++) {
                    n = m + i;     
                    
                    double Y     = 0.5 * (Sol->rhoY[n]/Sol->rho[n] + Sol0->rhoY[n]/Sol0->rho[n]); 
                    double coeff = Gammainv * Sol->rhoY[n] * Y;
                    double Nsqsc = 0.25*dt*dt * (g/Msq) * strat;                    
                    double gimpy = 1.0 / (1.0 + Nsqsc);
                                        
                    hplusx[n]    = coeff;
                    hplusy[n]    = coeff * gimpy;
				}
			}
									
			for(j = igy; j < icy - igy; j++) {m = j * icx;
				for(i = igx; i < icx - igx; i++) {n = m + i;
					hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
				}
			}
			break;
		}
		case 3: {			
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int igz = elem->igz;
			const int icz = elem->icz;
            
            const double dy = elem->dy;

            double Msq = ud.Msq;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hplusz  = hplus[2];
			double* hc      = hcenter;
            
			const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(dt*dt);
            
			const double cexp    = 1.0-th.gamm;
			int i, j, k, l, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = hplusz[i] = 0.0;
            
			for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                    
                    double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j]) \
                                        /(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;

                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        {             
                            double Y     = 0.5 * (Sol->rhoY[n]/Sol->rho[n] + Sol0->rhoY[n]/Sol0->rho[n]); 
                            double coeff = Gammainv * Sol->rhoY[n] * Y;
                            double Nsqsc = 0.25*dt*dt * (g/Msq) * strat;                    
                            double gimpy = 1.0 / (1.0 + Nsqsc);

                            hplusx[n]  = coeff;
                            hplusy[n]  = coeff * gimpy;
                            hplusz[n]  = coeff;
                        }
                    }
                }
			}
            
			for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
                    }
                }
			}
			break;
		}
		default: ERROR("ndim not in {1,2,3}");
	}
}

/* ========================================================================== */

void correction_nodes(
					  ConsVars* Sol,
					  const ElemSpaceDiscr* elem,
					  const NodeSpaceDiscr* node,
                      const double* hplus[3], 
					  double* p, 
					  const double t,
					  const double dt) {
	
	extern User_Data ud;
    extern MPV* mpv;

	const int ndim = elem->ndim;
	    
	switch(ndim) {
		case 1: {
			ERROR("function not available");
			break;
		}
			
		case 2: {
			
			const int igx = node->igx;
			const int icx = node->icx;
			const int igy = node->igy;
			const int icy = node->icy;
			
			const int icxe = node->icx - 1;
			
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];

			const double dx     = node->dx;
			const double dy     = node->dy;
			const double oodx   = 1.0 / dx;
			const double oody   = 1.0 / dy;
            const double dtowdx = 0.5*dt * oodx;
            const double dtowdy = 0.5*dt * oody;                
            			
			int i, j, m, me;
			
			for(j = igy; j < icy - igy - 1; j++) {
				m = j * icx; 
				me = j * icxe;
				for(i = igx; i < icx - igx - 1; i++) {
					const int n     = m + i;
					const int nicx  = n + icx;
					const int n1    = n + 1;
					const int n1icx = n + 1 + icx;
					const int ne    = me + i; 
					
					const double Dpx   = 0.5 * (p[n1]   - p[n] + p[n1icx] - p[nicx]);
					const double Dpy   = 0.5 * (p[nicx] - p[n] + p[n1icx] - p[n1]);
                    const double thinv = Sol->rho[ne] / Sol->rhoY[ne];
					                    
                    Sol->rhou[ne] += - dtowdx * thinv * hplusx[ne] * Dpx;
					Sol->rhov[ne] += - dtowdy * thinv * hplusy[ne] * Dpy;
				}
			} 
			
			break;
		}
        case 3: {
            
            const int igx = node->igx;
            const int icx = node->icx;
            const int igy = node->igy;
            const int icy = node->icy;
            const int igz = node->igz;
            const int icz = node->icz;
            
            const int icxe = node->icx - 1;
            const int icye = node->icy - 1;
                                                
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hplusz = hplus[2];
            
            const double dx     = node->dx;
            const double dy     = node->dy;
            const double dz     = node->dz;
            const double oodx   = 1.0 / dx;
            const double oody   = 1.0 / dy;
            const double oodz   = 1.0 / dz;            
            const double dtowdx = 0.5*dt * oodx;
            const double dtowdy = 0.5*dt * oody;                
            const double dtowdz = 0.5*dt * oodz;                

            for(int k = igz; k < icz-igz-1; k++) {
                int l  = k*icx*icy;
                int le = k*icxe*icye;
                for(int j = igy; j < icy - igy - 1; j++) {
                    int m  = l  + j*icx; 
                    int me = le + j*icxe;
                    for(int i = igx; i < icx - igx - 1; i++) {
                        int n000 = m + i;
                        int n010 = n000 + icx;
                        int n001 = n000 + 1;
                        int n011 = n000 + 1 + icx;
                        int n100 = m + i + icx*icy;
                        int n110 = n100 + icx;
                        int n101 = n100 + 1;
                        int n111 = n100 + 1 + icx;

                        int ne   = me + i; 
                        
                        double Dpx   = 0.25 * (p[n001] - p[n000] + p[n011] - p[n010] + p[n101] - p[n100] + p[n111] - p[n110]);
                        double Dpy   = 0.25 * (p[n010] - p[n000] + p[n011] - p[n001] + p[n110] - p[n100] + p[n111] - p[n101]);
                        double Dpz   = 0.25 * (p[n100] - p[n000] + p[n110] - p[n010] + p[n101] - p[n001] + p[n111] - p[n011]);
                        double thinv = Sol->rho[ne] / Sol->rhoY[ne];
                        
                        Sol->rhou[ne] += - dtowdx * thinv * hplusx[ne] * Dpx;
                        Sol->rhov[ne] += - dtowdy * thinv * hplusy[ne] * Dpy;
                        Sol->rhow[ne] += - dtowdz * thinv * hplusz[ne] * Dpz;
                    }
                } 
            }
            
            break;
        }
		default: ERROR("ndim not in {1, 2, 3}");
	}
}

/* ========================================================================== */

void euler_backward_gravity(ConsVars* Sol,
                            const MPV* mpv,
                            const ElemSpaceDiscr* elem,
                            const double dt)
{
    /* 
     evaluates Euler backward for the gravity term in a semi-implicit 
     trapezoidal or midpoint discretization of the pressure gradient and 
     gravity terms. 
     */
    extern User_Data ud;
    
    const double g   = ud.gravity_strength[1];
    const double Msq = ud.Msq;
    const double dy  = elem->dy;
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
        
    for (int k=0; k<icz; k++) {
        int l = k*icy*icx;
        for (int j=0; j<icy; j++) {
            int m = l + j*icx;
            double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j])/(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;

            for (int i=0; i<icx; i++) {
                int n        = m + i;
                double Nsqsc = dt*dt * (g/Msq) * strat;
                double v     = Sol->rhov[n]/Sol->rho[n];
                double dchi  = Sol->rhoX[BUOY][n]/Sol->rho[n];
                double chi   = Sol->rho[n]/Sol->rhoY[n];
                double dbuoy = -dchi/chi;    /* -dchi/chibar; */
                
                Sol->rhov[n] = Sol->rho[n] * (v + dt * (g/Msq) * dbuoy) / (1.0 + Nsqsc);
            }
        }
    }
    Set_Explicit_Boundary_Data(Sol, elem);
}

/* ========================================================================== */

void euler_forward_non_advective(ConsVars* Sol,
                                 const MPV* mpv,
                                 const ElemSpaceDiscr* elem,
                                 const NodeSpaceDiscr* node,
                                 const double dt,
                                 const enum EXPLICIT_PRESSURE with_pressure)
{
    /* 
     evaluates Euler forward for the pressure gradient, gravity, 
     and background stratification advection terms based on cell-
     centered data. 
     */
    extern User_Data ud;
    extern Thermodynamic th;

    double* p2n = mpv->p2_nodes;
    
    const double g    = ud.gravity_strength[1];
    const double Msq  = ud.Msq;
    const double Ginv = th.Gammainv; 
    
    switch (elem->ndim) {
        case 1:
        {
            const int icx = elem->icx;
            
            const int    igx = elem->igx;
            const double dx  = node->dx;
            
            for (int i=igx; i<icx-igx; i++) {
                int nc = i;
                int n0 = i;
                int n1 = i+1;
                
                double dpdx   = with_pressure*(p2n[n1]-p2n[n0])/dx;
                double rhoYovG = Ginv*Sol->rhoY[nc];
                
                Sol->rhou[nc]       += dt * ( - rhoYovG * dpdx);
            }
        }
            break;
        case 2:
        {
            const int icx = elem->icx;
            const int icy = elem->icy;
            
            const int igx = elem->igx;
            const int igy = elem->igy;
            
            const int inx = node->icx;
            
            const double dx = node->dx;
            const double dy = node->dy;
            
            for (int j=igy; j<icy-igy; j++) {
                int mc = j*icx;
                int mn = j*inx;
                double S0p    = mpv->HydroState_n->S0[j+1];
                double S0m    = mpv->HydroState_n->S0[j];
                
                for (int i=igx; i<icx-igx; i++) {
                    int nc   = mc + i;
                    
                    int n00 = mn  + i;
                    int n10 = n00 + inx;
                    int n01 = n00 + 1;
                    int n11 = n00 + 1 + inx;
                    
                    double dpdx   = with_pressure*0.5*(p2n[n01] - p2n[n00] + p2n[n11] - p2n[n10])/dx;
                    double dpdy   = with_pressure*0.5*(p2n[n10] - p2n[n00] + p2n[n11] - p2n[n01])/dy;
                    double dSdy   = (S0p-S0m) / dy;
                    
                    double rhoYovG = Ginv*Sol->rhoY[nc];
                    double v       = Sol->rhov[nc]/Sol->rho[nc];
                    double dchi    = Sol->rhoX[BUOY][nc]/Sol->rho[nc];
                    double chi     = Sol->rho[nc]/Sol->rhoY[nc];
                    double dbuoy   = -Sol->rho[nc]*dchi/chi;  /* -dchi/chibar; */
                    
                    Sol->rhou[nc]       += dt * ( - rhoYovG * dpdx);
                    Sol->rhov[nc]       += dt * ( - rhoYovG * dpdy + (g/Msq) * dbuoy);
                    Sol->rhoX[BUOY][nc] += dt * ( - v * dSdy) * Sol->rho[nc];
                }
            }
        }
            break;
        case 3: 
        {
            const int icx = elem->icx;
            const int icy = elem->icy;
            const int icz = elem->icz;
            
            const int igx = elem->igx;
            const int igy = elem->igy;
            const int igz = elem->igz;
            
            const int inx = node->icx;
            const int iny = node->icy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            
            for (int k=igz; k<icz-igz; k++) {
                int lc = k*icy*icx;
                int ln = k*iny*inx;
                for (int j=igy; j<icy-igy; j++) {
                    int mc = lc + j*icx;
                    int mn = ln + j*inx;
                    double S0p    = mpv->HydroState_n->S0[j+1];
                    double S0m    = mpv->HydroState_n->S0[j];
                    
                    for (int i=igx; i<icx-igx; i++) {
                        int nc        = mc + i;
                        
                        int n000 = mn   + i;
                        int n010 = n000 + inx;
                        int n001 = n000 + 1;
                        int n011 = n000 + 1 + inx;
                        int n100 = mn   + i + inx*iny;
                        int n110 = n100 + inx;
                        int n101 = n100 + 1;
                        int n111 = n100 + 1 + inx;
                        
                        double dpdx   = with_pressure*0.25*(p2n[n001]-p2n[n000]+p2n[n011]-p2n[n010]+p2n[n101]-p2n[n100]+p2n[n111]-p2n[n110])/dx;
                        double dpdy   = with_pressure*0.25*(p2n[n010]-p2n[n000]+p2n[n011]-p2n[n001]+p2n[n110]-p2n[n100]+p2n[n111]-p2n[n101])/dy;
                        double dpdz   = with_pressure*0.25*(p2n[n100]-p2n[n000]+p2n[n110]-p2n[n010]+p2n[n101]-p2n[n001]+p2n[n111]-p2n[n011])/dz;
                        double dSdy   = (S0p-S0m) / dy;
                        
                        double rhoYovG = Ginv*Sol->rhoY[nc];
                        double v       = Sol->rhov[nc]/Sol->rho[nc];
                        double dchi    = Sol->rhoX[BUOY][nc]/Sol->rho[nc];
                        double chi     = Sol->rho[nc]/Sol->rhoY[nc];
                        double dbuoy   = -Sol->rho[nc]*dchi/chi;  /* -dchi/chibar; */
                        
                        Sol->rhou[nc]       += dt * ( - rhoYovG * dpdx);
                        Sol->rhov[nc]       += dt * ( - rhoYovG * dpdy + (g/Msq) * dbuoy);
                        Sol->rhow[nc]       += dt * ( - rhoYovG * dpdz);
                        Sol->rhoX[BUOY][nc] += dt * ( - v * dSdy) * Sol->rho[nc];
                    }
                }
            }
        }
            break;

        default:
            break;
    }
    
    
    /* last half Euler backward step equals first half Euler forward step */
    if (ud.is_compressible) {
        for (int nn=0; nn<node->nc; nn++) {
            mpv->p2_nodes[nn] += mpv->dp2_nodes[nn];
        }
    }
    Set_Explicit_Boundary_Data(Sol, elem);
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
