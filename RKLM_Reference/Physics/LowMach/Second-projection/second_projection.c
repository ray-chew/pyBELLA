
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
#include "variable.h"
#include "mpv.h"
#include "error.h"
#include "warning.h"
#include "userdata.h"
#include "variable_coefficient_poisson_nodes.h"
#include "second_projection.h" 
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
#include "Hydrostatics.h"
#include "molecular_transport.h"

static double divergence_nodes(
							 double* rhs,
							 const ElemSpaceDiscr* elem,
							 const NodeSpaceDiscr* node,
							 const ConsVars* Sol,
							 const MPV* mpv,
							 const BDRY* bdry);

static void rhs_from_p_old(
                           double* rhs,
                           const ElemSpaceDiscr* elem,
                           const NodeSpaceDiscr* node,
                           const MPV* mpv,
                           const double* hcenter);

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
					  const double dt,
                      const enum CorrectionRange crange);

void bottom_grid_layer_set(double* obj, 
                           const double* src, 
                           const double weight,
                           const ElemSpaceDiscr* elem);

/* Functions needed for the hydrostatic case */
void hydrostatic_vertical_velo(ConsVars* Sol,
                               const ElemSpaceDiscr* elem,
                               const NodeSpaceDiscr* node);

void diss_to_rhs(double* rhs,
                 const double* W0,
                 const ElemSpaceDiscr* elem,
                 const NodeSpaceDiscr* node,
                 const double dt);

/* ========================================================================== */

#if OUTPUT_RHS_NODES
static int rhs_output_count = 0;
static int first_output_step = 0;
extern int step;  
#endif

void euler_backward_non_advective_impl_part(ConsVars* Sol,
                                            MPV* mpv,
                                            const ConsVars* Sol0,
                                            const ElemSpaceDiscr* elem,
                                            const NodeSpaceDiscr* node,
                                            const double t,
                                            const double dt,
                                            const double alpha_diff) {
    
    /* as of August 31, 2018, this routine is to be interpreted
     as carrying out an implicit Euler step for the linearized 
     fast system without advection. This changes in particular
     the interpretation of "dt" from having been the overall 
     time step of the total solver previously to now just being
     whatever time step the implicit Euler step is to be performed
     over. 
     */
    extern User_Data ud;
    extern BDRY* bdry;
    
    const int nc = node->nc;
    
    double** hplus   = mpv->wplus;
    double*  hcenter = mpv->wcenter;
    
    double* rhs      = mpv->rhs;
    double* p2       = mpv->dp2_nodes;
    
    double rhs_max;
    
    int x_periodic, y_periodic, z_periodic;
    int ii;
    
    /* switch for stopping loops in solver for periodic data before
     last grid line in each direction */
    x_periodic = 0;
    y_periodic = 0;
    z_periodic = 0;
    if(ud.bdrytype_min[0] == PERIODIC) x_periodic = 1;
    if(ud.bdrytype_min[1] == PERIODIC) y_periodic = 1;
    if(ud.bdrytype_min[2] == PERIODIC) z_periodic = 1;
    
    printf("\n\n====================================================");
    printf("\nSecond Projection");
    printf("\n====================================================\n");
    
    Set_Explicit_Boundary_Data(Sol, elem);

    /* There is the option of starting with the previous solution, 
     but that seems not to come with sizeable advantages */
    for(ii=0; ii<nc; ii++){
        /* p2[ii] = mpv->p2_nodes[ii];  */
        p2[ii]  = 0.0; 
        rhs[ii] = 0.0;
    }
    
    operator_coefficients_nodes(hplus, hcenter, elem, node, Sol, Sol0, mpv, dt);
    
    if (ud.mol_trans != NO_MOLECULAR_TRANSPORT) {
        extern double* diss;        
        molecular_transport(Sol, diss, elem, alpha_diff*dt);
        diss_to_rhs(rhs,diss,elem,node, dt); /* diss_to_rhs(rhs,W0,elem,node, alpha_diff*dt);  */        
    }

    /* loop for iterating on bottom topography boundary conditions */
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv, bdry);
    // assert(integral_condition_nodes(rhs, node, x_periodic, y_periodic, z_periodic) != VIOLATED);         

    /* rescaling of flux divergence for r.h.s. of the elliptic pressure equation */ 
    for (int i=0; i<node->nc; i++) {
        rhs[i] /= dt;
    }
    printf("\nrhsmax = %e\n", rhs_max);
    
    
#if OUTPUT_RHS_NODES
    FILE *prhsfile = NULL;
    char fn[120], fieldname[90];
    if (step >= first_output_step) {
        if (rhs_output_count < 10) {
            sprintf(fn, "%s/rhs_nodes/rhs_nodes_00%d.hdf", ud.file_name, rhs_output_count);
        } else if(rhs_output_count < 100) {
            sprintf(fn, "%s/rhs_nodes/rhs_nodes_0%d.hdf", ud.file_name, rhs_output_count);
        } else {
            sprintf(fn, "%s/rhs_nodes/rhs_nodes_%d.hdf", ud.file_name, rhs_output_count);
        }
        sprintf(fieldname, "rhs_nodes");    
        WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
        rhs_output_count++;
    }
#endif
    
    if (ud.is_compressible) {
        rhs_from_p_old(rhs, elem, node, mpv, hcenter);
    }
        
    variable_coefficient_poisson_nodes(p2, (const double **)hplus, hcenter, rhs, elem, node, x_periodic, y_periodic, z_periodic, dt);
    
    correction_nodes(Sol, elem, node, (const double**)hplus, p2, dt, FULL_FIELD);
    Set_Explicit_Boundary_Data(Sol, elem);
    
    for(ii=0; ii<nc; ii++) {
        double p2_new      = p2[ii];
        mpv->dp2_nodes[ii] = p2_new - mpv->p2_nodes[ii];
        mpv->p2_nodes[ii]  = p2_new;
    }
    
    set_ghostnodes_p2(mpv->p2_nodes, node, 2);       
    
#if OUTPUT_RHS_NODES
    memset(rhs, 0.0, node->nc*sizeof(double));
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv, bdry);
    /* catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);
     */
    
    if (step >= first_output_step) {
        if (rhs_output_count < 10) {
            sprintf(fn, "%s/rhs_nodes/rhs_nodes_00%d.hdf", ud.file_name, rhs_output_count);
        } else if(rhs_output_count < 100) {
            sprintf(fn, "%s/rhs_nodes/rhs_nodes_0%d.hdf", ud.file_name, rhs_output_count);
        } else {
            sprintf(fn, "%s/rhs_nodes/rhs_nodes_%d.hdf", ud.file_name, rhs_output_count);
        }
        sprintf(fieldname, "rhs_nodes");    
        WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
        rhs_output_count++;
    }
#endif
    
}

/* ========================================================================== */

static double divergence_nodes(
							   double* rhs,
							   const ElemSpaceDiscr* elem,
							   const NodeSpaceDiscr* node,
							   const ConsVars* Sol,
							   const MPV* mpv,
							   const BDRY* bdry) {
	
	extern User_Data ud;
	    
	const int ndim = node->ndim;

    double div_max = 0.0;
    
    int i, j, k, mn;
    
    switch(ndim) {
        case 1: {
            ERROR("divergence_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
    
            int is_x_periodic = 0;
            int is_y_periodic = 0;
            
            if(ud.bdrytype_min[0] == PERIODIC) is_x_periodic = 1;
            if(ud.bdrytype_min[1] == PERIODIC) is_y_periodic = 1;

            const int icxn = node->icx;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double dx      = node->dx;
            const double dy      = node->dy;
            const double oodx    = 1.0 / dx;
            const double oody    = 1.0 / dy;
            
            double Y;
            
            /* predicted time level divergence via scattering */
            for(j = igye - is_y_periodic; j < icye - igye + is_y_periodic; j++) {
                const int me = j * icxe;
                const int mn = j * icxn; 
                for(i = igxe - is_x_periodic; i < icxe - igxe + is_x_periodic; i++) {
                    const int n     = mn + i;
                    const int nicx  = n  + icxn;
                    const int n1    = n  + 1;
                    const int n1icx = n1 + icxn;
                    
                    const int ne = me + i; 
                    double tmpfx, tmpfy;
                    
                    Y = Sol->rhoY[ne] / Sol->rho[ne];
                    tmpfx = 0.5 * oodx * Y * Sol->rhou[ne];
                    tmpfy = 0.5 * oody * Y * Sol->rhov[ne];
                    
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
            for(i = igxe - is_x_periodic; i < icxe - igxe + is_x_periodic; i++) {
                const int n     = mn + i;
                const int n1    = n  + 1;
                
                double rhov_wall = bdry->wall_rhoYflux[i]; 
                double tmpy = 0.5 * oody * Y * rhov_wall; 
                
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
            
            
            int is_x_periodic = 0;
            int is_y_periodic = 0;
            int is_z_periodic = 0;
            
            if(ud.bdrytype_min[0] == PERIODIC) is_x_periodic = 1;
            if(ud.bdrytype_min[1] == PERIODIC) is_y_periodic = 1;
            if(ud.bdrytype_min[2] == PERIODIC) is_z_periodic = 1;

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
                        
            const int dixn = 1;
            const int diyn = icxn;
            const int dizn = icxn*icyn;
            
            double Y;
            
            /* predicted time level divergence via scattering */
            for(k = igze - is_z_periodic; k < icze - igze + is_z_periodic; k++) {
                const int le = k * icye * icxe;
                const int ln = k * icyn * icxn; 
                for(j = igye - is_y_periodic; j < icye - igye + is_y_periodic; j++) {
                    const int me = le + j * icxe;
                    const int mn = ln + j * icxn; 
                    for(i = igxe - is_x_periodic; i < icxe - igxe + is_x_periodic; i++) {
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
                        tmpfx = 0.25 * oodx * Y * Sol->rhou[ne]; 
                        tmpfy = 0.25 * oody * Y * Sol->rhov[ne];
                        tmpfz = 0.25 * oodz * Y * Sol->rhow[ne];
                        
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
            for(k = igze  - is_z_periodic; k < icze - igze + is_z_periodic; k++) {
                int ln = k * icyn*icxn;
                int mn = ln + j*icxn;
                for(i = igxe - is_x_periodic; i < icxe - igxe + is_x_periodic; i++) {
                    int nn   = mn + i;
                    int nn00 = nn;
                    int nn10 = nn + dixn;
                    int nn11 = nn + dixn + dizn;
                    int nn01 = nn +      + dizn;
                    
                    double rhov_wall = bdry->wall_rhoYflux[i]; 
                    double tmpy = 0.25 * oody * Y * rhov_wall;
                    
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
        
    return div_max;
    
}


/* ========================================================================== */

static void rhs_from_p_old(
                           double* rhs,
                           const ElemSpaceDiscr* elem,
                           const NodeSpaceDiscr* node,
                           const MPV* mpv,
                           const double* hcenter) {
        
    extern User_Data ud;
    extern double *W0;
    extern enum Boolean W0_in_use;

    int nodc = 0;
    
    assert(W0_in_use == WRONG);
    W0_in_use = CORRECT;

    double *rhs_hh = W0;
    
    const int ndim = node->ndim;
    
    const int igxn = node->igx;
    const int icxn = node->icx;
    const int igyn = node->igy;
    const int icyn = node->icy;
    
    
    switch(ndim) {
        case 1: {
            ERROR("rhs_from_p_old() not implemented for 1D\n");
            break;
        }
        case 2: {
            
            int is_x_periodic = 0;
            int is_y_periodic = 0;
            
            if(ud.bdrytype_min[0] == PERIODIC) is_x_periodic = 1;
            if(ud.bdrytype_min[1] == PERIODIC) is_y_periodic = 1;

            for (int ic=0; ic<node->nc; ic++) {
                rhs_hh[ic] = 0.0;
            }
            
            for(int j = igyn; j < icyn - igyn; j++) {
                const int mn = j * icxn; 
                for(int i = igxn; i < icxn - igxn; i++) {
                    const int nn = mn + i;
                    rhs_hh[nn] += hcenter[nn]*mpv->p2_nodes[nn];
                }
            } 
            
            for (int ic=0; ic<node->nc; ic++) {
                rhs[ic] += rhs_hh[ic];
            }

            break;
        }
            
        case 3: 
            ERROR("rhs_from_p_old() not implemented for 1D\n");
            break;
        default: ERROR("ndim not in {1, 2, 3}");
    }    
    
    W0_in_use = WRONG;
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
                rhs[n1] = rhs[n0];
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
                rhs[n1] = rhs[n0];
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
                rhs[n1] = rhs[n0];
            }
        }
    }
}

/* ========================================================================== */

void scale_wall_node_values(
                       double* rhs,  
                       const NodeSpaceDiscr* node, 
                       const ElemSpaceDiscr* elem,
                       const double factor) {
    
    /* nodal contol volumes at the boundaries are by factors of
     2, 4, 8, smaller than those within the domain. This is accounted
     for in the present routine for rigid wall boundary conditions.
     */
    extern User_Data ud;
    
    int i, j, k, l, m, l0, l1, m0, m1, n0, n1;
    
    const int igx = node->igx;
    const int icx = node->icx;
    const int igy = node->igy;
    const int icy = node->icy;
    const int igz = node->igz;
    const int icz = node->icz;
    const int icxicy = icx * icy;
    
    if(ud.bdrytype_min[0] == WALL){
        for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
            for(j = igy; j < icy - igy; j++) {m = l + j * icx; 
                n0 = m + igx;
                rhs[n0] *= factor;
            }
        }
    }

    if(ud.bdrytype_max[0] == WALL){
        for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
            for(j = igy; j < icy - igy; j++) {m = l + j * icx; 
                n1 = m + icx-igx-1;
                rhs[n1] *= factor;
            }
        }
    }

    if(node->ndim > 1){
        if(ud.bdrytype_min[1] == WALL){
            for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
                m0 = l + igy * icx;
                for(i = igx; i < icx - igx; i++) {
                    n0 = m0 + i;
                    rhs[n0] *= factor;
                }
            }
        }
        if(ud.bdrytype_max[1] == WALL){
            for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
                m1 = l + (icy-igy-1) * icx;
                for(i = igx; i < icx - igx; i++) {
                    n1 = m1 + i;
                    rhs[n1] *= factor;
                }
            }
        }
    }
    if(node->ndim > 2){
        if(ud.bdrytype_min[2] == WALL){
            l0 = igz * icxicy;
            for(j = igy; j < icy - igy; j++) {
                m0 = l0 + j * icx; 
                for(i = igx; i < icx - igx; i++) {
                    n0 = m0 + i;
                    rhs[n0] *= factor;
                }
            }
        }
        if(ud.bdrytype_max[2] == WALL){
            l1 = (icz-igz-1) * icxicy; 
            for(j = igy; j < icy - igy; j++) {
                m1 = l1 + j * icx; 
                for(i = igx; i < icx - igx; i++) {
                    n1 = m1 + i;
                    rhs[n1] *= factor;
                }
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
    
    double nonhydro = ud.nonhydrostasy;
    
    /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 */
    double time_offset = 3.0 - ud.acoustic_order; 
    // double time_offset = 1.0; 
    
    const double coriolis  = ud.coriolis_strength[0];

                    
	switch(ndim) {
		case 1: {    
			ERROR("operator_coefficients_nodes() not implemented for 1D\n");
			break;
		}
		case 2: {			
            
            int is_x_periodic = 0;
            int is_y_periodic = 0;
            
            if(ud.bdrytype_min[0] == PERIODIC) is_x_periodic = 1;
            if(ud.bdrytype_min[1] == PERIODIC) is_y_periodic = 1;
            
			const int igx = elem->igx;
            const int igy = elem->igy;

            const int icxe = elem->icx;
			const int icye = elem->icy;

            const double dy = elem->dy;
			
            const int icxn = node->icx;
            const int icyn = node->icy;
            
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hc      = hcenter;

			const double ccenter = - (ud.compressibility*ud.Msq)*th.gm1inv/(dt*dt)/time_offset;
            const double cexp    = 2.0-th.gamm;
            
            int i, j, m, n;
            

            int nodc = 0;
            
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = 0.0;
            
            for(j = igy - is_y_periodic - nodc; j < icye - igy + is_y_periodic + nodc; j++) {
                m = j * icxe;
                
                double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j]) \
                /(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;
                
                for(i = igx - is_x_periodic - nodc; i < icxe - igx + is_x_periodic + nodc; i++) {
                    n = m + i;     
                    double Y     = Sol->rhoY[n]/Sol->rho[n]; 
                    double coeff = Gammainv * Sol->rhoY[n] * Y;
                    double fsqsc = dt*dt * coriolis*coriolis;
                    double fimp  = 1.0 / (1.0 + fsqsc);
                    double Nsqsc = time_offset * dt*dt * (g/Msq) * strat;                    
                    double gimp  = 1.0 / (nonhydro + Nsqsc);
                    
                    hplusx[n]    = coeff * fimp;
                    hplusy[n]    = coeff * gimp;
                }
            }
            
            for(j = igy; j < icyn - igy; j++) {
                int mn = j * icxn;
                int me = j * icxe;
                for(int i = igx; i < icxn - igx; i++) {
                    int nn = mn + i;
                    int ne00 = me + i;
                    int ne10 = me + i - 1;
                    int ne01 = me + i - icxe;
                    int ne11 = me + i -icxe - 1;
                    hc[nn] = ccenter * 0.25*(pow(Sol->rhoY[ne00],cexp) + pow(Sol->rhoY[ne01],cexp) + pow(Sol->rhoY[ne10],cexp) + pow(Sol->rhoY[ne11],cexp));
                }
            }
            
            scale_wall_node_values(hc, node, elem, 0.5);
            
            break;
		}
		case 3: {			

            int is_x_periodic = 0;
            int is_y_periodic = 0;
            int is_z_periodic = 0;
            
            if(ud.bdrytype_min[0] == PERIODIC) is_x_periodic = 1;
            if(ud.bdrytype_min[1] == PERIODIC) is_y_periodic = 1;
            if(ud.bdrytype_min[2] == PERIODIC) is_z_periodic = 1;
            

			const int igx = elem->igx;
            const int igy = elem->igy;
            const int igz = elem->igz;

            const int icxe = elem->icx;
			const int icye = elem->icy;
			const int icze = elem->icz;
            
            const int icxn = node->icx;
            const int icyn = node->icy;
            const int iczn = node->icz;
            
            const double dy = elem->dy;

            double Msq = ud.Msq;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hplusz  = hplus[2];
			double* hc      = hcenter;
            
			const double ccenter = - (ud.compressibility*ud.Msq)*th.gamminv/(dt*dt)/time_offset;
            
			const double cexp    = 1.0-th.gamm;
            
            for(int i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = hplusz[i] = 0.0;
            
			for(int k = igz - is_z_periodic; k < icze - igz + is_z_periodic; k++) {
                int l = k * icxe*icye;
                
                for(int j = igy - is_y_periodic; j < icye - igy + is_y_periodic; j++) {
                    int m = l + j * icxe;
                    
                    double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j]) \
                                        /(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;

                    for(int i = igx - is_x_periodic; i < icxe - igx + is_x_periodic; i++) {
                        int n = m + i;
                        double Y     = Sol->rhoY[n]/Sol->rho[n]; 
                        double coeff = Gammainv * Sol->rhoY[n] * Y;
                        double fsqsc = dt*dt * coriolis*coriolis;
                        double fimp  = 1.0 / (1.0 + fsqsc);
                        double Nsqsc = time_offset * dt*dt * (g/Msq) * strat;                    
                        double gimp  = 1.0 / (nonhydro + Nsqsc);
                        
                        hplusx[n]  = coeff * fimp;
                        hplusy[n]  = coeff * gimp;
                        hplusz[n]  = coeff * fimp;
                    }
                }
			}
            
            
            for(int k = igz; k < iczn - igz; k++) {
                int ln = k * icxn*icyn;
                int le = k * icxe*icye;
                for(int j = igy; j < icyn - igy; j++) {
                    int mn = ln + j * icxn;
                    int me = le + j * icxe;
                    for(int i = igx; i < icxn - igx; i++) {
                        int nn = mn + i;
                        int ne000 = me + i;
                        int ne100 = me + i - 1;
                        int ne010 = me + i - icxe;
                        int ne110 = me + i -icxe - 1;
                        int ne001 = me + i;
                        int ne101 = me + i - 1;
                        int ne011 = me + i - icxe;
                        int ne111 = me + i -icxe - 1;
                        hc[nn] = ccenter * 0.125*(  pow(Sol->rhoY[ne000],cexp) \
                                                  + pow(Sol->rhoY[ne010],cexp) \
                                                  + pow(Sol->rhoY[ne100],cexp) \
                                                  + pow(Sol->rhoY[ne110],cexp) \
                                                  + pow(Sol->rhoY[ne001],cexp) \
                                                  + pow(Sol->rhoY[ne011],cexp) \
                                                  + pow(Sol->rhoY[ne101],cexp) \
                                                  + pow(Sol->rhoY[ne111],cexp));
                    }
                }
            }
            
            scale_wall_node_values(hc, node, elem, 0.5);

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
                      const double dt,
                      const enum CorrectionRange crange) {
	
	extern User_Data ud;
    extern MPV* mpv;

	const int ndim = elem->ndim;
    	    
    const double coriolis  = ud.coriolis_strength[0];
    
    /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 */
    double time_offset = 3.0 - ud.acoustic_order; 
    // double time_offset = 1.0; 


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
			
			const double dx     = node->dx;
			const double dy     = node->dy;
			const double oodx   = 1.0 / dx;
			const double oody   = 1.0 / dy;
                        			
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            
            int i, j, m, me;
			
            int jmax  = icy - igy - 1;
            int vcorr = 1;
            
            if (crange == BOTTOM_ONLY) {
                jmax  = igy+1;
                vcorr = 0;
            }
            
			for(j = igy; j < jmax; j++) {
				m = j * icx; 
				me = j * icxe;
                
                double dSdy = (mpv->HydroState_n->S0[j+1] - mpv->HydroState_n->S0[j]) * oody;

				for(i = igx; i < icx - igx - 1; i++) {
					const int n     = m + i;
					const int nicx  = n + icx;
					const int n1    = n + 1;
					const int n1icx = n + 1 + icx;
					const int ne    = me + i; 
					
					const double Dpx   = 0.5 * oodx * (p[n1]   - p[n] + p[n1icx] - p[nicx]);
					const double Dpy   = 0.5 * oody * (p[nicx] - p[n] + p[n1icx] - p[n1]);
                    const double thinv = Sol->rho[ne] / Sol->rhoY[ne];
					                    
                    Sol->rhou[ne] += - dt * thinv * hplusx[ne] * Dpx;
					Sol->rhov[ne] += - dt * thinv * hplusy[ne] * Dpy * vcorr;
                    /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 
                     the following line was not there on Oct. 24. */
                    Sol->rhoX[BUOY][ne] += - time_offset * dt * dSdy * Sol->rhov[ne] * vcorr;
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
            
            const double dx     = node->dx;
            const double dy     = node->dy;
            const double dz     = node->dz;
            const double oodx   = 1.0 / dx;
            const double oody   = 1.0 / dy;
            const double oodz   = 1.0 / dz;
                        
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hplusz = hplus[2];
                        
            int jmax = icy - igy - 1;
            if (crange == BOTTOM_ONLY) {
                jmax = igy+1;
            }

            for(int k = igz; k < icz-igz-1; k++) {
                int l  = k*icx*icy;
                int le = k*icxe*icye;
                
                for(int j = igy; j < jmax; j++) {
                    int m  = l  + j*icx; 
                    int me = le + j*icxe;

                    double dSdy = (mpv->HydroState_n->S0[j+1] - mpv->HydroState_n->S0[j]) * oody;

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
                        
                        double Dpx   = 0.25 * oodx * (p[n001] - p[n000] + p[n011] - p[n010] + p[n101] - p[n100] + p[n111] - p[n110]);
                        double Dpy   = 0.25 * oody * (p[n010] - p[n000] + p[n011] - p[n001] + p[n110] - p[n100] + p[n111] - p[n101]);
                        double Dpz   = 0.25 * oodz * (p[n100] - p[n000] + p[n110] - p[n010] + p[n101] - p[n001] + p[n111] - p[n011]);
                        double thinv = Sol->rho[ne] / Sol->rhoY[ne];
                        
                        Sol->rhou[ne] += - dt * thinv * hplusx[ne] * (Dpx + dt * coriolis * Dpz);
                        Sol->rhov[ne] += - dt * thinv * hplusy[ne] * Dpy;
                        Sol->rhow[ne] += - dt * thinv * hplusz[ne] * (Dpz - dt * coriolis * Dpx);
                        /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 
                         the following line was not there on October 24 */
                        Sol->rhoX[BUOY][ne] += - time_offset * dt * dSdy * Sol->rhov[ne];
                    }
                } 
            }
            
            break;
        }
		default: ERROR("ndim not in {1, 2, 3}");
	}
}



/* ========================================================================== */

void euler_backward_non_advective_expl_part(ConsVars* Sol,
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
    
    double nonhydro = ud.nonhydrostasy;
    
    const double g   = ud.gravity_strength[1];
    const double Msq = ud.Msq;
    const double dy  = elem->dy;
    
    /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 */
    const double time_offset = 3.0 - ud.acoustic_order;
    // const double time_offset = 1.0;
    
    const double coriolis  = ud.coriolis_strength[0];
    const double u0        = ud.wind_speed;
    const double fsqsc     = dt*dt*coriolis*coriolis;
    const double ooopfsqsc = 1.0 / (1.0 + fsqsc);
    
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
                
                /* implicit gravity */
                double Nsqsc = time_offset * dt*dt * (g/Msq) * strat;
                double dbuoy = -Sol->rhoY[n]*Sol->rhoX[BUOY][n]/Sol->rho[n];    
                double rhov  = (nonhydro * Sol->rhov[n] + dt * (g/Msq) * dbuoy) / (nonhydro + Nsqsc);
                
                /* implicit Coriolis */
                double drhou = Sol->rhou[n] - u0*Sol->rho[n];
                Sol->rhou[n] = u0*Sol->rho[n] + ooopfsqsc * (drhou + dt * coriolis * Sol->rhow[n]);
                Sol->rhov[n] = rhov;
                Sol->rhow[n] = ooopfsqsc * (Sol->rhow[n] - dt * coriolis * drhou);
            }
        }
    }
    Set_Explicit_Boundary_Data(Sol, elem);
}

/* ========================================================================== */

void euler_forward_non_advective(ConsVars* Sol,
                                 MPV* mpv,
                                 const ConsVars* Sol0,
                                 const ElemSpaceDiscr* elem,
                                 const NodeSpaceDiscr* node,
                                 const double dt,
                                 const enum EXPLICIT_PRESSURE wp)
{
    /* 
     evaluates Euler forward for the pressure gradient, gravity, 
     and background stratification advection terms based on cell-
     centered data. 
     */
    extern User_Data ud;
    extern Thermodynamic th;
    extern BDRY* bdry;

    extern double *W0;
    extern enum Boolean W0_in_use;

    double nonhydro = ud.nonhydrostasy;

    /*
    double* dp2n = mpv->dp2_nodes;
     */
    assert(W0_in_use == WRONG);
    W0_in_use = CORRECT;
    double* dp2n = W0;

    double* p2n  = mpv->p2_nodes;
    const double g        = ud.gravity_strength[1];
    const double Msq      = ud.Msq;
    const double Ginv     = th.Gammainv;     
    const double coriolis = ud.coriolis_strength[0];
    const double u0       = ud.wind_speed;

    double *div = mpv->rhs;
    double div_max;
            
    memset(dp2n, 0.0, node->nc*sizeof(double));
    memset(div, 0.0, node->nc*sizeof(double));

    /* TODO: call to divergence_nodes() might be unnecessary at least after the first
     step, because mpv->dp2_nodes[] should contain the relevant information already
     from the last time step 
     */
    /* last two arguments to  divergence_nodes() tuned to pure divergence calculation */
    int x_periodic, y_periodic, z_periodic;    
    x_periodic = 0;
    y_periodic = 0;
    z_periodic = 0;
    if(ud.bdrytype_min[0] == PERIODIC) x_periodic = 1;
    if(ud.bdrytype_min[1] == PERIODIC) y_periodic = 1;
    if(ud.bdrytype_min[2] == PERIODIC) z_periodic = 1;

    div_max = divergence_nodes(div, elem, node, (const ConsVars*)Sol, mpv, bdry);
    // catch_periodic_directions(div, node, elem, x_periodic, y_periodic, z_periodic);
    scale_wall_node_values(div, node, elem, 2.0);

#if 0
    div_max = 0.0;
    for (int nn=0; nn<node->nc; nn++) {
        div_max = MAX_own(div_max, div[nn]);
    }
    printf("\n div_max = %e -- euler_forward_non_advective()", div_max);
#endif
    
    switch (elem->ndim) {
        case 1:
        {
            const int icx   = elem->icx;
            const int igx   = elem->igx;
            const double dx = node->dx;
                        
            for (int i=igx; i<icx-igx+1; i++) {
                int nc = i;
                int nn0 = i;
                int nn1 = i+1;
                
                double dpdx    = wp*(p2n[nn1]-p2n[nn0])/dx;
                double rhoYovG = Ginv*Sol->rhoY[nc];
                double drhou   = Sol->rhou[nc] - u0*Sol->rho[nc];
                double dpidP   = (th.gm1 / ud.Msq) * \
                                0.5 * (pow(Sol->rhoY[nc], th.gamm - 2.0) + pow(Sol->rhoY[nc-1], th.gamm - 2.0));

                Sol->rhou[nc]  = u0*Sol->rho[nc] + dt * (- rhoYovG * dpdx + coriolis * Sol->rhow[nc]);
                Sol->rhow[nc]  = Sol->rhow[nc] - dt * coriolis * drhou;
                
                dp2n[nn0] -= dt * dpidP * div[nn0];
            }
            ERROR("boundary fix in  euler_forward_non_advective()  not implemented in 1D yet\n");           
        }
            break;
        case 2:
        {
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;

            const int icxn = node->icx;
            
            const double dx = node->dx;
            const double dy = node->dy;
            
            for (int j=igye; j<icye-igye+1; j++) {
                int mc = j*icxe;
                int mn = j*icxn;
                
                double S0p = mpv->HydroState_n->S0[j+1];
                double S0c = mpv->HydroState->S0[j];
                double S0m = mpv->HydroState_n->S0[j];
                
                for (int i=igxe; i<icxe-igxe+1; i++) {
                    int nc  = mc + i;
                    
                    int nn00 = mn  + i;
                    int nn10 = nn00 + icxn;
                    int nn01 = nn00 + 1;
                    int nn11 = nn00 + 1 + icxn;
   
                    double dpdx    = wp*0.5*(p2n[nn01]-p2n[nn00]+p2n[nn11]-p2n[nn10])/dx;
                    double dpdy    = wp*0.5*(p2n[nn10]-p2n[nn00]+p2n[nn11]-p2n[nn01])/dy;
                    double dSdy    = (S0p-S0m) / dy;
                    double rhoYovG = Ginv*Sol->rhoY[nc];
                    double v       = Sol->rhov[nc]/Sol->rho[nc];
                    double dchi    = Sol->rhoX[BUOY][nc]/Sol->rho[nc];
                    double dbuoy   = -Sol->rhoY[nc]*dchi;  
                    double drhou   = Sol->rhou[nc] - u0*Sol->rho[nc];

                    double dpidP   = (th.gm1 / ud.Msq) * \
                                        0.25 * (pow(Sol->rhoY[nc], th.gamm - 2.0)      + pow(Sol->rhoY[nc-1], th.gamm - 2.0) + \
                                                pow(Sol->rhoY[nc-icxe], th.gamm - 2.0) + pow(Sol->rhoY[nc-icxe-1], th.gamm - 2.0));

                    double time_offset_expl = ud.acoustic_order - 1.0;
                    Sol->rhou[nc]  = Sol->rhou[nc] + dt * ( - rhoYovG * dpdx + coriolis * Sol->rhow[nc]);
                    Sol->rhov[nc]  = Sol->rhov[nc] + dt * ( - rhoYovG * dpdy + (g/Msq) * dbuoy) * nonhydro; 
                    Sol->rhow[nc]  = Sol->rhow[nc] - dt * coriolis * drhou;
                    Sol->rhoX[BUOY][nc] = (Sol->rho[nc] * ( Sol->rho[nc]/Sol->rhoY[nc] - S0c)) +  time_offset_expl * dt * ( - v * dSdy) * Sol->rho[nc];

                    
                    dp2n[nn00] -= dt * dpidP * div[nn00];
                }
            }
        }
            break;
        case 3: 
        {
            const int icx = elem->icx;
            const int icy = elem->icy;
            const int icz = elem->icz;
            const int icxy = icx*icy;
            
            const int igx = elem->igx;
            const int igy = elem->igy;
            const int igz = elem->igz;
            
            const int inx = node->icx;
            const int iny = node->icy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            
            for (int k=igz; k<icz-igz+1; k++) {
                int lc = k*icy*icx;
                int ln = k*iny*inx;
                for (int j=igy; j<icy-igy+1; j++) {
                    int mc = lc + j*icx;
                    int mn = ln + j*inx;
                    
                    double S0p = mpv->HydroState_n->S0[j+1];
                    double S0c = mpv->HydroState->S0[j];
                    double S0m = mpv->HydroState_n->S0[j];
                    
                    for (int i=igx; i<icx-igx+1; i++) {
                        int nc        = mc + i;
                        
                        int nn000 = mn   + i;
                        int nn010 = nn000 + inx;
                        int nn001 = nn000 + 1;
                        int nn011 = nn000 + 1 + inx;
                        int nn100 = mn   + i + inx*iny;
                        int nn110 = nn100 + inx;
                        int nn101 = nn100 + 1;
                        int nn111 = nn100 + 1 + inx;
                        
                        double dpdx   = wp*0.25*(p2n[nn001]-p2n[nn000]+p2n[nn011]-p2n[nn010]+p2n[nn101]-p2n[nn100]+p2n[nn111]-p2n[nn110])/dx;
                        double dpdy   = wp*0.25*(p2n[nn010]-p2n[nn000]+p2n[nn011]-p2n[nn001]+p2n[nn110]-p2n[nn100]+p2n[nn111]-p2n[nn101])/dy;
                        double dpdz   = wp*0.25*(p2n[nn100]-p2n[nn000]+p2n[nn110]-p2n[nn010]+p2n[nn101]-p2n[nn001]+p2n[nn111]-p2n[nn011])/dz;
                        double dSdy   = (S0p-S0m) / dy;
                        
                        double rhoYovG = Ginv*Sol->rhoY[nc];
                        double v       = Sol->rhov[nc]/Sol->rho[nc];
                        double dchi    = Sol->rhoX[BUOY][nc]/Sol->rho[nc];
                        double chi     = Sol->rho[nc]/Sol->rhoY[nc];
                        double dbuoy   = -Sol->rho[nc]*dchi/chi;  /* -dchi/chibar; */
                        double drhou   = Sol->rhou[nc] - u0*Sol->rho[nc];

                        double dpidP   = (th.gm1 / ud.Msq) * \
                        0.125 * (pow(Sol->rhoY[nc], th.gamm - 2.0)          + pow(Sol->rhoY[nc-1], th.gamm - 2.0) + \
                                 pow(Sol->rhoY[nc-icx], th.gamm - 2.0)      + pow(Sol->rhoY[nc-icx-1], th.gamm - 2.0) + \
                                 pow(Sol->rhoY[nc-icxy], th.gamm - 2.0)     + pow(Sol->rhoY[nc-1-icxy], th.gamm - 2.0) + \
                                 pow(Sol->rhoY[nc-icx-icxy], th.gamm - 2.0) + pow(Sol->rhoY[nc-icx-1-icxy], th.gamm - 2.0));
                        
                        double time_offset_expl = ud.acoustic_order - 1.0;
                        Sol->rhou[nc]  = Sol->rhou[nc] + dt * ( - rhoYovG * dpdx + coriolis * Sol->rhow[nc]);
                        Sol->rhov[nc]  = Sol->rhov[nc] + dt * ( - rhoYovG * dpdy + (g/Msq) * dbuoy) * nonhydro; 
                        Sol->rhow[nc]  = Sol->rhow[nc] + dt * ( - rhoYovG * dpdz - coriolis * drhou);
                        Sol->rhoY[nc]  = Sol->rhoY[nc] - dt * div[nc];
                        Sol->rhoX[BUOY][nc] = (Sol->rho[nc] * ( Sol->rho[nc]/Sol->rhoY[nc] - S0c)) + time_offset_expl * dt * ( - v * dSdy) * Sol->rho[nc];
                        dp2n[nn000] -= dt * dpidP * div[nn000];
                    }
                }
            }
            ERROR("boundary fix in  euler_forward_non_advective()  not implemented in 3D yet\n");
        }
            break;

        default:
            break;
    }
    
    if (ud.is_compressible){
        double weight = (ud.acoustic_order - 1.0);
        for (int nn=0; nn<node->nc; nn++) {
            mpv->p2_nodes[nn]  += weight*dp2n[nn];
            /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 
             the following line did not exist in the October 24 version */            
            // mpv->dp2_nodes[nn]  = weight*dp2n[nn]; 
        }
    }
     
    W0_in_use = WRONG;

    set_ghostnodes_p2(mpv->p2_nodes, node, 2);       
    Set_Explicit_Boundary_Data(Sol, elem);
    
}


/* ========================================================================== */

void bottom_grid_layer_set(double* obj, 
                           const double* src, 
                           const double weight,
                           const ElemSpaceDiscr* elem)
{
    /* copy bottom grid layer (i.e., j=elem->igy) data from src to obj */
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    const int igy = elem->igy;
    
    for (int k=0; k<icz; k++) {
        int nk  = k*icx*icy;
        int njk = nk + igy*icx;
        for (int i=0; i<icx; i++) {
            int nijk = njk + i;
            obj[nijk] = weight*src[nijk] + (1.0-weight)*obj[nijk];
        }
    }
}

/* ========================================================================== */

void hydrostatic_vertical_velo(ConsVars* Sol,
                               const ElemSpaceDiscr* elem,
                               const NodeSpaceDiscr* node)
{
    extern double* W0;
    extern enum Boolean W0_in_use;
    
    /* In the hydrostatic case, the vertical rhoY flux can be computed
     directly from the divergence of the horizontal one
     */
    
    /* implemented thus far only for a vertical slice */
    assert(elem->ndim == 2);
    
    assert(W0_in_use == WRONG);
    W0_in_use = CORRECT;
    
    double *rhoYv = W0;
    
    /* get div-controlled vertical velocity on grid edges first */
    const int igxe = elem->igx;
    const int igye = elem->igy;
    
    const int icxe = elem->icx;
    const int icye = elem->icy;

    const int igxn = node->igx;
    const int igyn = node->igy;
    
    const int icxn = node->icx;
    const int ifyn = node->ify;

    const double dyovdx = elem->dy/elem->dx;
    
    for (int i=igxn; i<icxn-igxn; i++) {
        int mc0  = i;
        int mfy  = i*ifyn;
        int nc00 = mc0 + igye*icxe;
        int nc0m = nc00-1;
        int nfy  = mfy + igyn + 1;
    
        int ncm0, ncmm;
        
        double diverr;
        
        double rhoYup = 0.5*Sol->rhoY[nc00]*Sol->rhou[nc00]/Sol->rho[nc00];
        double rhoYum = 0.5*Sol->rhoY[nc0m]*Sol->rhou[nc0m]/Sol->rho[nc0m];
        
        rhoYv[nfy] = - dyovdx * (rhoYup - rhoYum);
        
        for (int j=igyn+2; j<ifyn-igyn-1; j++) {
            nfy = mfy + j;
            nc00 = mc0 + (j-1)*icxe;
            nc0m = mc0 + (j-1)*icxe - 1;
            ncm0 = mc0 + (j-1)*icxe - icxe;
            ncmm = mc0 + (j-1)*icxe - 1 - icxe;
            rhoYup = 0.5*( Sol->rhoY[nc00]*Sol->rhou[nc00]/Sol->rho[nc00]
                          +Sol->rhoY[ncm0]*Sol->rhou[ncm0]/Sol->rho[ncm0]);
            rhoYum = 0.5*( Sol->rhoY[nc0m]*Sol->rhou[nc0m]/Sol->rho[nc0m]
                          +Sol->rhoY[ncmm]*Sol->rhou[ncmm]/Sol->rho[ncmm]);
            rhoYv[nfy] = rhoYv[nfy-1] - dyovdx * (rhoYup - rhoYum);
        }
        /* test whether the vertically integrated divergence vanishes */
        nfy  = mfy + ifyn-igyn-2;
        ncm0 = mc0 + (icye-igxe-1)*icxe;
        ncmm = ncm0 - 1;
        rhoYup = 0.5*Sol->rhoY[ncm0]*Sol->rhou[ncm0]/Sol->rho[ncm0];
        rhoYum = 0.5*Sol->rhoY[ncmm]*Sol->rhou[ncmm]/Sol->rho[ncmm];
        diverr = rhoYv[nfy] - dyovdx * (rhoYup - rhoYum);
        //  printf("diverr = %e\n", diverr);
    }
    
    /* cell-centered momenta from vertical rhoY fluxes on the vertical edges */
    for (int i=igxe; i<icxe-igxe; i++) {
        int mc = i;
        int mf = i*ifyn;
        for (int j=igye; j<icye-igxe; j++) {
            int nc = mc + j*icxe;
            int nf = mf + j+1;
            Sol->rhov[nc] = 0.5*(rhoYv[nf] + rhoYv[nf+ifyn])*Sol->rho[nc]/Sol->rhoY[nc];
        }
    }
    
    W0_in_use = WRONG;
}

/* ========================================================================== */

void diss_to_rhs(double* rhs,
                 const double* diss,
                 const ElemSpaceDiscr* elem,
                 const NodeSpaceDiscr* node,
                 const double dt)
{
    /* here we construct the right hand side of the potential temperature 
     balance equation from energy dissipation */
    extern User_Data ud;
    
    if (ud.mol_trans == FULL_MOLECULAR_TRANSPORT) {
        
        assert(0); /* not implemented yet */
    
    } else if (ud.mol_trans == STRAKA_DIFFUSION_MODEL) {
    
        switch (elem->ndim) {
            case 1:
                assert(0);  /* not implemented yet */ 
                break;

            case 2: {
                int icxe = elem->icx;
                int icye = elem->icy;
                int igxe = elem->igx;
                int igye = elem->igy;
                
                int icxn = node->icx;
                int igxn = node->igx;
                
                /* scatter the elem-based  diss/dt  to  the node-based  rhs */ 
                for (int j=igye; j<icye-igye; j++) {
                    int ncj = j*icxe;
                    int nnj = j*icxn;
                    for (int i=igxe; i<icxe-igxe; i++) {
                        int ncij = ncj + i;
                        int nnij = nnj + i;
                        int nnijmm = nnij;
                        int nnijpm = nnij + 1;
                        int nnijpp = nnij + 1 + icxn;
                        int nnijmp = nnij + icxn;
                        
                        double drhs = 0.25*diss[ncij]/dt;
                        
                        rhs[nnijmm] -= drhs;
                        rhs[nnijpm] -= drhs;
                        rhs[nnijpp] -= drhs;
                        rhs[nnijmp] -= drhs;
                    }
                }
                if (ud.bdrytype_max[0] == PERIODIC) {
                    for (int j=igye; j<icye-igye; j++) {
                        int ncj    = j*icxe;
                        int ncijl  = ncj + igxe;
                        int ncijr  = ncj + icxe-igxe-1;
                        
                        int nnj    = j*icxn;
                        int nnijml = nnj + igxn;
                        int nnijpl = nnj + igxn + icxn;
                        int nnijmr = nnj + (icxn-igxn-1);
                        int nnijpr = nnj + (icxn-igxn-1) + icxn;
                        
                        double drhsl = 0.25*diss[ncijl];
                        double drhsr = 0.25*diss[ncijr];
                        
                        rhs[nnijml] += drhsl;
                        rhs[nnijpl] += drhsl;
                        rhs[nnijmr] += drhsr;
                        rhs[nnijpr] += drhsr;
                    }
                }
                if (ud.bdrytype_max[1] == PERIODIC) {
                    assert(0);
                }
            }
                break;

            case 3:
                assert(0);  /* not implemented yet */
                break;

            default:
                break;
        }  
    } 
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

