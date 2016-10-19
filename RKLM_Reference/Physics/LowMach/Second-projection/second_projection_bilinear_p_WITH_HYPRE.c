
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
#include "set_ghostnodes_p.h"
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
#include "set_ghostcells_p.h"

static double divergence_nodes(
							 double* rhs,
							 const ElemSpaceDiscr* elem,
							 const NodeSpaceDiscr* node,
							 const ConsVars* Sol,
                             double* eta,
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
										double* wgrav,
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
                      const double* hgrav,
					  double* p2, 
					  const double t,
					  const double dt);

/* ========================================================================== */

void second_projection(
                       ConsVars* Sol,
                       MPV* mpv,
                       const ConsVars* Sol0,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       const double p_update,
                       const double t,
                       const double dt) 
{
	
	extern User_Data ud;
	extern BDRY* bdry;
        
	const int nc = node->nc;
	
	double** hplus   = mpv->Level[0]->wplus;
	double*  hcenter = mpv->Level[0]->wcenter;
	double*  hgrav   = mpv->Level[0]->wgrav;
	
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
		p2[ii] = mpv->dp2_nodes[ii];
		rhs[ii] = 0.0;
	}
	    
    double rhs_weight_new = 2.0;
    double rhs_weight_old = 2.0*ud.compressibility;
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, rhs_weight_new);
    if (ud.compressibility) {
        divergence_nodes(rhs, elem, node, Sol0, mpv->eta0, mpv, bdry, dt, rhs_weight_old); /*  for psinc, this can be commented out */
    }
        
	catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);

    /* test */
#if 1
    FILE *prhsfile = NULL;
    char fn[120], fieldname[90];
    sprintf(fn, "%s/rhs_nodes/rhs_nodes_000.hdf", ud.file_name);
    sprintf(fieldname, "rhs-nodes");
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
     
    
	assert(integral_condition_nodes(rhs, node, x_periodic, y_periodic, z_periodic) != VIOLATED); 
	operator_coefficients_nodes(hplus, hcenter, hgrav, elem, node, Sol, Sol0, mpv, dt);
	
	variable_coefficient_poisson_nodes(p2, (const double **)hplus, hcenter, hgrav, rhs, x_periodic, y_periodic, z_periodic, dt);
    
#if 1
    rhs_max = 0.0;
    double rhs_min = 100000.0;
    double hx_max  = 0.0;
    double hy_max  = 0.0;
    double p2_max  = 0.0;
    double hx_min  = 100000.0;
    double hy_min  = 100000.0;
    double p2_min  = 100000.0;
    for (int in = 0; in<node->nc; in++) {
        rhs_max = MAX_own(rhs_max, fabs(rhs[in]));
        rhs_min = MIN_own(rhs_min, fabs(rhs[in]));
    }
    for (int jn = elem->igy; jn < elem->icy-elem->igy; jn++) {
        for (int in = elem->igx; in < elem->icx-elem->igx; in++) {
            int nn = jn*elem->icx + in;
            hx_max = MAX_own(hx_max, fabs(hplus[0][nn]));
            hx_min = MIN_own(hx_min, fabs(hplus[0][nn]));
            hy_max = MAX_own(hy_max, fabs(hplus[1][nn]));
            hy_min = MIN_own(hy_min, fabs(hplus[1][nn]));
        }
    }
    for (int jn = node->igy; jn < node->icy-node->igy; jn++) {
        for (int in = node->igx; in < node->icx-node->igx; in++) {
            int nn = jn*node->icx + in;
            p2_max = MAX_own(p2_max, p2[nn]);
            p2_min = MIN_own(p2_min, p2[nn]);
        }
    }
    printf("rhs_max, rhs_min, h1_max, h1_min, h2_max, h2_max, p2_min, p2_max = %e, %e, %e, %e, %e, %e, %e, %e\n", rhs_max, rhs_min, hx_max, hx_min, hy_max, hy_min, p2_min, p2_max);

    sprintf(fn, "%s/rhs_nodes/rhs_nodes_001.hdf", ud.file_name);
    sprintf(fieldname, "rhs-nodes");
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);

#endif
        
#if 1
    extern double *W0;
    int imax, jmax;
    rhs_max = 0.0;
    for (int jn=node->igy; jn<node->icy-node->igy; jn++) {
        int mn = jn*node->icx;
        for (int in=node->igx; in<node->icx-node->igx; in++) {
            int nn = mn+in;
            if (rhs_max < fabs(rhs[nn])) {
                rhs_max = fabs(rhs[nn]);
                imax = in;
                jmax = jn;
            } 
        }
    }
    printf("rhs_max before = %e\n", rhs_max);
    correction_nodes(Sol, elem, node, (const double **)hplus, hgrav, p2, t, dt);
    EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, p2, (const double **)hplus, hcenter, hgrav, x_periodic, y_periodic, z_periodic, W0);
    precon_apply(W0,W0,node);
    for(ii=0; ii<nc; ii++) rhs[ii] = 0.0;
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, rhs_weight_new);
    catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);

#if 1
    FILE *prhs2file = NULL;
    char fn2[120], fieldname2[90];
    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_002.hdf", ud.file_name);
    sprintf(fieldname2, "rhs-nodes-post");
    WriteHDF(prhs2file, node->icx, node->icy, node->icz, node->ndim, rhs, fn2, fieldname2);
    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_003.hdf", ud.file_name);
    sprintf(fieldname2, "lap-final");
    WriteHDF(prhs2file, node->icx, node->icy, node->icz, node->ndim, W0, fn2, fieldname2);
#endif

    rhs_max = 0.0;
    for (int jn=node->igy+1; jn<node->icy-node->igy-1; jn++) {
        int mn = jn*node->icx;
        for (int in=node->igx+1; in<node->icx-node->igx-1; in++) {
            int nn = mn+in;
            if (rhs_max < fabs(rhs[nn])) {
                rhs_max = fabs(rhs[nn]);
                imax = in;
                jmax = jn;
            } 
        }
    }
    printf("rhs_max after  = %e\n", rhs_max);
#else
    correction_nodes(Sol, elem, node, (const double**)hplus, hgrav, p2, t, dt);
#endif
    
    for(ii=0; ii<nc; ii++) {
        /* hunt down factor of 0.5 in the next line !! */
        mpv->p2_nodes[ii]  = SCND_PROJ_OLDP_WEIGHT*mpv->p2_nodes[ii] + p_update*SCND_PROJ_DELP_WEIGHT*0.5*p2[ii];
        mpv->dp2_nodes[ii] = p2[ii];
    }

    for(ii=0; ii<nc; ii++) {
        /* hunt down factor of 0.5 in the next line !! */
        mpv->p2_nodes[ii]  = SCND_PROJ_OLDP_WEIGHT*mpv->p2_nodes[ii] + p_update*SCND_PROJ_DELP_WEIGHT*0.5*p2[ii];        
        mpv->dp2_nodes[ii] = p2[ii];
    }
}

/* ========================================================================== */

static double divergence_nodes(
							 double* rhs,
							 const ElemSpaceDiscr* elem,
							 const NodeSpaceDiscr* node,
							 const ConsVars* Sol,
                             double* eta,
							 const MPV* mpv,
							 const BDRY* bdry,
							 const double dt,
							 const double weight) {
	
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
				
            /*
			const int igxn = node->igx;
			const int igyn = node->igy;
			const int icyn = node->icy;
             */
			const int icxn = node->icx;
			
			const int igxe = elem->igx;
			const int icxe = elem->icx;
			const int igye = elem->igy;
			const int icye = elem->icy;
			
			const double dx = node->dx;
			const double dy = node->dy;
			const double oodxdt = 1.0 / (dx * dt);
			const double oodydt = 1.0 / (dy * dt);
			const double oow = 1.0 / 2.0;
            const double todt = 2.0 / dt;
			const double oowdxdt = weight * oow * oodxdt;
			const double oowdydt = weight * oow * oodydt;
			
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
					double tmpfx, tmpfy, tmpfe;
					
					Y = Sol->rhoY[ne] / Sol->rho[ne];
					tmpfx = oowdxdt * Y * Sol->rhou[ne];
					tmpfy = oowdydt * Y * Sol->rhov[ne];
                    tmpfe = todt * eta[ne];
					
					rhs[n]           += + tmpfx + tmpfy - tmpfe;
					rhs[n1]          += - tmpfx + tmpfy + tmpfe;
					rhs[n1icx]       += - tmpfx - tmpfy - tmpfe;
					rhs[nicx]        += + tmpfx - tmpfy + tmpfe;
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
						
			break;
		}
		case 3: {
            
            /*
			const int igxn = node->igx;
			const int igyn = node->igy;
			const int igzn = node->igz;
			const int iczn = node->icz;
             */
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

			const double oodxdt = 1.0 / (dx * dt);
			const double oodydt = 1.0 / (dy * dt);
			const double oodzdt = 1.0 / (dz * dt);
            
			const double oow = 1.0 / 4.0;
			const double oowdxdt = weight * oow * oodxdt;
			const double oowdydt = weight * oow * oodydt;
			const double oowdzdt = weight * oow * oodzdt;
            
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
                        const int nn000  = mn + i;
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
			
			/* set ghost values of rhs to zero 
			j = igyn - 1;
			l = j * icxn;
			for(i = 0; i < icxn; i++) 
				rhs[l + i] = 0.0;
			
			j = icyn - igyn;
			l = j * icxn;
			for(i = 0; i < icxn; i++)
				rhs[l + i] = 0.0;
			
			i = igxn - 1;
			l = i;
			for(j = 0; j < icyn; j++)
				rhs[l + j * icxn] = 0.0;
			
			i = icxn - igxn;
			l = i;
			for(j = 0; j < icyn; j++) 
				rhs[l + j * icxn] = 0.0;
			*/
			break;
		}
		default: ERROR("ndim not in {1, 2, 3}");
	}
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_nodes.hdf");
    sprintf(fieldname, "rhs-nodes");    
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
    
    for (int in = 0; in<node->nc; in++) {
        div_max = MAX_own(div_max, fabs(rhs[in]));
    }
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
    
    const double del = DBL_EPSILON * node->nc;
    
	double tmp = 0.0, rhs_max=0.0;
	
	for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icxicy;
		for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx; 
			for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
				tmp += rhs[n];
                if(fabs(rhs[n])>rhs_max) rhs_max = fabs(rhs[n]);
			}
		}
	}
    
    /* if(fabs(tmp) > 100*del) { */
    if(fabs(tmp) > 10000*del*rhs_max) {
        printf("CHEATING:  tmp = %e,  rhs_max = %e, (tmp/rhs_max)/del = %e\n", tmp/rhs_max, rhs_max, (tmp/rhs_max)/del);
        /* return VIOLATED; */
        return VIOLATED;
    }
    else if(fabs(tmp) > 100*del*rhs_max) {
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
	
	if(y_periodic == 1){
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
	
	if(z_periodic == 1){
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
										double* hgrav,
										const ElemSpaceDiscr* elem,
										const NodeSpaceDiscr* node,
										ConsVars* Sol,
										const ConsVars* Sol0,
										const MPV* mpv,
                                        const double dt) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
    const double Msq = ud.Msq;
    const double Gammainv = th.Gammainv;
    
    const int ndim = node->ndim;
    
    const int impl_grav_th2 = ud.implicit_gravity_theta2;
    const int impl_grav_pr2 = ud.implicit_gravity_press2;
    
    assert(impl_grav_pr2 == 0);
            
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

            const double dx = elem->dx;
            const double dy = elem->dy;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hc      = hcenter;
			double* hg      = hgrav;

			const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);

            const double cexp    = 1.0-th.gamm;
            int i, j, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = 0.0;
						                        
			for(j = igy; j < icy - igy; j++) {
                m = j * icx;
				
                for(i = igx; i < icx - igx; i++) {
                    n = m + i;     
                    
#ifdef THERMCON
                    double theta   = Gammainv * Sol->rhoY[n]*Sol->rhoY[n] / Sol->rho[n] ;
#else
                    double theta   = Sol->rhoY[n] / Sol->rho[n] ;
#endif               
                    double dthetax = (mpv->HydroState->Y0[j] - mpv->HydroState->Y0[j]) / (2.0*dx);
                    double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                    
                    double dthetay = (mpv->HydroState->Y0[j+1] - mpv->HydroState->Y0[j-1]) / (2.0*dy);
                    double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);
                                        
                    hplusx[n]  = theta * gimpx;
                    hplusy[n]  = theta * gimpy;
                    hg[n]      = 0.0;
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
            
            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;

            const int dix = 1;
			const int diy = icx;
			const int diz = icx*icy;

            double Msq = ud.Msq;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hplusz  = hplus[2];
			double* hc      = hcenter;
            
			const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
            
			const double cexp    = 1.0-th.gamm;
			int i, j, k, l, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = hplusz[i] = 0.0;
            
			for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        {
#ifdef THERMCON
                            double theta   = Gammainv * Sol->rhoY[n]*Sol->rhoY[n] / Sol->rho[n] ;
#else
                            double theta   = Sol->rhoY[n] / Sol->rho[n] ;
#endif               
                            
                            double dthetax = ((Sol->rhoY[n+dix]  / Sol->rho[n+dix]) - (Sol->rhoY[n-dix]  / Sol->rho[n-dix])) / (2.0*dx);
                            double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                            
                            double dthetay = ((Sol->rhoY[n+diy]  / Sol->rho[n+diy]) - (Sol->rhoY[n-diy]  / Sol->rho[n-diy])) / (2.0*dy);
                            double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);

                            double dthetaz = ((Sol->rhoY[n+diz]  / Sol->rho[n+diz]) - (Sol->rhoY[n-diz]  / Sol->rho[n-diz])) / (2.0*dz);
                            double gimpz    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetaz/theta);

                            hplusx[n]  = theta * gimpx;
                            hplusy[n]  = theta * gimpy;
                            hplusz[n]  = theta * gimpz;
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


#if 0
/* ========================================================================== */

static void diagonal_preconditioner(
									double* diaginv, 
									const ElemSpaceDiscr* elem,
									const NodeSpaceDiscr* node,
									ConsVars* Sol) {
	
	extern User_Data ud;
	
	const int ndim = node->ndim;
	
	switch(ndim) {
		case 1: {    
			ERROR("interface_enthalpy_nodes() not implemented for 1D\n");
			break;
		}
		case 2: {
			
			int n;
			
#ifdef PRECON
			int i, j, m;
			const int igxn = node->igx;
			const int icxn = node->icx;
			const int igyn = node->igy;
			const int icyn = node->icy;
			
			const int icxe = elem->icx;
			
			for(j = igyn; j < icyn - igyn; j++) {m = j * icxn;
				for(i = igxn; i < icxn - igxn; i++) {n = m + i;
					{
						int nne, nnw, nse, nsw;
						
						nne  =   i   +   j   * icxe;
						nnw  = (i-1) +   j   * icxe;
						nsw  = (i-1) + (j-1) * icxe;
						nse  =   i   + (j-1) * icxe;
						
						/* diaginv[n] = 4.0 / (Sol->rhoY[nne] + Sol->rhoY[nnw] + Sol->rhoY[nsw] + Sol->rhoY[nse]); */

                        diaginv[n] = 4.0 / (  Sol->rhoY[nne] / Sol->rho[nne]
											+ Sol->rhoY[nnw] / Sol->rho[nnw] 
											+ Sol->rhoY[nsw] / Sol->rho[nsw]
											+ Sol->rhoY[nse] / Sol->rho[nse]);
					}
				}
			}
#else
			for(n=0; n<node->nc; n++) {
				diaginv[n] = 1.0;
			}	
#endif
			
			break;
		}
		case 3: {
						
#ifdef PRECON
			int i, j, k, l, m, n;
			const int igxn = node->igx;
			const int icxn = node->icx;
			const int igyn = node->igy;
			const int icyn = node->icy;
			const int igzn = node->igz;
			const int iczn = node->icz;
			
			const int icxe = elem->icx;
			const int icye = elem->icy;
            
            const int dixe = 1;
            const int diye = icxe;
            const int dize = icxe*icye;
			
			for(k = igzn; k < iczn - igzn; k++) {l = k * icxn*icyn;
                for(j = igyn; j < icyn - igyn; j++) {m = l + j * icxn;
                    for(i = igxn; i < icxn - igxn; i++) {n = m + i;
						int nc     = (k-1)*dize + (j-1)*diye + (i-1)*dixe;
						int nc000  = nc;
						int nc010  = nc        + diye;
						int nc011  = nc        + diye + dize;
						int nc001  = nc               + dize;
						int nc100  = nc + dixe;
						int nc110  = nc + dixe + diye;
						int nc111  = nc + dixe + diye + dize;
						int nc101  = nc + dixe        + dize;
						
						/* diaginv[n] = 4.0 / (Sol->rhoY[nne] + Sol->rhoY[nnw] + Sol->rhoY[nsw] + Sol->rhoY[nse]); */
						diaginv[n] = 8.0 / (  Sol->rhoY[nc000] / Sol->rho[nc000]
											+ Sol->rhoY[nc010] / Sol->rho[nc010] 
											+ Sol->rhoY[nc011] / Sol->rho[nc011]
											+ Sol->rhoY[nc001] / Sol->rho[nc001]
                                            + Sol->rhoY[nc100] / Sol->rho[nc100]
											+ Sol->rhoY[nc110] / Sol->rho[nc110] 
											+ Sol->rhoY[nc111] / Sol->rho[nc111]
											+ Sol->rhoY[nc101] / Sol->rho[nc101]);
                    }
                }
			}
#else /* PRECON */
            int n;

			for(n=0; n<node->nc; n++) {
				diaginv[n] = 1.0;
			}	
#endif
			
			break;
		}
		default: ERROR("ndim not in {1,2,3}");
	}
}
#endif

/* ========================================================================== */

void correction_nodes(
					  ConsVars* Sol,
					  const ElemSpaceDiscr* elem,
					  const NodeSpaceDiscr* node,
                      const double* hplus[3], 
                      const double* hgrav,
					  double* p, 
					  const double t,
					  const double dt) {
	
	extern User_Data ud;
	
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
			
			const double dx = node->dx;
			const double dy = node->dy;
			const double oodx = 1.0 / dx;
			const double oody = 1.0 / dy;
			const double dtowdx = 0.5*dt * oodx;
			const double dtowdy = 0.5*dt * oody;
            
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
			
			int i, j, m, me;
			
			for(j = igy; j < icy - igy - 1; j++) {
				m = j * icx; 
				me = j * icxe;
				
				for(i = igx; i < icx - igx - 1; i++) {
					const int n = m + i;
					const int nicx = n + icx;
					const int n1 = n + 1;
					const int n1icx = n1 + icx;
					
					const int ne = me + i; 
					
					const double Dpx   = 0.5 * (p[n1]   - p[n] + p[n1icx] - p[nicx]);
					const double Dpy   = 0.5 * (p[nicx] - p[n] + p[n1icx] - p[n1]);
                    const double thinv = Sol->rho[ne] / Sol->rhoY[ne] ;

					Sol->rhou[ne] += - dtowdx * thinv * hplusx[ne] * Dpx;
					Sol->rhov[ne] += - dtowdy * thinv * hplusy[ne] * Dpy;
				}
			} 
			
			break;
		}
		case 3: {
			
			extern User_Data ud;
			
			const int igxe = elem->igx;
			const int icxe = elem->icx;
			const int igye = elem->igy;
			const int icye = elem->icy;
			const int igze = elem->igz;
			const int icze = elem->icz;
            
			const int dixe = 1;
			const int diye = icxe;
			const int dize = icxe*icye;
			
			const int icxn = node->icx;
			const int icyn = node->icy;
            /*
			const int iczn = node->icz;    
             */
            const int dixn = 1;
			const int diyn = icxn;
			const int dizn = icxn*icyn;
            
			const double dx = node->dx;
			const double dy = node->dy;
			const double dz = node->dz;
			const double oodx = 1.0 / dx;
			const double oody = 1.0 / dy;
			const double oodz = 1.0 / dz;
			const double dtowdx = 0.5*dt * oodx;
			const double dtowdy = 0.5*dt * oody;
			const double dtowdz = 0.5*dt * oodz;
			
			int i, j, k, ln, mn, nn, le, me, ne;
			
			for(k = igze; k < icze - igze; k++) {
				le = k * dize; 
				ln = k * dizn;
								
				for(j = igye; j < icye - igye; j++) {
					me = le + j * diye; 
					mn = ln + j * diyn;
                    
					for(i = igxe; i < icxe - igxe; i++) {
						ne = me + i * dixe;
						nn = mn + i * dixn;
                        
						const int nn000 = nn;
						const int nn010 = nn        + diyn;
						const int nn011 = nn        + diyn + dizn;
						const int nn001 = nn               + dizn;
						const int nn100 = nn + dixn;
						const int nn110 = nn + dixn + diyn;
						const int nn111 = nn + dixn + diyn + dizn;
						const int nn101 = nn + dixn        + dizn;
						
						const double Dpx = 0.25 * (  p[nn100] - p[nn000]
                                                   + p[nn110] - p[nn010]
                                                   + p[nn111] - p[nn011]
                                                   + p[nn101] - p[nn001]
                                                   );
						
						const double Dpy = 0.25 * (  p[nn010] - p[nn000]
                                                   + p[nn110] - p[nn100]
                                                   + p[nn111] - p[nn101]
                                                   + p[nn011] - p[nn001]
                                                   );
						
						const double Dpz = 0.25 * (  p[nn001] - p[nn000]
                                                   + p[nn101] - p[nn100]
						                           + p[nn111] - p[nn110]
						                           + p[nn011] - p[nn010]
                                                   );
						Sol->rhou[ne] += - dtowdx * Dpx;
						Sol->rhov[ne] += - dtowdy * Dpy;
						Sol->rhow[ne] += - dtowdz * Dpz;
					}
				} 
			} 
			break;
		}
		default: ERROR("ndim not in {1, 2,3}");
	}
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
