/*
 *  set_ghostcells_p2_S2.c
 *  LowMach.π
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 * Rupert:  This can be simplified tremendously, as I am never
 using anything but simple zero-slope conditions
 across the domain boundaries 
 */

#include <assert.h>

#include "Common.h"
#include "error.h"
#include "kgrid.h"
#include "set_ghostcells_p.h"
#include "enum_bdry.h"
#include "variable.h"
#include "userdata.h"
#include "thermodynamic.h"
#include "ProjectionType.h"

/*------------------------------------------------------------------------------
 left p
 ------------------------------------------------------------------------------*/
static void set_ghostcells_p2_left_void(
										double* p,
										const ElemSpaceDiscr* elem,
										const int ig) {}

static void set_ghostcells_p2_left_wall(
										double* p,
										const ElemSpaceDiscr* elem,
										const int ig) {
	
	const int igx = elem->igx;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igx;
	const int bobj = igx - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		msrc = (bsrc + i);
		mobj = (bobj - i);
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {n = l + j * icx;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void set_ghostcells_p2_left_inflow(
										  double* p,
										  const ElemSpaceDiscr* elem, 
										  const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_left_outflow(
										   double* p,
										   const ElemSpaceDiscr* elem, 
										   const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_left_periodic(
											double* p,
											const ElemSpaceDiscr* elem, 
											const int ig) {
	
	const int igx = elem->igx;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icx - igx - 1;
	const int bobj = igx - 1;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc - i;
		nobj = bobj - i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				p[m + nobj] = p[m + nsrc];
			}
		}
	}
}

static void set_ghostcells_p2_left_neumann(
										   double* p,
										   const ElemSpaceDiscr* elem, 
										   const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_left_open(
										double* p,
										const ElemSpaceDiscr* elem,
										const int ig) {
	
	/* For the open bc I implement homogeneous Neumann conditions
     as I am making the necessary average modifications to the
	 energy fluxes due to global predictor-divergence in the
	 source term for the Poisson problem. 
	 
	 This approach is ONLY VIABLE for the five-point-stencil.
	 The diagonal effects in the nine-point would invalidate it.
	 Therefore, in interface_coefficiencts_moments() all the 
	 moments need to be set to zero, or the associated contri-
	 butions in the Laplacian and the flux corrections have 
	 to be set identically equal to zero. 
	 */
	
	const int igx = elem->igx;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igx;
	const int bobj = igx - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		msrc = (bsrc + i);
		mobj = (bobj - i);
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {n = l + j * icx;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void (*set_ghostcells_p2_left[])(
										double* p,
										const ElemSpaceDiscr* elem, 
										const int ig) = {
set_ghostcells_p2_left_void,
set_ghostcells_p2_left_wall,
set_ghostcells_p2_left_inflow,
set_ghostcells_p2_left_outflow,
set_ghostcells_p2_left_periodic,
set_ghostcells_p2_left_neumann,
set_ghostcells_p2_left_open
};

/*------------------------------------------------------------------------------
 right p
 ------------------------------------------------------------------------------*/
static void set_ghostcells_p2_right_void(
										 double* p,
										 const ElemSpaceDiscr* elem, 
										 const int ig) {}

static void set_ghostcells_p2_right_wall(
										 double* p,
										 const ElemSpaceDiscr* elem, 
										 const int ig) {
	
	const int igx = elem->igx;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icx - igx - 1;
	const int bobj = icx - igx;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		msrc = (bsrc - i);
		mobj = (bobj + i);
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {n = l + j * icx;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void set_ghostcells_p2_right_inflow(
										   double* p,
										   const ElemSpaceDiscr* elem, 
										   const int ig) {  
	ERROR("function not available");
}

static void set_ghostcells_p2_right_outflow(
											double* p,
											const ElemSpaceDiscr* elem, 
											const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_right_periodic(
											 double* p,
											 const ElemSpaceDiscr* elem, 
											 const int ig) {
	
	const int igx = elem->igx;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igx;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				p[m + nobj] = p[m + nsrc];
			}
		}
	}
}

static void set_ghostcells_p2_right_neumann(
											double* p,
											const ElemSpaceDiscr* elem, 
											const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_right_open(
										 double* p,
										 const ElemSpaceDiscr* elem, 
										 const int ig) {
	
	/* See comment in 
     set_ghostcells_p2_right_open() 
	 above
	 */
	
	const int igx = elem->igx;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icx - igx - 1;
	const int bobj = icx - igx;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		msrc = (bsrc - i);
		mobj = (bobj + i);
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {n = l + j * icx;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void (*set_ghostcells_p2_right[])(
										 double* p,
										 const ElemSpaceDiscr* elem, 
										 const int ig) = {
set_ghostcells_p2_right_void,
set_ghostcells_p2_right_wall,
set_ghostcells_p2_right_inflow,
set_ghostcells_p2_right_outflow,
set_ghostcells_p2_right_periodic,
set_ghostcells_p2_right_neumann,
set_ghostcells_p2_right_open
};

/*------------------------------------------------------------------------------
 bottom p
 ------------------------------------------------------------------------------*/
static void set_ghostcells_p2_bottom_void(
										  double* p,
                                          const double* hplus[3],
                                          const double* hgrav,
										  const ElemSpaceDiscr* elem, 
										  const int ig) {}

static void set_ghostcells_p2_bottom_wall(
										  double* p,
                                          const double* hplus[3],
                                          const double* hgrav,
										  const ElemSpaceDiscr* elem, 
										  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
    const int ify = elem->ify;
    
    const double dy = elem->dy;

	const int bsrc = igy;
	const int bobj = igy - 1;
    
    const double* hplusy = hplus[1];
	
	int i, j, k, l, msrc, mobj, n, nh;
    double xsi;
	
	assert(ig <= igy);
    assert(elem->ndim < 3);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {
            l = k * icx * icy;
			for(i = 0; i < icx; i++) {
                n   = l + i;
                nh  = i*ify + (bsrc + j);
                xsi = 0.5 * dy * hgrav[nh] / hplusy[nh];
				p[n + mobj] = p[n + msrc] * (1.0+xsi) / (1.0-xsi);
			}
		}
	}
}

static void set_ghostcells_p2_bottom_inflow(
											double* p,
                                            const double* hplus[3],
                                            const double* hgrav,
											const ElemSpaceDiscr* elem, 
											const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_bottom_outflow(
											 double* p,
                                             const double* hplus[3],
                                             const double* hgrav,
											 const ElemSpaceDiscr* elem, 
											 const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_bottom_periodic(
											  double* p,
                                              const double* hplus[3],
                                              const double* hgrav,
											  const ElemSpaceDiscr* elem, 
											  const int ig) {
	
	const int icx = elem->icx;
	const int igy = elem->igy;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = igy - 1;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j)*icx;
		mobj = (bobj - j)*icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void set_ghostcells_p2_bottom_neumann(
											 double* p,
                                             const double* hplus[3],
                                             const double* hgrav,
											 const ElemSpaceDiscr* elem, 
											 const int ig) {
	ERROR("function not available");
}

static void (*set_ghostcells_p2_bottom[])(
										  double* p,
                                          const double* hplus[3],
                                          const double* hgrav,
										  const ElemSpaceDiscr* elem, 
										  const int ig) = {
set_ghostcells_p2_bottom_void,
set_ghostcells_p2_bottom_wall,
set_ghostcells_p2_bottom_inflow,
set_ghostcells_p2_bottom_outflow,
set_ghostcells_p2_bottom_periodic,
set_ghostcells_p2_bottom_neumann
};

/*------------------------------------------------------------------------------
 top p
 ------------------------------------------------------------------------------*/
static void set_ghostcells_p2_top_void(
									   double* p,
                                       const double* hplus[3],
                                       const double* hgrav,
									   const ElemSpaceDiscr* elem, 
									   const int ig) {}

static void set_ghostcells_p2_top_wall(
									   double* p,
                                       const double* hplus[3],
                                       const double* hgrav,
									   const ElemSpaceDiscr* elem, 
									   const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
    const int ify = elem->ify;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
    
    const double dy = elem->dy;

    const double* hplusy = hplus[1];

    int i, j, k, l, msrc, mobj, n, nh;
    double xsi;
	
	assert(ig <= igy);
    assert(elem->ndim < 3);

	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {
            l = k * icx * icy;
			for(i = 0; i < icx; i++) {
                n   = l + i;
                nh  = i*ify + (bsrc-j) + 1;
                xsi = 0.5 * dy * hgrav[nh] / hplusy[nh];
				p[n + mobj] = p[n + msrc] * (1.0-xsi) / (1.0+xsi);
			}
		}
	}
}

static void set_ghostcells_p2_top_inflow(
										 double* p,
                                         const double* hplus[3],
                                         const double* hgrav,
										 const ElemSpaceDiscr* elem, 
										 const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_top_outflow(
										  double* p,
                                          const double* hplus[3],
                                          const double* hgrav,
										  const ElemSpaceDiscr* elem, 
										  const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_top_periodic(
										   double* p,
                                           const double* hplus[3],
                                           const double* hgrav,
										   const ElemSpaceDiscr* elem, 
										   const int ig) {
	
	const int icx = elem->icx;
	const int igy = elem->igy;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j)*icx;
		mobj = (bobj + j)*icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void set_ghostcells_p2_top_neumann(
										  double* p,
                                          const double* hplus[3],
                                          const double* hgrav,
										  const ElemSpaceDiscr* elem, 
										  const int ig) {
	ERROR("function not available");
}

static void (*set_ghostcells_p2_top[])(
									   double* p,
                                       const double* hplus[3],
                                       const double* hgrav,
									   const ElemSpaceDiscr* elem, 
									   const int ig) = {
set_ghostcells_p2_top_void,
set_ghostcells_p2_top_wall,
set_ghostcells_p2_top_inflow,
set_ghostcells_p2_top_outflow,
set_ghostcells_p2_top_periodic,
set_ghostcells_p2_top_neumann
};

/*------------------------------------------------------------------------------
 back p
 ------------------------------------------------------------------------------*/
static void set_ghostcells_p2_back_void(
										double* p,
										const ElemSpaceDiscr* elem, 
										const int ig) {}

static void set_ghostcells_p2_back_wall(
										double* p,
										const ElemSpaceDiscr* elem, 
										const int ig) {
	
	const int igz = elem->igz;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icxicy = icx * icy;
	
	const int bsrc = igz;
	const int bobj = igz - 1;
	int i, j, k, lsrc, lobj, m, n;
	
	assert(ig <= igz);
	
	for(k = 0; k < ig; k++) {
		lsrc = (bsrc + k) * icxicy;
		lobj = (bobj - k) * icxicy;
		
		for(j = 0; j < icy; j++) {m = j * icx;
			for(i = 0; i < icx; i++) {n = m + i;
				p[n + lobj] = p[n + lsrc];
			}
		}
	}
}

static void set_ghostcells_p2_back_inflow(
										  double* p,
										  const ElemSpaceDiscr* elem, 
										  const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_back_outflow(
										   double* p,
										   const ElemSpaceDiscr* elem, 
										   const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_back_periodic(
											double* p,
											const ElemSpaceDiscr* elem, 
											const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_back_neumann(
										   double* p,
										   const ElemSpaceDiscr* elem, 
										   const int ig) {
	ERROR("function not available");
}

static void (*set_ghostcells_p2_back[])(
										double* p,
										const ElemSpaceDiscr* elem, 
										const int ig) = {
set_ghostcells_p2_back_void,
set_ghostcells_p2_back_wall,
set_ghostcells_p2_back_inflow,
set_ghostcells_p2_back_outflow,
set_ghostcells_p2_back_periodic,
set_ghostcells_p2_back_neumann
};

/*------------------------------------------------------------------------------
 front p
 ------------------------------------------------------------------------------*/
static void set_ghostcells_p2_front_void(
										 double* p,
										 const ElemSpaceDiscr* elem, 
										 const int ig) {}

static void set_ghostcells_p2_front_wall(
										 double* p,
										 const ElemSpaceDiscr* elem, 
										 const int ig) {
	
	const int igz = elem->igz;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int icxicy = icx * icy;
	
	const int bsrc = icz - igz - 1;
	const int bobj = icz - igz;
	int i, j, k, lsrc, lobj, m, n;
	
	assert(ig <= igz);
	
	for(k = 0; k < ig; k++) {
		lsrc = (bsrc - k) * icxicy;
		lobj = (bobj + k) * icxicy;
		
		for(j = 0; j < icy; j++) {m = j * icx;
			for(i = 0; i < icx; i++) {n = m + i;
				p[n + lobj] = p[n + lsrc];
			}
		}
	}
}

static void set_ghostcells_p2_front_inflow(
										   double* p,
										   const ElemSpaceDiscr* elem, 
										   const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_front_outflow(
											double* p,
											const ElemSpaceDiscr* elem, 
											const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_front_periodic(
											 double* p,
											 const ElemSpaceDiscr* elem, 
											 const int ig) {
	ERROR("function not available");
}

static void set_ghostcells_p2_front_neumann(
											double* p,
											const ElemSpaceDiscr* elem, 
											const int ig) {
	ERROR("function not available");
}

static void (*set_ghostcells_p2_front[])(
										 double* p,
										 const ElemSpaceDiscr* elem, 
										 const int ig) = {
set_ghostcells_p2_front_void,
set_ghostcells_p2_front_wall,
set_ghostcells_p2_front_inflow,
set_ghostcells_p2_front_outflow,
set_ghostcells_p2_front_periodic,
set_ghostcells_p2_front_neumann
};


void set_ghostcells_p2(
					   double* p,
                       const double* hplus[3],
                       const double* hgrav,
					   const ElemSpaceDiscr* elem, 
					   const int ig) {
	
	
	(*set_ghostcells_p2_left[elem->left])(p, elem, ig);
	(*set_ghostcells_p2_right[elem->right])(p, elem, ig);
	(*set_ghostcells_p2_bottom[elem->bottom])(p, hplus, hgrav, elem, ig);
	(*set_ghostcells_p2_top[elem->top])(p, hplus, hgrav, elem, ig); 
	(*set_ghostcells_p2_back[elem->back])(p, elem, ig);
	(*set_ghostcells_p2_front[elem->front])(p, elem, ig);
}
