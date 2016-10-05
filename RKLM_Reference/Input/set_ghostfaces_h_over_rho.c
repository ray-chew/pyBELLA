/*
 *  set_ghostfaces_h_over_rho.c
 *  LowMach.π
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <assert.h>

#include "error.h"
#include "kgrid.h"
#include "set_ghostfaces_h_over_rho.h"
#include "enum_bdry.h"



/*------------------------------------------------------------------------------
 left p
 ------------------------------------------------------------------------------*/
static void set_ghostfaces_h_over_rho_left_void(
												double** hplus, 
												double** hminus, 
												double** lhplus, 
												double** lhminus,  
												const NodeSpaceDiscr* node,
												const int ig) {}

static void set_ghostfaces_h_over_rho_left_wall(
												double** hplus, 
												double** hminus, 
												double** lhplus, 
												double** lhminus,  
												const NodeSpaceDiscr* node,
												const int ig) {
	
	double* hplusx   = hplus[0];
	double* hminusx  = hminus[0];
	double* lhplusx  = lhplus[0];
	double* lhminusx = lhminus[0];
	
	double* hplusy   = hplus[1];
	double* hminusy  = hminus[1];
	double* lhplusy  = lhplus[1];
	double* lhminusy = lhminus[1];
	
	const int igx = node->igx;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	
	const int ifx = node->ifx;
	const int ify = node->ify;
	
	const int bsrc_face_x = igx + 1;
	const int bobj_face_x = igx;
	
	const int bsrc_node_y = igx * ify;
	const int bobj_node_y = igx * ify;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc_face_x + i;
		nobj = bobj_face_x - i;
		for(k = 0; k < icz; k++) {l = k * ifx * icy;
			for(j = 0; j < icy; j++) {m = l + j * ifx;
				hplusx[m + nobj]   = 0.0;
				hminusx[m + nobj]  = 0.0;
				lhplusx[m + nobj]  = 0.0;
				lhminusx[m + nobj] = 0.0;
			}
		}
	}
	
    i = 0;
    nsrc = bsrc_node_y + i*ify;
    nobj = bobj_node_y - i*ify;
    for(k = 0; k < icz; k++) {l = k * icx * ify;
		for(j = 0; j < ify; j++) {m = l + j;
			hminusy[m + nobj]  = 0.0;
			lhminusy[m + nobj] = 0.0;
		}
    }
	
	
	for(i = 1; i < ig; i++) {
		nsrc = bsrc_node_y + i*ify;
		nobj = bobj_node_y - i*ify;
		for(k = 0; k < icz; k++) {l = k * icx * ify;
			for(j = 0; j < ify; j++) {m = l + j;
				hplusy[m + nobj]   = 0.0;
				lhplusy[m + nobj]  = 0.0;
				hminusy[m + nobj]  = 0.0;
				lhminusy[m + nobj] = 0.0;
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_left_inflow(
												  double** hplus, 
												  double** hminus, 
												  double** lhplus, 
												  double** lhminus,  
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_left_outflow(
												   double** hplus, 
												   double** hminus, 
												   double** lhplus, 
												   double** lhminus,  
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_left_periodic(
													double** hplus, 
													double** hminus, 
													double** lhplus, 
													double** lhminus,  
													const NodeSpaceDiscr* node, 
													const int ig) {
	
	double* hplusx   = hplus[0];
	double* hminusx  = hminus[0];
	double* lhplusx  = lhplus[0];
	double* lhminusx = lhminus[0];
	
	double* hminusy  = hminus[1];
	double* lhminusy = lhminus[1];
	
	const int igx = node->igx;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	
	const int ifx = node->ifx;
	const int ify = node->ify;
	
	const int bsrc_face_x = ifx - igx - 2;
	const int bobj_face_x = igx;
	
	const int bsrc_node_y = (icx - igx - 1) * ify;
	const int bobj_node_y = igx * ify;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc_face_x - i;
		nobj = bobj_face_x - i;
		for(k = 0; k < icz; k++) {l = k * ifx * icy;
			for(j = 0; j < icy; j++) {m = l + j * ifx;
				hplusx[m + nobj]   = hplusx[m + nsrc];
				hminusx[m + nobj]  = hminusx[m + nsrc];
				lhplusx[m + nobj]  = lhplusx[m + nsrc];
				lhminusx[m + nobj] = lhminusx[m + nsrc];
			}
		}
	}
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc_node_y - i*ify;
		nobj = bobj_node_y - i*ify;
		for(k = 0; k < icz; k++) {l = k * icx * ify;
			for(j = 0; j < ify; j++) {m = l + j;
				hminusy[m + nobj]  = hminusy[m + nsrc];
				lhminusy[m + nobj] = lhminusy[m + nsrc];
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_left_neumann(
												   double** hplus, 
												   double** hminus, 
												   double** lhplus, 
												   double** lhminus,  
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*set_ghostfaces_h_over_rho_left[])(
												double** hplus, 
												double** hminus, 
												double** lhplus, 
												double** lhminus,  
												const NodeSpaceDiscr* node, 
												const int ig) = {
set_ghostfaces_h_over_rho_left_void,
set_ghostfaces_h_over_rho_left_wall,
set_ghostfaces_h_over_rho_left_inflow,
set_ghostfaces_h_over_rho_left_outflow,
set_ghostfaces_h_over_rho_left_periodic,
set_ghostfaces_h_over_rho_left_neumann
};

/*------------------------------------------------------------------------------
 right p
 ------------------------------------------------------------------------------*/
static void set_ghostfaces_h_over_rho_right_void(
												 double** hplus, 
												 double** hminus, 
												 double** lhplus, 
												 double** lhminus,  
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void set_ghostfaces_h_over_rho_right_wall(
												 double** hplus, 
												 double** hminus, 
												 double** lhplus, 
												 double** lhminus,  
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	
	double* hplusx   = hplus[0];
	double* hminusx  = hminus[0];
	double* lhplusx  = lhplus[0];
	double* lhminusx = lhminus[0];
	
	double* hplusy   = hplus[1];
	double* hminusy  = hminus[1];
	double* lhplusy  = lhplus[1];
	double* lhminusy = lhminus[1];
	
	const int igx = node->igx;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	
	const int ifx = node->ifx;
	const int ify = node->ify;
	
	const int bsrc_face_x = ifx - igx - 2;
	const int bobj_face_x = ifx - igx - 1;
	
	const int bsrc_node_y = (icx - igx - 1) * ify;
	const int bobj_node_y = (icx - igx - 1) * ify;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc_face_x - i;
		nobj = bobj_face_x + i;
		for(k = 0; k < icz; k++) {l = k * ifx * icy;
			for(j = 0; j < icy; j++) {m = l + j * ifx;
				hplusx[m + nobj]   = 0.0;
				hminusx[m + nobj]  = 0.0;
				lhplusx[m + nobj]  = 0.0;
				lhminusx[m + nobj] = 0.0;
			}
		}
	}
	
    i = 0; 
    nsrc = bsrc_node_y - i*ify;
    nobj = bobj_node_y + i*ify;
    for(k = 0; k < icz; k++) {l = k * icx * ify;
		for(j = 0; j < icy; j++) {m = l + j;
			hplusy[m + nobj]   = 0.0;
			lhplusy[m + nobj]  = 0.0;
		}
    }
	
	
	for(i = 1; i < ig; i++) {
		nsrc = bsrc_node_y - i*ify;
		nobj = bobj_node_y + i*ify;
		for(k = 0; k < icz; k++) {l = k * icx * ify;
			for(j = 0; j < icy; j++) {m = l + j;
				hplusy[m + nobj]    = 0.0;
				lhplusy[m + nobj]   = 0.0;
				hminusy[m + nobj]   = 0.0;
				lhminusy[m + nobj]  = 0.0;
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_right_inflow(
												   double** hplus, 
												   double** hminus, 
												   double** lhplus, 
												   double** lhminus,  
												   const NodeSpaceDiscr* node, 
												   const int ig) {  
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_right_outflow(
													double** hplus, 
													double** hminus, 
													double** lhplus, 
													double** lhminus,  
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_right_periodic(
													 double** hplus, 
													 double** hminus, 
													 double** lhplus, 
													 double** lhminus,  
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	double* hplusx   = hplus[0];
	double* hminusx  = hminus[0];
	double* lhplusx  = lhplus[0];
	double* lhminusx = lhminus[0];
	
	double* hplusy   = hplus[1];
	double* lhplusy  = lhplus[1];
	
	const int igx = node->igx;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	
	const int ifx = node->ifx;
	const int ify = node->ify;
	
	const int bsrc_face_x = igx + 1;
	const int bobj_face_x = ifx - igx - 1;
	
	const int bsrc_node_y = igx * ify;
	const int bobj_node_y = (icx - igx - 1) * ify;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc_face_x + i;
		nobj = bobj_face_x + i;
		for(k = 0; k < icz; k++) {l = k * ifx * icy;
			for(j = 0; j < icy; j++) {m = l + j * ifx;
				hplusx[m + nobj]   = hplusx[m + nsrc];
				hminusx[m + nobj]  = hminusx[m + nsrc];
				lhplusx[m + nobj]  = lhplusx[m + nsrc];
				lhminusx[m + nobj] = lhminusx[m + nsrc];
			}
		}
	}
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc_node_y + i*ify;
		nobj = bobj_node_y + i*ify;
		for(k = 0; k < icz; k++) {l = k * icx * ify;
			for(j = 0; j < icy; j++) {m = l + j;
				hplusy[m + nobj]   = hplusy[m + nsrc];
				lhplusy[m + nobj]  = lhplusy[m + nsrc];
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_right_neumann(
													double** hplus, 
													double** hminus, 
													double** lhplus, 
													double** lhminus,  
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*set_ghostfaces_h_over_rho_right[])(
												 double** hplus, 
												 double** hminus, 
												 double** lhplus, 
												 double** lhminus,  
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
set_ghostfaces_h_over_rho_right_void,
set_ghostfaces_h_over_rho_right_wall,
set_ghostfaces_h_over_rho_right_inflow,
set_ghostfaces_h_over_rho_right_outflow,
set_ghostfaces_h_over_rho_right_periodic,
set_ghostfaces_h_over_rho_right_neumann
};

/*------------------------------------------------------------------------------
 bottom p
 ------------------------------------------------------------------------------*/
static void set_ghostfaces_h_over_rho_bottom_void(
												  double** hplus, 
												  double** hminus, 
												  double** lhplus, 
												  double** lhminus,  
												  const NodeSpaceDiscr* node, 
												  const int ig) {}

static void set_ghostfaces_h_over_rho_bottom_wall(
												  double** hplus, 
												  double** hminus, 
												  double** lhplus, 
												  double** lhminus,  
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	
	double* hplusx   = hplus[0];
	double* hminusx  = hminus[0];
	double* lhplusx  = lhplus[0];
	double* lhminusx = lhminus[0];
	
	double* hplusy   = hplus[1];
	double* hminusy  = hminus[1];
	double* lhplusy  = lhplus[1];
	double* lhminusy = lhminus[1];
	
	const int igy = node->igy;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	
	const int ifx = node->ifx;
	const int ify = node->ify;
	
	const int bsrc_face_y = igy + 1;
	const int bobj_face_y = igy;
	
	const int bsrc_node_x = igy * ifx;
	const int bobj_node_x = igy * ifx;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		nsrc = bsrc_face_y + j;
		nobj = bobj_face_y - j;
		for(k = 0; k < icz; k++) {l = k * ify * icx;
			for(i = 0; i < icx; i++) {m = l + i * ify;
				hplusy[m + nobj]   = 0.0;
				hminusy[m + nobj]  = 0.0;
				lhplusy[m + nobj]  = 0.0;
				lhminusy[m + nobj] = 0.0;
			}
		}
	}
	
    j = 0;
    nsrc = bsrc_node_x + j*ifx;
    nobj = bobj_node_x - j*ifx;
    for(k = 0; k < icz; k++) {l = k * icy * ifx;
		for(i = 0; i < ifx; i++) {m = l + i;
			hminusx[m + nobj]  = 0.0;
			lhminusx[m + nobj] = 0.0;
		}
    }
	
	
	for(j = 1; j < ig; j++) {
		nsrc = bsrc_node_x + j*ifx;
		nobj = bobj_node_x - j*ifx;
		for(k = 0; k < icz; k++) {l = k * icy * ifx;
			for(i = 0; i < ifx; i++) {m = l + i;
				hplusx[m + nobj]   = 0.0;
				lhplusx[m + nobj]  = 0.0;
				hminusx[m + nobj]  = 0.0;
				lhminusx[m + nobj] = 0.0;
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_bottom_inflow(
													double** hplus, 
													double** hminus, 
													double** lhplus, 
													double** lhminus,  
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_bottom_outflow(
													 double** hplus, 
													 double** hminus, 
													 double** lhplus, 
													 double** lhminus,  
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_bottom_periodic(
													  double** hplus, 
													  double** hminus, 
													  double** lhplus, 
													  double** lhminus,  
													  const NodeSpaceDiscr* node, 
													  const int ig) {
	
	double* hminusx  = hminus[0];
	double* lhminusx = lhminus[0];
	
	double* hplusy   = hplus[1];
	double* hminusy  = hminus[1];
	double* lhplusy  = lhplus[1];
	double* lhminusy = lhminus[1];
	
	const int igy = node->igy;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	
	const int ifx = node->ifx;
	const int ify = node->ify;
	
	const int bsrc_face_y = ify - igy - 2;
	const int bobj_face_y = igy;
	
	const int bsrc_node_x = (icy - igy - 1) * ifx;
	const int bobj_node_x = igy * ifx;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		nsrc = bsrc_face_y - j;
		nobj = bobj_face_y - j;
		for(k = 0; k < icz; k++) {l = k * ify * icx;
			for(i = 0; i < icx; i++) {m = l + i * ify;
				hplusy[m + nobj]   = hplusy[m + nsrc];
				hminusy[m + nobj]  = hminusy[m + nsrc];
				lhplusy[m + nobj]  = lhplusy[m + nsrc];
				lhminusy[m + nobj] = lhminusy[m + nsrc];
			}
		}
	}
	
	for(j = 0; j < ig; j++) {
		nsrc = bsrc_node_x - j*ifx;
		nobj = bobj_node_x - j*ifx;
		for(k = 0; k < icz; k++) {l = k * icy * ifx;
			for(i = 0; i < ifx; i++) {m = l + i;
				hminusx[m + nobj]  = hminusx[m + nsrc];
				lhminusx[m + nobj] = lhminusx[m + nsrc];
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_bottom_neumann(
													 double** hplus, 
													 double** hminus, 
													 double** lhplus, 
													 double** lhminus,  
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void (*set_ghostfaces_h_over_rho_bottom[])(
												  double** hplus, 
												  double** hminus, 
												  double** lhplus, 
												  double** lhminus,  
												  const NodeSpaceDiscr* node, 
												  const int ig) = {
set_ghostfaces_h_over_rho_bottom_void,
set_ghostfaces_h_over_rho_bottom_wall,
set_ghostfaces_h_over_rho_bottom_inflow,
set_ghostfaces_h_over_rho_bottom_outflow,
set_ghostfaces_h_over_rho_bottom_periodic,
set_ghostfaces_h_over_rho_bottom_neumann
};

/*------------------------------------------------------------------------------
 top p
 ------------------------------------------------------------------------------*/
static void set_ghostfaces_h_over_rho_top_void(
											   double** hplus, 
											   double** hminus, 
											   double** lhplus, 
											   double** lhminus,  
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void set_ghostfaces_h_over_rho_top_wall(
											   double** hplus, 
											   double** hminus, 
											   double** lhplus, 
											   double** lhminus,  
											   const NodeSpaceDiscr* node, 
											   const int ig) {
	
	double* hplusx   = hplus[0];
	double* hminusx  = hminus[0];
	double* lhplusx  = lhplus[0];
	double* lhminusx = lhminus[0];
	
	double* hplusy   = hplus[1];
	double* hminusy  = hminus[1];
	double* lhplusy  = lhplus[1];
	double* lhminusy = lhminus[1];
	
	const int igy = node->igy;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int ifx = node->ifx;
	const int ify = node->ify;
    
	const int bsrc_face_y = ify - igy - 2;
	const int bobj_face_y = ify - igy - 1;
	
	const int bsrc_node_x = (icy - igy - 1) * ifx;
	const int bobj_node_x = (icy - igy - 1) * ifx;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		nsrc = bsrc_face_y - j;
		nobj = bobj_face_y + j;
		for(k = 0; k < icz; k++) {l = k * ify * icx;
			for(i = 0; i < icx; i++) {m = l + i * ify;
				hplusy[m + nobj]   = 0.0;
				hminusy[m + nobj]  = 0.0;
				lhplusy[m + nobj]  = 0.0;
				lhminusy[m + nobj] = 0.0;
			}
		}
	}
	
    j = 0;
    nsrc = bsrc_node_x - j*ifx;
    nobj = bobj_node_x + j*ifx;
    for(k = 0; k < icz; k++) {l = k * icy * ifx;
		for(i = 0; i < icx; i++) {m = l + i;
			hplusx[m + nobj]   = 0.0;
			lhplusx[m + nobj]  = 0.0;
		}
    }
	
	
	for(j = 1; j < ig; j++) {
		nsrc = bsrc_node_x - j*ifx;
		nobj = bobj_node_x + j*ifx;
		for(k = 0; k < icz; k++) {l = k * icy * ifx;
			for(i = 0; i < icx; i++) {m = l + i;
				hminusx[m + nobj]   = 0.0;
				lhminusx[m + nobj]  = 0.0;
				hplusx[m + nobj]    = 0.0;
				lhplusx[m + nobj]   = 0.0;
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_top_inflow(
												 double** hplus, 
												 double** hminus, 
												 double** lhplus, 
												 double** lhminus,  
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_top_outflow(
												  double** hplus, 
												  double** hminus, 
												  double** lhplus, 
												  double** lhminus,  
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_top_periodic(
												   double** hplus, 
												   double** hminus, 
												   double** lhplus, 
												   double** lhminus,  
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	double* hplusx   = hplus[0];
	double* lhplusx  = lhplus[0];
	
	double* hplusy   = hplus[1];
	double* hminusy  = hminus[1];
	double* lhplusy  = lhplus[1];
	double* lhminusy = lhminus[1];
	
	const int igy = node->igy;
	
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int ifx = node->ifx;
	const int ify = node->ify;
    
	const int bsrc_face_y = igy + 1;
	const int bobj_face_y = ify - igy - 1;
	
	const int bsrc_node_x = igy * ifx;
	const int bobj_node_x = (icy - igy - 1) * ifx;
	
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		nsrc = bsrc_face_y + j;
		nobj = bobj_face_y + j;
		for(k = 0; k < icz; k++) {l = k * ify * icx;
			for(i = 0; i < icx; i++) {m = l + i * ify;
				hplusy[m + nobj]   = hplusy[m + nsrc];
				hminusy[m + nobj]  = hminusy[m + nsrc];
				lhplusy[m + nobj]  = lhplusy[m + nsrc];
				lhminusy[m + nobj] = lhminusy[m + nsrc];
			}
		}
	}
	
	for(j = 0; j < ig; j++) {
		nsrc = bsrc_node_x + j*ifx;
		nobj = bobj_node_x + j*ifx;
		for(k = 0; k < icz; k++) {l = k * icy * ifx;
			for(i = 0; i < icx; i++) {m = l + i;
				hplusx[m + nobj]   = hplusx[m + nsrc];
				lhplusx[m + nobj]  = lhplusx[m + nsrc];
			}
		}
	}
}

static void set_ghostfaces_h_over_rho_top_neumann(
												  double** hplus, 
												  double** hminus, 
												  double** lhplus, 
												  double** lhminus,  
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*set_ghostfaces_h_over_rho_top[])(
											   double** hplus, 
											   double** hminus, 
											   double** lhplus, 
											   double** lhminus,  
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
set_ghostfaces_h_over_rho_top_void,
set_ghostfaces_h_over_rho_top_wall,
set_ghostfaces_h_over_rho_top_inflow,
set_ghostfaces_h_over_rho_top_outflow,
set_ghostfaces_h_over_rho_top_periodic,
set_ghostfaces_h_over_rho_top_neumann
};

/*------------------------------------------------------------------------------
 back p
 ------------------------------------------------------------------------------*/
static void set_ghostfaces_h_over_rho_back_void(
												double** hplus, 
												double** hminus, 
												double** lhplus, 
												double** lhminus,  
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void set_ghostfaces_h_over_rho_back_wall(
												double** hplus, 
												double** hminus, 
												double** lhplus, 
												double** lhminus,  
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_back_inflow(
												  double** hplus, 
												  double** hminus, 
												  double** lhplus, 
												  double** lhminus,  
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_back_outflow(
												   double** hplus, 
												   double** hminus, 
												   double** lhplus, 
												   double** lhminus,  
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_back_periodic(
													double** hplus, 
													double** hminus, 
													double** lhplus, 
													double** lhminus,  
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_back_neumann(
												   double** hplus, 
												   double** hminus, 
												   double** lhplus, 
												   double** lhminus,  
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*set_ghostfaces_h_over_rho_back[])(
												double** hplus, 
												double** hminus, 
												double** lhplus, 
												double** lhminus,  
												const NodeSpaceDiscr* node, 
												const int ig) = {
set_ghostfaces_h_over_rho_back_void,
set_ghostfaces_h_over_rho_back_wall,
set_ghostfaces_h_over_rho_back_inflow,
set_ghostfaces_h_over_rho_back_outflow,
set_ghostfaces_h_over_rho_back_periodic,
set_ghostfaces_h_over_rho_back_neumann
};

/*------------------------------------------------------------------------------
 front p
 ------------------------------------------------------------------------------*/
static void set_ghostfaces_h_over_rho_front_void(
												 double** hplus, 
												 double** hminus, 
												 double** lhplus, 
												 double** lhminus,  
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void set_ghostfaces_h_over_rho_front_wall(
												 double** hplus, 
												 double** hminus, 
												 double** lhplus, 
												 double** lhminus,  
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_front_inflow(
												   double** hplus, 
												   double** hminus, 
												   double** lhplus, 
												   double** lhminus,  
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_front_outflow(
													double** hplus, 
													double** hminus, 
													double** lhplus, 
													double** lhminus,  
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_front_periodic(
													 double** hplus, 
													 double** hminus, 
													 double** lhplus, 
													 double** lhminus,  
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void set_ghostfaces_h_over_rho_front_neumann(
													double** hplus, 
													double** hminus, 
													double** lhplus, 
													double** lhminus,  
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*set_ghostfaces_h_over_rho_front[])(
												 double** hplus, 
												 double** hminus, 
												 double** lhplus, 
												 double** lhminus,  
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
set_ghostfaces_h_over_rho_front_void,
set_ghostfaces_h_over_rho_front_wall,
set_ghostfaces_h_over_rho_front_inflow,
set_ghostfaces_h_over_rho_front_outflow,
set_ghostfaces_h_over_rho_front_periodic,
set_ghostfaces_h_over_rho_front_neumann
};


void set_ghostfaces_h_over_rho(
							   double** hplus, 
							   double** hminus, 
							   double** lhplus, 
							   double** lhminus, 
							   const NodeSpaceDiscr* node, 
							   const int ig) {
	
	(*set_ghostfaces_h_over_rho_left[node->left])(hplus, hminus, lhplus, lhminus, node, ig);
	(*set_ghostfaces_h_over_rho_right[node->right])(hplus, hminus, lhplus, lhminus, node, ig);
	(*set_ghostfaces_h_over_rho_bottom[node->bottom])(hplus, hminus, lhplus, lhminus, node, ig);
	(*set_ghostfaces_h_over_rho_top[node->top])(hplus, hminus, lhplus, lhminus, node, ig); 
	(*set_ghostfaces_h_over_rho_back[node->back])(hplus, hminus, lhplus, lhminus, node, ig);
	(*set_ghostfaces_h_over_rho_front[node->front])(hplus, hminus, lhplus, lhminus, node, ig);
}

