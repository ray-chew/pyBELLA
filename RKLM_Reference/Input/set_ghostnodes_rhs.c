/*
 *  set_ghostnodes_rhs.c
 *  LowMach.¹
 *
 *  Created by WorkAccount on Sun Feb 22 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#include <assert.h>

#include "error.h"
#include "kgrid.h"
#include "set_ghostnodes_rhs.h"
#include "enum_bdry.h"


/*------------------------------------------------------------------------------
 left p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhs_left_void(
											   double* rhs, 
											   const NodeSpaceDiscr* node,
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhs_left_wall(
											   double* rhs, 
											   const NodeSpaceDiscr* node,
											   const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bobj = igx - 1;
	int i, j, k, l, m, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nobj = bobj - i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhs[m + nobj] = 0.0;
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_left_inflow(
												 double* rhs, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_left_outflow(
												  double* rhs, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_left_periodic(
												   double* rhs, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = icx - igx - 2;
	const int bobj = igx - 1;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc - i;
		nobj = bobj - i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhs[m + nobj] = rhs[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_left_neumann(
												  double* rhs, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhs_left[])(
											   double* rhs, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhs_left_void,
NodeSpaceDiscr_ghost_rhs_left_wall,
NodeSpaceDiscr_ghost_rhs_left_inflow,
NodeSpaceDiscr_ghost_rhs_left_outflow,
NodeSpaceDiscr_ghost_rhs_left_periodic,
NodeSpaceDiscr_ghost_rhs_left_neumann
};


/*------------------------------------------------------------------------------
 right p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhs_right_void(
												double* rhs, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhs_right_wall(
												double* rhs, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bobj = icx - igx;
	int i, j, k, l, m, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i <= ig; i++) {      /* nodes on right included !  "i < ig" replaced */
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhs[m + nobj] = 0.0;
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_right_inflow(
												  double* rhs, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_right_outflow(
												   double* rhs, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_right_periodic(
													double* rhs, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx;           /* nodes on right included !  "igx + 1" replaced */
	const int bobj = icx - igx - 1; /* nodes on right included !  "icx - igx" replaced */
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i <= ig; i++) {      /* nodes on right included !  "i < ig" replaced */
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhs[m + nobj] = rhs[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_right_neumann(
												   double* rhs, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhs_right[])(
												double* rhs, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhs_right_void,
NodeSpaceDiscr_ghost_rhs_right_wall,
NodeSpaceDiscr_ghost_rhs_right_inflow,
NodeSpaceDiscr_ghost_rhs_right_outflow,
NodeSpaceDiscr_ghost_rhs_right_periodic,
NodeSpaceDiscr_ghost_rhs_right_neumann
};

/*------------------------------------------------------------------------------
 bottom p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhs_bottom_void(
												 double* rhs, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhs_bottom_wall(
												 double* rhs, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	
	const int igy = node->igy;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bobj = igy - 1;
	int i, j, k, l, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhs[n + mobj] = 0.0;
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_bottom_inflow(
												   double* rhs, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_bottom_outflow(
													double* rhs, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_bottom_periodic(
													 double* rhs, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = icy - igy - 2;
	const int bobj = igy - 1;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhs[n + mobj] = rhs[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_bottom_neumann(
													double* rhs, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhs_bottom[])(
												 double* rhs, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhs_bottom_void,
NodeSpaceDiscr_ghost_rhs_bottom_wall,
NodeSpaceDiscr_ghost_rhs_bottom_inflow,
NodeSpaceDiscr_ghost_rhs_bottom_outflow,
NodeSpaceDiscr_ghost_rhs_bottom_periodic,
NodeSpaceDiscr_ghost_rhs_bottom_neumann
};

/*------------------------------------------------------------------------------
 top p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhs_top_void(
											  double* rhs, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {}

static void NodeSpaceDiscr_ghost_rhs_top_wall(
											  double* rhs, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {
	
	const int igy = node->igy;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bobj = icy - igy;
	int i, j, k, l, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhs[n + mobj] = 0.0;
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_top_inflow(
												double* rhs, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_top_outflow(
												 double* rhs, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_top_periodic(
												  double* rhs, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy;           /* nodes on top bdry included! "igy + 1" replaced */
	const int bobj = icy - igy - 1; /* nodes on top bdry included! "icy - igy" replaced */
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j <= ig; j++) {      /* nodes on top bdry included! "j < ig" replaced */
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhs[n + mobj] = rhs[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_top_neumann(
												 double* rhs, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhs_top[])(
											  double* rhs, 
											  const NodeSpaceDiscr* node, 
											  const int ig) = {
NodeSpaceDiscr_ghost_rhs_top_void,
NodeSpaceDiscr_ghost_rhs_top_wall,
NodeSpaceDiscr_ghost_rhs_top_inflow,
NodeSpaceDiscr_ghost_rhs_top_outflow,
NodeSpaceDiscr_ghost_rhs_top_periodic,
NodeSpaceDiscr_ghost_rhs_top_neumann
};

/*------------------------------------------------------------------------------
 back p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhs_back_void(
											   double* rhs, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhs_back_wall(
											   double* rhs, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {
	
	const int igz = node->igz;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icxicy = icx * icy;
	
	const int bobj = igz - 1;
	int i, j, k, lobj, m, n;
	
	assert(ig <= igz);
	
	for(k = 0; k < ig; k++) {
		lobj = (bobj - k) * icxicy;
		
		for(j = 0; j < icy; j++) {m = j * icx;
			for(i = 0; i < icx; i++) {n = m + i;
				rhs[n + lobj] = 0.0;
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_back_inflow(
												 double* rhs, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_back_outflow(
												  double* rhs, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_back_periodic(
												   double* rhs, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_back_neumann(
												  double* rhs, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhs_back[])(
											   double* rhs, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhs_back_void,
NodeSpaceDiscr_ghost_rhs_back_wall,
NodeSpaceDiscr_ghost_rhs_back_inflow,
NodeSpaceDiscr_ghost_rhs_back_outflow,
NodeSpaceDiscr_ghost_rhs_back_periodic,
NodeSpaceDiscr_ghost_rhs_back_neumann
};

/*------------------------------------------------------------------------------
 front p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhs_front_void(
												double* rhs, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhs_front_wall(
												double* rhs, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	
	const int igz = node->igz;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int icxicy = icx * icy;
	
	const int bobj = icz - igz;
	int i, j, k, lobj, m, n;
	
	assert(ig <= igz);
	
	for(k = 0; k < ig; k++) {
		lobj = (bobj + k) * icxicy;
		
		for(j = 0; j < icy; j++) {m = j * icx;
			for(i = 0; i < icx; i++) {n = m + i;
				rhs[n + lobj] = 0.0;
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhs_front_inflow(
												  double* rhs, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_front_outflow(
												   double* rhs, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_front_periodic(
													double* rhs, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhs_front_neumann(
												   double* rhs, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhs_front[])(
												double* rhs, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhs_front_void,
NodeSpaceDiscr_ghost_rhs_front_wall,
NodeSpaceDiscr_ghost_rhs_front_inflow,
NodeSpaceDiscr_ghost_rhs_front_outflow,
NodeSpaceDiscr_ghost_rhs_front_periodic,
NodeSpaceDiscr_ghost_rhs_front_neumann
};


void set_ghostnodes_rhs(
						double* rhs, 
						const NodeSpaceDiscr* node, 
						const int ig) {
	
	(*NodeSpaceDiscr_ghost_rhs_left[node->left])(rhs, node, ig);
	(*NodeSpaceDiscr_ghost_rhs_right[node->right])(rhs, node, ig);
	(*NodeSpaceDiscr_ghost_rhs_bottom[node->bottom])(rhs, node, ig);
	(*NodeSpaceDiscr_ghost_rhs_top[node->top])(rhs, node, ig); 
	(*NodeSpaceDiscr_ghost_rhs_back[node->back])(rhs, node, ig);
	(*NodeSpaceDiscr_ghost_rhs_front[node->front])(rhs, node, ig);
}


