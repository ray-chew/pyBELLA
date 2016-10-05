/*
 *  set_ghostnodes_p.c
 *  LowMach.¹
 *
 *  Created by WorkAccount on Sun Feb 22 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#include <assert.h>

#include "error.h"
#include "kgrid.h"
#include "set_ghostnodes_p.h"
#include "enum_bdry.h"


/*------------------------------------------------------------------------------
 left p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_p_left_void(
											 double* p, 
											 const NodeSpaceDiscr* node,
											 const int ig) {}

static void NodeSpaceDiscr_ghost_p_left_wall(
											 double* p, 
											 const NodeSpaceDiscr* node,
											 const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = igx - 1;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj - i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				p[m + nobj] = p[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_left_inflow(
											   double* p, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_left_outflow(
												double* p, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_left_periodic(
												 double* p, 
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
				p[m + nobj] = p[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_left_neumann(
												double* p, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_p_left[])(
											 double* p, 
											 const NodeSpaceDiscr* node, 
											 const int ig) = {
NodeSpaceDiscr_ghost_p_left_void,
NodeSpaceDiscr_ghost_p_left_wall,
NodeSpaceDiscr_ghost_p_left_inflow,
NodeSpaceDiscr_ghost_p_left_outflow,
NodeSpaceDiscr_ghost_p_left_periodic,
NodeSpaceDiscr_ghost_p_left_neumann
};


/*------------------------------------------------------------------------------
 right p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_p_right_void(
											  double* p, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {}

static void NodeSpaceDiscr_ghost_p_right_wall(
											  double* p, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = icx - igx - 2;      
	const int bobj = icx - igx; 
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {      /* What the hell did I want with "i <= ig" ? */
		nsrc = bsrc - i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				p[m + nobj] = p[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_right_inflow(
												double* p, 
												const NodeSpaceDiscr* node, 
												const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_right_outflow(
												 double* p, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_right_periodic(
												  double* p, 
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
				p[m + nobj] = p[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_right_neumann(
												 double* p, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_p_right[])(
											  double* p, 
											  const NodeSpaceDiscr* node, 
											  const int ig) = {
NodeSpaceDiscr_ghost_p_right_void,
NodeSpaceDiscr_ghost_p_right_wall,
NodeSpaceDiscr_ghost_p_right_inflow,
NodeSpaceDiscr_ghost_p_right_outflow,
NodeSpaceDiscr_ghost_p_right_periodic,
NodeSpaceDiscr_ghost_p_right_neumann
};

/*------------------------------------------------------------------------------
 bottom p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_p_bottom_void(
											   double* p, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_p_bottom_wall(
											   double* p, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {
	
	const int igy = node->igy;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_bottom_inflow(
												 double* p, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_bottom_outflow(
												  double* p, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_bottom_periodic(
												   double* p, 
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
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_bottom_neumann(
												  double* p, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_p_bottom[])(
											   double* p, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_p_bottom_void,
NodeSpaceDiscr_ghost_p_bottom_wall,
NodeSpaceDiscr_ghost_p_bottom_inflow,
NodeSpaceDiscr_ghost_p_bottom_outflow,
NodeSpaceDiscr_ghost_p_bottom_periodic,
NodeSpaceDiscr_ghost_p_bottom_neumann
};

/*------------------------------------------------------------------------------
 top p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_p_top_void(
											double* p, 
											const NodeSpaceDiscr* node, 
											const int ig) {}

static void NodeSpaceDiscr_ghost_p_top_wall(
											double* p, 
											const NodeSpaceDiscr* node, 
											const int ig) {
	
	const int igy = node->igy;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = icy - igy - 2;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_top_inflow(
											  double* p, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_top_outflow(
											   double* p, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_top_periodic(
												double* p, 
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
				p[n + mobj] = p[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_p_top_neumann(
											   double* p, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_p_top[])(
											double* p, 
											const NodeSpaceDiscr* node, 
											const int ig) = {
NodeSpaceDiscr_ghost_p_top_void,
NodeSpaceDiscr_ghost_p_top_wall,
NodeSpaceDiscr_ghost_p_top_inflow,
NodeSpaceDiscr_ghost_p_top_outflow,
NodeSpaceDiscr_ghost_p_top_periodic,
NodeSpaceDiscr_ghost_p_top_neumann
};

/*------------------------------------------------------------------------------
 back p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_p_back_void(
											 double* p, 
											 const NodeSpaceDiscr* node, 
											 const int ig) {}

static void NodeSpaceDiscr_ghost_p_back_wall(
											 double* p, 
											 const NodeSpaceDiscr* node, 
											 const int ig) {
	
	const int igz = node->igz;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icxicy = icx * icy;
	
	const int bsrc = igz + 1;
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

static void NodeSpaceDiscr_ghost_p_back_inflow(
											   double* p, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_back_outflow(
												double* p, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_back_periodic(
												 double* p, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_back_neumann(
												double* p, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_p_back[])(
											 double* p, 
											 const NodeSpaceDiscr* node, 
											 const int ig) = {
NodeSpaceDiscr_ghost_p_back_void,
NodeSpaceDiscr_ghost_p_back_wall,
NodeSpaceDiscr_ghost_p_back_inflow,
NodeSpaceDiscr_ghost_p_back_outflow,
NodeSpaceDiscr_ghost_p_back_periodic,
NodeSpaceDiscr_ghost_p_back_neumann
};

/*------------------------------------------------------------------------------
 front p
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_p_front_void(
											  double* p, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {}

static void NodeSpaceDiscr_ghost_p_front_wall(
											  double* p, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {
	
	const int igz = node->igz;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int icxicy = icx * icy;
	
	const int bsrc = icz - igz - 2;
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

static void NodeSpaceDiscr_ghost_p_front_inflow(
												double* p, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_front_outflow(
												 double* p, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_front_periodic(
												  double* p, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_p_front_neumann(
												 double* p, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_p_front[])(
											  double* p, 
											  const NodeSpaceDiscr* node, 
											  const int ig) = {
NodeSpaceDiscr_ghost_p_front_void,
NodeSpaceDiscr_ghost_p_front_wall,
NodeSpaceDiscr_ghost_p_front_inflow,
NodeSpaceDiscr_ghost_p_front_outflow,
NodeSpaceDiscr_ghost_p_front_periodic,
NodeSpaceDiscr_ghost_p_front_neumann
};


void set_ghostnodes_p(
					  double* p, 
					  const NodeSpaceDiscr* node, 
					  const int ig) {
	
	(*NodeSpaceDiscr_ghost_p_left[node->left])(p, node, ig);
	(*NodeSpaceDiscr_ghost_p_right[node->right])(p, node, ig);
	(*NodeSpaceDiscr_ghost_p_bottom[node->bottom])(p, node, ig);
	(*NodeSpaceDiscr_ghost_p_top[node->top])(p, node, ig); 
	(*NodeSpaceDiscr_ghost_p_back[node->back])(p, node, ig);
	(*NodeSpaceDiscr_ghost_p_front[node->front])(p, node, ig);
}


