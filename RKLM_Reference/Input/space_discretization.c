/*******************************************************************************
 File:   space_discretization.c
 Author: Nicola
 Date:   Fri Feb 27 20:28:55 CET 1998
 *******************************************************************************/
#include "Common.h"
#include "kgrid.h"
#include "space_discretization.h"
#include "enum_bdry.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "error.h"
#include "variable.h"


/*------------------------------------------------------------------------------
 left rho
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rho_left_void(
											   double* rho, 
											   const ElemSpaceDiscr* elem,
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rho_left_wall(
											   double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_left_inflow(
												 double* rho, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_left_outflow(
												  double* rho, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_left_periodic(
												   double* rho, 
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
				rho[m + nobj] = rho[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_left_neumann(
												  double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rho_left[])(
											   double* rho, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rho_left_void,
ElemSpaceDiscr_ghost_rho_left_wall,
ElemSpaceDiscr_ghost_rho_left_inflow,
ElemSpaceDiscr_ghost_rho_left_outflow,
ElemSpaceDiscr_ghost_rho_left_periodic,
ElemSpaceDiscr_ghost_rho_left_neumann,
ElemSpaceDiscr_ghost_rho_left_neumann
};


/*------------------------------------------------------------------------------
 left rhou
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhou_left_void(
												double* rhou, 
												const ElemSpaceDiscr* elem,
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhou_left_wall(
												double* rhou, 
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
				rhou[n + mobj] = -rhou[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_left_inflow(
												  double* rhou, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_left_outflow(
												   double* rhou, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_left_periodic(
													double* rhou, 
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
				rhou[m + nobj] = rhou[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_left_neumann(
												   double* rhou, 
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
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhou_left[])(
												double* rhou, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhou_left_void,
ElemSpaceDiscr_ghost_rhou_left_wall,
ElemSpaceDiscr_ghost_rhou_left_inflow,
ElemSpaceDiscr_ghost_rhou_left_outflow,
ElemSpaceDiscr_ghost_rhou_left_periodic,
ElemSpaceDiscr_ghost_rhou_left_neumann,
ElemSpaceDiscr_ghost_rhou_left_neumann
};


/*------------------------------------------------------------------------------
 left rhov
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhov_left_void(
												double* rhov, 
												const ElemSpaceDiscr* elem,
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhov_left_wall(
												double* rhov, 
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
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_left_inflow(
												  double* rhov, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_left_outflow(
												   double* rhov, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_left_periodic(
													double* rhov, 
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
				rhov[m + nobj] = rhov[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_left_neumann(
												   double* rhov, 
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
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhov_left[])(
												double* rhov, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhov_left_void,
ElemSpaceDiscr_ghost_rhov_left_wall,
ElemSpaceDiscr_ghost_rhov_left_inflow,
ElemSpaceDiscr_ghost_rhov_left_outflow,
ElemSpaceDiscr_ghost_rhov_left_periodic,
ElemSpaceDiscr_ghost_rhov_left_neumann,
ElemSpaceDiscr_ghost_rhov_left_neumann
};


/*------------------------------------------------------------------------------
 left rhow
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhow_left_void(
												double* rhow, 
												const ElemSpaceDiscr* elem,
												const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_left_wall(
												double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_left_inflow(
												  double* rhow, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_left_outflow(
												   double* rhow, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_left_periodic(
													double* rhow, 
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
				rhow[m + nobj] = rhow[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_left_neumann(
												   double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhow_left[])(
												double* rhow, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhow_left_void,
ElemSpaceDiscr_ghost_rhow_left_wall,
ElemSpaceDiscr_ghost_rhow_left_inflow,
ElemSpaceDiscr_ghost_rhow_left_outflow,
ElemSpaceDiscr_ghost_rhow_left_periodic,
ElemSpaceDiscr_ghost_rhow_left_neumann,
ElemSpaceDiscr_ghost_rhow_left_neumann
};


/*------------------------------------------------------------------------------
 left rhoe
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoe_left_void(
												double* rhoe, 
												const ElemSpaceDiscr* elem,
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhoe_left_wall(
												double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_left_inflow(
												  double* rhoe, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_left_outflow(
												   double* rhoe, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_left_periodic(
													double* rhoe, 
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
				rhoe[m + nobj] = rhoe[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_left_neumann(
												   double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhoe_left[])(
												double* rhoe, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhoe_left_void,
ElemSpaceDiscr_ghost_rhoe_left_wall,
ElemSpaceDiscr_ghost_rhoe_left_inflow,
ElemSpaceDiscr_ghost_rhoe_left_outflow,
ElemSpaceDiscr_ghost_rhoe_left_periodic,
ElemSpaceDiscr_ghost_rhoe_left_neumann,
ElemSpaceDiscr_ghost_rhoe_left_neumann
};


/*------------------------------------------------------------------------------
 left rhoY
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoY_left_void(
												double* rhoY, 
												const ElemSpaceDiscr* elem,
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhoY_left_wall(
												double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_left_inflow(
												  double* rhoY, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_left_outflow(
												   double* rhoY, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_left_periodic(
													double* rhoY, 
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
				rhoY[m + nobj] = rhoY[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_left_neumann(
												   double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhoY_left[])(
												double* rhoY, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhoY_left_void,
ElemSpaceDiscr_ghost_rhoY_left_wall,
ElemSpaceDiscr_ghost_rhoY_left_inflow,
ElemSpaceDiscr_ghost_rhoY_left_outflow,
ElemSpaceDiscr_ghost_rhoY_left_periodic,
ElemSpaceDiscr_ghost_rhoY_left_neumann,
ElemSpaceDiscr_ghost_rhoY_left_neumann
};


/*------------------------------------------------------------------------------
 left rhoZ
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoZ_left_void(
												double** rhoZ, 
												const ElemSpaceDiscr* elem,
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhoZ_left_wall(
												double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_left_inflow(
												  double** rhoZ, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_left_outflow(
												   double** rhoZ, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_left_periodic(
													double** rhoZ, 
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
				rhoZ[PRES][m + nobj] = rhoZ[PRES][m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_left_neumann(
												   double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhoZ_left[])(
												double** rhoZ, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhoZ_left_void,
ElemSpaceDiscr_ghost_rhoZ_left_wall,
ElemSpaceDiscr_ghost_rhoZ_left_inflow,
ElemSpaceDiscr_ghost_rhoZ_left_outflow,
ElemSpaceDiscr_ghost_rhoZ_left_periodic,
ElemSpaceDiscr_ghost_rhoZ_left_neumann,
ElemSpaceDiscr_ghost_rhoZ_left_neumann
};


/*------------------------------------------------------------------------------
 right rho
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rho_right_void(
												double* rho, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rho_right_wall(
												double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_right_inflow(
												  double* rho, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_right_outflow(
												   double* rho, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_right_periodic(
													double* rho, 
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
				rho[m + nobj] = rho[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_right_neumann(
												   double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rho_right[])(
												double* rho, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rho_right_void,
ElemSpaceDiscr_ghost_rho_right_wall,
ElemSpaceDiscr_ghost_rho_right_inflow,
ElemSpaceDiscr_ghost_rho_right_outflow,
ElemSpaceDiscr_ghost_rho_right_periodic,
ElemSpaceDiscr_ghost_rho_right_neumann,
ElemSpaceDiscr_ghost_rho_right_neumann
};


/*------------------------------------------------------------------------------
 right rhou
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhou_right_void(
												 double* rhou, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhou_right_wall(
												 double* rhou, 
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
				rhou[n + mobj] = -rhou[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_right_inflow(
												   double* rhou, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_right_outflow(
													double* rhou, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_right_periodic(
													 double* rhou, 
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
				rhou[m + nobj] = rhou[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_right_neumann(
													double* rhou, 
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
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhou_right[])(
												 double* rhou, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhou_right_void,
ElemSpaceDiscr_ghost_rhou_right_wall,
ElemSpaceDiscr_ghost_rhou_right_inflow,
ElemSpaceDiscr_ghost_rhou_right_outflow,
ElemSpaceDiscr_ghost_rhou_right_periodic,
ElemSpaceDiscr_ghost_rhou_right_neumann,
ElemSpaceDiscr_ghost_rhou_right_neumann
};


/*------------------------------------------------------------------------------
 right rhov
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhov_right_void(
												 double* rhov, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhov_right_wall(
												 double* rhov, 
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
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_right_inflow(
												   double* rhov, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_right_outflow(
													double* rhov, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_right_periodic(
													 double* rhov, 
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
				rhov[m + nobj] = rhov[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_right_neumann(
													double* rhov, 
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
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhov_right[])(
												 double* rhov, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhov_right_void,
ElemSpaceDiscr_ghost_rhov_right_wall,
ElemSpaceDiscr_ghost_rhov_right_inflow,
ElemSpaceDiscr_ghost_rhov_right_outflow,
ElemSpaceDiscr_ghost_rhov_right_periodic,
ElemSpaceDiscr_ghost_rhov_right_neumann,
ElemSpaceDiscr_ghost_rhov_right_neumann
};


/*------------------------------------------------------------------------------
 right rhow
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhow_right_void(
												 double* rhow, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhow_right_wall(
												 double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_right_inflow(
												   double* rhow, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_right_outflow(
													double* rhow, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_right_periodic(
													 double* rhow, 
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
				rhow[m + nobj] = rhow[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_right_neumann(
													double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhow_right[])(
												 double* rhow, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhow_right_void,
ElemSpaceDiscr_ghost_rhow_right_wall,
ElemSpaceDiscr_ghost_rhow_right_inflow,
ElemSpaceDiscr_ghost_rhow_right_outflow,
ElemSpaceDiscr_ghost_rhow_right_periodic,
ElemSpaceDiscr_ghost_rhow_right_neumann,
ElemSpaceDiscr_ghost_rhow_right_neumann
};


/*------------------------------------------------------------------------------
 right rhoe
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoe_right_void(
												 double* rhoe, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhoe_right_wall(
												 double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_right_inflow(
												   double* rhoe, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_right_outflow(
													double* rhoe, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_right_periodic(
													 double* rhoe, 
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
				rhoe[m + nobj] = rhoe[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_right_neumann(
													double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhoe_right[])(
												 double* rhoe, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhoe_right_void,
ElemSpaceDiscr_ghost_rhoe_right_wall,
ElemSpaceDiscr_ghost_rhoe_right_inflow,
ElemSpaceDiscr_ghost_rhoe_right_outflow,
ElemSpaceDiscr_ghost_rhoe_right_periodic,
ElemSpaceDiscr_ghost_rhoe_right_neumann,
ElemSpaceDiscr_ghost_rhoe_right_neumann
};


/*------------------------------------------------------------------------------
 right rhoY
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoY_right_void(
												 double* rhoY, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhoY_right_wall(
												 double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_right_inflow(
												   double* rhoY, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_right_outflow(
													double* rhoY, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_right_periodic(
													 double* rhoY, 
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
				rhoY[m + nobj] = rhoY[m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_right_neumann(
													double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhoY_right[])(
												 double* rhoY, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhoY_right_void,
ElemSpaceDiscr_ghost_rhoY_right_wall,
ElemSpaceDiscr_ghost_rhoY_right_inflow,
ElemSpaceDiscr_ghost_rhoY_right_outflow,
ElemSpaceDiscr_ghost_rhoY_right_periodic,
ElemSpaceDiscr_ghost_rhoY_right_neumann,
ElemSpaceDiscr_ghost_rhoY_right_neumann
};


/*------------------------------------------------------------------------------
 right rhoZ
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoZ_right_void(
												 double** rhoZ, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhoZ_right_wall(
												 double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_right_inflow(
												   double** rhoZ, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {  
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_right_outflow(
													double** rhoZ, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_right_periodic(
													 double** rhoZ, 
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
				rhoZ[PRES][m + nobj] = rhoZ[PRES][m + nsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_right_neumann(
													double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void (*ElemSpaceDiscr_ghost_rhoZ_right[])(
												 double** rhoZ, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhoZ_right_void,
ElemSpaceDiscr_ghost_rhoZ_right_wall,
ElemSpaceDiscr_ghost_rhoZ_right_inflow,
ElemSpaceDiscr_ghost_rhoZ_right_outflow,
ElemSpaceDiscr_ghost_rhoZ_right_periodic,
ElemSpaceDiscr_ghost_rhoZ_right_neumann,
ElemSpaceDiscr_ghost_rhoZ_right_neumann
};


/*------------------------------------------------------------------------------
 bottom rho
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rho_bottom_void(
												 double* rho, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rho_bottom_wall(
												 double* rho, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_bottom_inflow(
												   double* rho, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_bottom_outflow(
													double* rho, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_bottom_periodic(
													 double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_bottom_neumann(
													double* rho, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rho_bottom[])(
												 double* rho, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rho_bottom_void,
ElemSpaceDiscr_ghost_rho_bottom_wall,
ElemSpaceDiscr_ghost_rho_bottom_inflow,
ElemSpaceDiscr_ghost_rho_bottom_outflow,
ElemSpaceDiscr_ghost_rho_bottom_periodic,
ElemSpaceDiscr_ghost_rho_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhou
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhou_bottom_void(
												  double* rhou, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {}

static void ElemSpaceDiscr_ghost_rhou_bottom_wall(
												  double* rhou, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_bottom_inflow(
													double* rhou, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_bottom_outflow(
													 double* rhou, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_bottom_periodic(
													  double* rhou, 
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
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_bottom_neumann(
													 double* rhou, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhou_bottom[])(
												  double* rhou, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) = {
ElemSpaceDiscr_ghost_rhou_bottom_void,
ElemSpaceDiscr_ghost_rhou_bottom_wall,
ElemSpaceDiscr_ghost_rhou_bottom_inflow,
ElemSpaceDiscr_ghost_rhou_bottom_outflow,
ElemSpaceDiscr_ghost_rhou_bottom_periodic,
ElemSpaceDiscr_ghost_rhou_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhov
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhov_bottom_void(
												  double* rhov, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {}

static void ElemSpaceDiscr_ghost_rhov_bottom_wall(
												  double* rhov, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhov[n + mobj] = - rhov[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_bottom_inflow(
													double* rhov, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_bottom_outflow(
													 double* rhov, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_bottom_periodic(
													  double* rhov, 
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
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_bottom_neumann(
													 double* rhov, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhov_bottom[])(
												  double* rhov, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) = {
ElemSpaceDiscr_ghost_rhov_bottom_void,
ElemSpaceDiscr_ghost_rhov_bottom_wall,
ElemSpaceDiscr_ghost_rhov_bottom_inflow,
ElemSpaceDiscr_ghost_rhov_bottom_outflow,
ElemSpaceDiscr_ghost_rhov_bottom_periodic,
ElemSpaceDiscr_ghost_rhov_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhow
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhow_bottom_void(
												  double* rhow, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {}

static void ElemSpaceDiscr_ghost_rhow_bottom_wall(
												  double* rhow, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_bottom_inflow(
													double* rhow, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_bottom_outflow(
													 double* rhow, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_bottom_periodic(
													  double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_bottom_neumann(
													 double* rhow, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhow_bottom[])(
												  double* rhow, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) = {
ElemSpaceDiscr_ghost_rhow_bottom_void,
ElemSpaceDiscr_ghost_rhow_bottom_wall,
ElemSpaceDiscr_ghost_rhow_bottom_inflow,
ElemSpaceDiscr_ghost_rhow_bottom_outflow,
ElemSpaceDiscr_ghost_rhow_bottom_periodic,
ElemSpaceDiscr_ghost_rhow_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhoe
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoe_bottom_void(
												  double* rhoe, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {}

static void ElemSpaceDiscr_ghost_rhoe_bottom_wall(
												  double* rhoe, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_bottom_inflow(
													double* rhoe, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_bottom_outflow(
													 double* rhoe, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_bottom_periodic(
													  double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_bottom_neumann(
													 double* rhoe, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoe_bottom[])(
												  double* rhoe, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) = {
ElemSpaceDiscr_ghost_rhoe_bottom_void,
ElemSpaceDiscr_ghost_rhoe_bottom_wall,
ElemSpaceDiscr_ghost_rhoe_bottom_inflow,
ElemSpaceDiscr_ghost_rhoe_bottom_outflow,
ElemSpaceDiscr_ghost_rhoe_bottom_periodic,
ElemSpaceDiscr_ghost_rhoe_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhoY
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoY_bottom_void(
												  double* rhoY, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {}

static void ElemSpaceDiscr_ghost_rhoY_bottom_wall(
												  double* rhoY, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_bottom_inflow(
													double* rhoY, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_bottom_outflow(
													 double* rhoY, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_bottom_periodic(
													  double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_bottom_neumann(
													 double* rhoY, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoY_bottom[])(
												  double* rhoY, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) = {
ElemSpaceDiscr_ghost_rhoY_bottom_void,
ElemSpaceDiscr_ghost_rhoY_bottom_wall,
ElemSpaceDiscr_ghost_rhoY_bottom_inflow,
ElemSpaceDiscr_ghost_rhoY_bottom_outflow,
ElemSpaceDiscr_ghost_rhoY_bottom_periodic,
ElemSpaceDiscr_ghost_rhoY_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhoZ
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoZ_bottom_void(
												  double** rhoZ, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {}

static void ElemSpaceDiscr_ghost_rhoZ_bottom_wall(
												  double** rhoZ, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = igy;
	const int bobj = igy - 1;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj - j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_bottom_inflow(
													double** rhoZ, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_bottom_outflow(
													 double** rhoZ, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_bottom_periodic(
													  double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_bottom_neumann(
													 double** rhoZ, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoZ_bottom[])(
												  double** rhoZ, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) = {
ElemSpaceDiscr_ghost_rhoZ_bottom_void,
ElemSpaceDiscr_ghost_rhoZ_bottom_wall,
ElemSpaceDiscr_ghost_rhoZ_bottom_inflow,
ElemSpaceDiscr_ghost_rhoZ_bottom_outflow,
ElemSpaceDiscr_ghost_rhoZ_bottom_periodic,
ElemSpaceDiscr_ghost_rhoZ_bottom_neumann
};

/*------------------------------------------------------------------------------
 top rho
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rho_top_void(
											  double* rho, 
											  const ElemSpaceDiscr* elem, 
											  const int ig) {}

static void ElemSpaceDiscr_ghost_rho_top_wall(
											  double* rho, 
											  const ElemSpaceDiscr* elem, 
											  const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_top_inflow(
												double* rho, 
												const ElemSpaceDiscr* elem, 
												const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_top_outflow(
												 double* rho, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_top_periodic(
												  double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_top_neumann(
												 double* rho, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rho_top[])(
											  double* rho, 
											  const ElemSpaceDiscr* elem, 
											  const int ig) = {
ElemSpaceDiscr_ghost_rho_top_void,
ElemSpaceDiscr_ghost_rho_top_wall,
ElemSpaceDiscr_ghost_rho_top_inflow,
ElemSpaceDiscr_ghost_rho_top_outflow,
ElemSpaceDiscr_ghost_rho_top_periodic,
ElemSpaceDiscr_ghost_rho_top_neumann
};


/*------------------------------------------------------------------------------
 top rhou
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhou_top_void(
											   double* rhou, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rhou_top_wall(
											   double* rhou, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_top_inflow(
												 double* rhou, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_top_outflow(
												  double* rhou, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_top_periodic(
												   double* rhou, 
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
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}


static void ElemSpaceDiscr_ghost_rhou_top_neumann(
												  double* rhou, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhou_top[])(
											   double* rhou, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rhou_top_void,
ElemSpaceDiscr_ghost_rhou_top_wall,
ElemSpaceDiscr_ghost_rhou_top_inflow,
ElemSpaceDiscr_ghost_rhou_top_outflow,
ElemSpaceDiscr_ghost_rhou_top_periodic,
ElemSpaceDiscr_ghost_rhou_top_neumann
};


/*------------------------------------------------------------------------------
 top rhov
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhov_top_void(
											   double* rhov, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rhov_top_wall(
											   double* rhov, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhov[n + mobj] = - rhov[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_top_inflow(
												 double* rhov, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_top_outflow(
												  double* rhov, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_top_periodic(
												   double* rhov, 
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
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_top_neumann(
												  double* rhov, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhov_top[])(
											   double* rhov, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rhov_top_void,
ElemSpaceDiscr_ghost_rhov_top_wall,
ElemSpaceDiscr_ghost_rhov_top_inflow,
ElemSpaceDiscr_ghost_rhov_top_outflow,
ElemSpaceDiscr_ghost_rhov_top_periodic,
ElemSpaceDiscr_ghost_rhov_top_neumann
};


/*------------------------------------------------------------------------------
 top rhow
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhow_top_void(
											   double* rhow, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rhow_top_wall(
											   double* rhow, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_top_inflow(
												 double* rhow, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_top_outflow(
												  double* rhow, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_top_periodic(
												   double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_top_neumann(
												  double* rhow, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhow_top[])(
											   double* rhow, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rhow_top_void,
ElemSpaceDiscr_ghost_rhow_top_wall,
ElemSpaceDiscr_ghost_rhow_top_inflow,
ElemSpaceDiscr_ghost_rhow_top_outflow,
ElemSpaceDiscr_ghost_rhow_top_periodic,
ElemSpaceDiscr_ghost_rhow_top_neumann
};


/*------------------------------------------------------------------------------
 top rhoe
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoe_top_void(
											   double* rhoe, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rhoe_top_wall(
											   double* rhoe, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_top_inflow(
												 double* rhoe, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_top_outflow(
												  double* rhoe, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_top_periodic(
												   double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_top_neumann(
												  double* rhoe, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoe_top[])(
											   double* rhoe, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rhoe_top_void,
ElemSpaceDiscr_ghost_rhoe_top_wall,
ElemSpaceDiscr_ghost_rhoe_top_inflow,
ElemSpaceDiscr_ghost_rhoe_top_outflow,
ElemSpaceDiscr_ghost_rhoe_top_periodic,
ElemSpaceDiscr_ghost_rhoe_top_neumann
};


/*------------------------------------------------------------------------------
 top rhoY
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoY_top_void(
											   double* rhoY, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rhoY_top_wall(
											   double* rhoY, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_top_inflow(
												 double* rhoY, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_top_outflow(
												  double* rhoY, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_top_periodic(
												   double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_top_neumann(
												  double* rhoY, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoY_top[])(
											   double* rhoY, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rhoY_top_void,
ElemSpaceDiscr_ghost_rhoY_top_wall,
ElemSpaceDiscr_ghost_rhoY_top_inflow,
ElemSpaceDiscr_ghost_rhoY_top_outflow,
ElemSpaceDiscr_ghost_rhoY_top_periodic,
ElemSpaceDiscr_ghost_rhoY_top_neumann
};


/*------------------------------------------------------------------------------
 top rhoZ
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoZ_top_void(
											   double** rhoZ, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rhoZ_top_wall(
											   double** rhoZ, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {
	
	const int igy = elem->igy;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int bsrc = icy - igy - 1;
	const int bobj = icy - igy;
	int i, j, k, l, msrc, mobj, n;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc - j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_top_inflow(
												 double** rhoZ, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_top_outflow(
												  double** rhoZ, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_top_periodic(
												   double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_top_neumann(
												  double** rhoZ, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoZ_top[])(
											   double** rhoZ, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rhoZ_top_void,
ElemSpaceDiscr_ghost_rhoZ_top_wall,
ElemSpaceDiscr_ghost_rhoZ_top_inflow,
ElemSpaceDiscr_ghost_rhoZ_top_outflow,
ElemSpaceDiscr_ghost_rhoZ_top_periodic,
ElemSpaceDiscr_ghost_rhoZ_top_neumann
};

/*------------------------------------------------------------------------------
 back rho
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rho_back_void(
											   double* rho, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) {}

static void ElemSpaceDiscr_ghost_rho_back_wall(
											   double* rho, 
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
				rho[n + lobj] = rho[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_back_inflow(
												 double* rho, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_back_outflow(
												  double* rho, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_back_periodic(
												   double* rho, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_back_neumann(
												  double* rho, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rho_back[])(
											   double* rho, 
											   const ElemSpaceDiscr* elem, 
											   const int ig) = {
ElemSpaceDiscr_ghost_rho_back_void,
ElemSpaceDiscr_ghost_rho_back_wall,
ElemSpaceDiscr_ghost_rho_back_inflow,
ElemSpaceDiscr_ghost_rho_back_outflow,
ElemSpaceDiscr_ghost_rho_back_periodic,
ElemSpaceDiscr_ghost_rho_back_neumann
};


/*------------------------------------------------------------------------------
 back rhou
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhou_back_void(
												double* rhou, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhou_back_wall(
												double* rhou, 
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
				rhou[n + lobj] = rhou[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_back_inflow(
												  double* rhou, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_back_outflow(
												   double* rhou, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_back_periodic(
													double* rhou, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_back_neumann(
												   double* rhou, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhou_back[])(
												double* rhou, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhou_back_void,
ElemSpaceDiscr_ghost_rhou_back_wall,
ElemSpaceDiscr_ghost_rhou_back_inflow,
ElemSpaceDiscr_ghost_rhou_back_outflow,
ElemSpaceDiscr_ghost_rhou_back_periodic,
ElemSpaceDiscr_ghost_rhou_back_neumann
};


/*------------------------------------------------------------------------------
 back rhov
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhov_back_void(
												double* rhov, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhov_back_wall(
												double* rhov, 
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
				rhov[n + lobj] = rhov[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_back_inflow(
												  double* rhov, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_back_outflow(
												   double* rhov, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_back_periodic(
													double* rhov, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_back_neumann(
												   double* rhov, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhov_back[])(
												double* rhov, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhov_back_void,
ElemSpaceDiscr_ghost_rhov_back_wall,
ElemSpaceDiscr_ghost_rhov_back_inflow,
ElemSpaceDiscr_ghost_rhov_back_outflow,
ElemSpaceDiscr_ghost_rhov_back_periodic,
ElemSpaceDiscr_ghost_rhov_back_neumann
};

/*------------------------------------------------------------------------------
 back rhow
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhow_back_void(
												double* rhow, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhow_back_wall(
												double* rhow, 
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
				rhow[n + lobj] = - rhow[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_back_inflow(
												  double* rhow, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_back_outflow(
												   double* rhow, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_back_periodic(
													double* rhow, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_back_neumann(
												   double* rhow, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhow_back[])(
												double* rhow, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhow_back_void,
ElemSpaceDiscr_ghost_rhow_back_wall,
ElemSpaceDiscr_ghost_rhow_back_inflow,
ElemSpaceDiscr_ghost_rhow_back_outflow,
ElemSpaceDiscr_ghost_rhow_back_periodic,
ElemSpaceDiscr_ghost_rhow_back_neumann
};


/*------------------------------------------------------------------------------
 back rhoe
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoe_back_void(
												double* rhoe, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhoe_back_wall(
												double* rhoe, 
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
				rhoe[n + lobj] = rhoe[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_back_inflow(
												  double* rhoe, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_back_outflow(
												   double* rhoe, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_back_periodic(
													double* rhoe, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_back_neumann(
												   double* rhoe, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoe_back[])(
												double* rhoe, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhoe_back_void,
ElemSpaceDiscr_ghost_rhoe_back_wall,
ElemSpaceDiscr_ghost_rhoe_back_inflow,
ElemSpaceDiscr_ghost_rhoe_back_outflow,
ElemSpaceDiscr_ghost_rhoe_back_periodic,
ElemSpaceDiscr_ghost_rhoe_back_neumann
};


/*------------------------------------------------------------------------------
 back rhoY
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoY_back_void(
												double* rhoY, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhoY_back_wall(
												double* rhoY, 
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
				rhoY[n + lobj] = rhoY[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_back_inflow(
												  double* rhoY, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_back_outflow(
												   double* rhoY, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_back_periodic(
													double* rhoY, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_back_neumann(
												   double* rhoY, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoY_back[])(
												double* rhoY, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhoY_back_void,
ElemSpaceDiscr_ghost_rhoY_back_wall,
ElemSpaceDiscr_ghost_rhoY_back_inflow,
ElemSpaceDiscr_ghost_rhoY_back_outflow,
ElemSpaceDiscr_ghost_rhoY_back_periodic,
ElemSpaceDiscr_ghost_rhoY_back_neumann
};


/*------------------------------------------------------------------------------
 back rhoZ
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoZ_back_void(
												double** rhoZ, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rhoZ_back_wall(
												double** rhoZ, 
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
				rhoZ[PRES][n + lobj] = rhoZ[PRES][n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_back_inflow(
												  double** rhoZ, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_back_outflow(
												   double** rhoZ, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_back_periodic(
													double** rhoZ, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_back_neumann(
												   double** rhoZ, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoZ_back[])(
												double** rhoZ, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rhoZ_back_void,
ElemSpaceDiscr_ghost_rhoZ_back_wall,
ElemSpaceDiscr_ghost_rhoZ_back_inflow,
ElemSpaceDiscr_ghost_rhoZ_back_outflow,
ElemSpaceDiscr_ghost_rhoZ_back_periodic,
ElemSpaceDiscr_ghost_rhoZ_back_neumann
};

/*------------------------------------------------------------------------------
 front rho
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rho_front_void(
												double* rho, 
												const ElemSpaceDiscr* elem, 
												const int ig) {}

static void ElemSpaceDiscr_ghost_rho_front_wall(
												double* rho, 
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
				rho[n + lobj] = rho[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rho_front_inflow(
												  double* rho, 
												  const ElemSpaceDiscr* elem, 
												  const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_front_outflow(
												   double* rho, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_front_periodic(
													double* rho, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rho_front_neumann(
												   double* rho, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rho_front[])(
												double* rho, 
												const ElemSpaceDiscr* elem, 
												const int ig) = {
ElemSpaceDiscr_ghost_rho_front_void,
ElemSpaceDiscr_ghost_rho_front_wall,
ElemSpaceDiscr_ghost_rho_front_inflow,
ElemSpaceDiscr_ghost_rho_front_outflow,
ElemSpaceDiscr_ghost_rho_front_periodic,
ElemSpaceDiscr_ghost_rho_front_neumann
};


/*------------------------------------------------------------------------------
 front rhou
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhou_front_void(
												 double* rhou, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhou_front_wall(
												 double* rhou, 
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
				rhou[n + lobj] = rhou[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhou_front_inflow(
												   double* rhou, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_front_outflow(
													double* rhou, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_front_periodic(
													 double* rhou, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhou_front_neumann(
													double* rhou, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhou_front[])(
												 double* rhou, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhou_front_void,
ElemSpaceDiscr_ghost_rhou_front_wall,
ElemSpaceDiscr_ghost_rhou_front_inflow,
ElemSpaceDiscr_ghost_rhou_front_outflow,
ElemSpaceDiscr_ghost_rhou_front_periodic,
ElemSpaceDiscr_ghost_rhou_front_neumann
};


/*------------------------------------------------------------------------------
 front rhov
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhov_front_void(
												 double* rhov, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhov_front_wall(
												 double* rhov, 
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
				rhov[n + lobj] = rhov[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhov_front_inflow(
												   double* rhov, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_front_outflow(
													double* rhov, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_front_periodic(
													 double* rhov, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhov_front_neumann(
													double* rhov, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhov_front[])(
												 double* rhov, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhov_front_void,
ElemSpaceDiscr_ghost_rhov_front_wall,
ElemSpaceDiscr_ghost_rhov_front_inflow,
ElemSpaceDiscr_ghost_rhov_front_outflow,
ElemSpaceDiscr_ghost_rhov_front_periodic,
ElemSpaceDiscr_ghost_rhov_front_neumann
};


/*------------------------------------------------------------------------------
 front rhow
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhow_front_void(
												 double* rhow, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhow_front_wall(
												 double* rhow, 
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
				rhow[n + lobj] = - rhow[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhow_front_inflow(
												   double* rhow, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_front_outflow(
													double* rhow, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_front_periodic(
													 double* rhow, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhow_front_neumann(
													double* rhow, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhow_front[])(
												 double* rhow, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhow_front_void,
ElemSpaceDiscr_ghost_rhow_front_wall,
ElemSpaceDiscr_ghost_rhow_front_inflow,
ElemSpaceDiscr_ghost_rhow_front_outflow,
ElemSpaceDiscr_ghost_rhow_front_periodic,
ElemSpaceDiscr_ghost_rhow_front_neumann
};


/*------------------------------------------------------------------------------
 front rhoe
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoe_front_void(
												 double* rhoe, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhoe_front_wall(
												 double* rhoe, 
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
				rhoe[n + lobj] = rhoe[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoe_front_inflow(
												   double* rhoe, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_front_outflow(
													double* rhoe, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_front_periodic(
													 double* rhoe, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoe_front_neumann(
													double* rhoe, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoe_front[])(
												 double* rhoe, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhoe_front_void,
ElemSpaceDiscr_ghost_rhoe_front_wall,
ElemSpaceDiscr_ghost_rhoe_front_inflow,
ElemSpaceDiscr_ghost_rhoe_front_outflow,
ElemSpaceDiscr_ghost_rhoe_front_periodic,
ElemSpaceDiscr_ghost_rhoe_front_neumann
};


/*------------------------------------------------------------------------------
 front rhoY
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoY_front_void(
												 double* rhoY, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhoY_front_wall(
												 double* rhoY, 
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
				rhoY[n + lobj] = rhoY[n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoY_front_inflow(
												   double* rhoY, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_front_outflow(
													double* rhoY, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_front_periodic(
													 double* rhoY, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoY_front_neumann(
													double* rhoY, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoY_front[])(
												 double* rhoY, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhoY_front_void,
ElemSpaceDiscr_ghost_rhoY_front_wall,
ElemSpaceDiscr_ghost_rhoY_front_inflow,
ElemSpaceDiscr_ghost_rhoY_front_outflow,
ElemSpaceDiscr_ghost_rhoY_front_periodic,
ElemSpaceDiscr_ghost_rhoY_front_neumann
};


/*------------------------------------------------------------------------------
 front rhoZ
 ------------------------------------------------------------------------------*/
static void ElemSpaceDiscr_ghost_rhoZ_front_void(
												 double** rhoZ, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) {}

static void ElemSpaceDiscr_ghost_rhoZ_front_wall(
												 double** rhoZ, 
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
				rhoZ[PRES][n + lobj] = rhoZ[PRES][n + lsrc];
			}
		}
	}
}

static void ElemSpaceDiscr_ghost_rhoZ_front_inflow(
												   double** rhoZ, 
												   const ElemSpaceDiscr* elem, 
												   const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_front_outflow(
													double** rhoZ, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_front_periodic(
													 double** rhoZ, 
													 const ElemSpaceDiscr* elem, 
													 const int ig) {
	ERROR("function not available");
}

static void ElemSpaceDiscr_ghost_rhoZ_front_neumann(
													double** rhoZ, 
													const ElemSpaceDiscr* elem, 
													const int ig) {
	ERROR("function not available");
}

static void (*ElemSpaceDiscr_ghost_rhoZ_front[])(
												 double** rhoZ, 
												 const ElemSpaceDiscr* elem, 
												 const int ig) = {
ElemSpaceDiscr_ghost_rhoZ_front_void,
ElemSpaceDiscr_ghost_rhoZ_front_wall,
ElemSpaceDiscr_ghost_rhoZ_front_inflow,
ElemSpaceDiscr_ghost_rhoZ_front_outflow,
ElemSpaceDiscr_ghost_rhoZ_front_periodic,
ElemSpaceDiscr_ghost_rhoZ_front_neumann
};


/*------------------------------------------------------------------------------
 conservative variables
 ------------------------------------------------------------------------------*/

void ElemSpaceDiscr_ghost_rho(
							  double* rho, 
							  const ElemSpaceDiscr* elem, 
							  const int ig) {
	
	(*ElemSpaceDiscr_ghost_rho_left[elem->left])(rho, elem, ig);
	(*ElemSpaceDiscr_ghost_rho_right[elem->right])(rho, elem, ig);
	(*ElemSpaceDiscr_ghost_rho_bottom[elem->bottom])(rho, elem, ig);
	(*ElemSpaceDiscr_ghost_rho_top[elem->top])(rho, elem, ig);
	(*ElemSpaceDiscr_ghost_rho_back[elem->back])(rho, elem, ig);
	(*ElemSpaceDiscr_ghost_rho_front[elem->front])(rho, elem, ig);
}


void ElemSpaceDiscr_ghost_rhou(
							   double* rhou, 
							   const ElemSpaceDiscr* elem, 
							   const int ig) {
	
	(*ElemSpaceDiscr_ghost_rhou_left[elem->left])(rhou, elem, ig);
	(*ElemSpaceDiscr_ghost_rhou_right[elem->right])(rhou, elem, ig);
	(*ElemSpaceDiscr_ghost_rhou_bottom[elem->bottom])(rhou, elem, ig);
	(*ElemSpaceDiscr_ghost_rhou_top[elem->top])(rhou, elem, ig);
	(*ElemSpaceDiscr_ghost_rhou_back[elem->back])(rhou, elem, ig);
	(*ElemSpaceDiscr_ghost_rhou_front[elem->front])(rhou, elem, ig);
}


void ElemSpaceDiscr_ghost_rhov(
							   double* rhov, 
							   const ElemSpaceDiscr* elem, 
							   const int ig) {
	
	(*ElemSpaceDiscr_ghost_rhov_left[elem->left])(rhov, elem, ig);
	(*ElemSpaceDiscr_ghost_rhov_right[elem->right])(rhov, elem, ig);
	(*ElemSpaceDiscr_ghost_rhov_bottom[elem->bottom])(rhov, elem, ig);
	(*ElemSpaceDiscr_ghost_rhov_top[elem->top])(rhov, elem, ig);
	(*ElemSpaceDiscr_ghost_rhov_back[elem->back])(rhov, elem, ig);
	(*ElemSpaceDiscr_ghost_rhov_front[elem->front])(rhov, elem, ig);
}


void ElemSpaceDiscr_ghost_rhow(
							   double* rhow, 
							   const ElemSpaceDiscr* elem, 
							   const int ig) {
	
	(*ElemSpaceDiscr_ghost_rhow_left[elem->left])(rhow, elem, ig);
	(*ElemSpaceDiscr_ghost_rhow_right[elem->right])(rhow, elem, ig);
	(*ElemSpaceDiscr_ghost_rhow_bottom[elem->bottom])(rhow, elem, ig);
	(*ElemSpaceDiscr_ghost_rhow_top[elem->top])(rhow, elem, ig);
	(*ElemSpaceDiscr_ghost_rhow_back[elem->back])(rhow, elem, ig);
	(*ElemSpaceDiscr_ghost_rhow_front[elem->front])(rhow, elem, ig);
}


void ElemSpaceDiscr_ghost_rhoe(
							   double* rhoe, 
							   const ElemSpaceDiscr* elem, 
							   const int ig) {
	
	(*ElemSpaceDiscr_ghost_rhoe_left[elem->left])(rhoe, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoe_right[elem->right])(rhoe, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoe_bottom[elem->bottom])(rhoe, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoe_top[elem->top])(rhoe, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoe_back[elem->back])(rhoe, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoe_front[elem->front])(rhoe, elem, ig);
}


void ElemSpaceDiscr_ghost_rhoY(
							   double* rhoY, 
							   const ElemSpaceDiscr* elem, 
							   const int ig) {
	
	(*ElemSpaceDiscr_ghost_rhoY_left[elem->left])(rhoY, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoY_right[elem->right])(rhoY, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoY_bottom[elem->bottom])(rhoY, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoY_top[elem->top])(rhoY, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoY_back[elem->back])(rhoY, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoY_front[elem->front])(rhoY, elem, ig);
}


void ElemSpaceDiscr_ghost_rhoZ(
							   double** rhoZ, 
							   const ElemSpaceDiscr* elem, 
							   const int ig) {
	
	(*ElemSpaceDiscr_ghost_rhoZ_left[elem->left])(rhoZ, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoZ_right[elem->right])(rhoZ, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoZ_bottom[elem->bottom])(rhoZ, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoZ_top[elem->top])(rhoZ, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoZ_back[elem->back])(rhoZ, elem, ig);
	(*ElemSpaceDiscr_ghost_rhoZ_front[elem->front])(rhoZ, elem, ig);
}


void ElemSpaceDiscr_ghost(
						  ConsVars* Sol, 
						  const ElemSpaceDiscr* elem, 
						  const int ig) {
	
	ElemSpaceDiscr_ghost_rho(Sol->rho, elem, ig);
	ElemSpaceDiscr_ghost_rhou(Sol->rhou, elem, ig);
	ElemSpaceDiscr_ghost_rhov(Sol->rhov, elem, ig);
	ElemSpaceDiscr_ghost_rhow(Sol->rhow, elem, ig);
	ElemSpaceDiscr_ghost_rhoe(Sol->rhoe, elem, ig);
	ElemSpaceDiscr_ghost_rhoY(Sol->rhoY, elem, ig);
	ElemSpaceDiscr_ghost_rhoZ(Sol->rhoZ, elem, ig);
}


/*------------------------------------------------------------------------------
 left rho
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rho_left_void(
											   double* rho, 
											   const NodeSpaceDiscr* node,
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rho_left_wall(
											   double* rho, 
											   const NodeSpaceDiscr* node,
											   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_left_inflow(
												 double* rho, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_left_outflow(
												  double* rho, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_left_periodic(
												   double* rho, 
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
				rho[m + nobj] = rho[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_left_neumann(
												  double* rho, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rho_left[])(
											   double* rho, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rho_left_void,
NodeSpaceDiscr_ghost_rho_left_wall,
NodeSpaceDiscr_ghost_rho_left_inflow,
NodeSpaceDiscr_ghost_rho_left_outflow,
NodeSpaceDiscr_ghost_rho_left_periodic,
NodeSpaceDiscr_ghost_rho_left_neumann,
NodeSpaceDiscr_ghost_rho_left_neumann
};


/*------------------------------------------------------------------------------
 left rhou
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhou_left_void(
												double* rhou, 
												const NodeSpaceDiscr* node,
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhou_left_wall(
												double* rhou, 
												const NodeSpaceDiscr* node,
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_left_inflow(
												  double* rhou, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_left_outflow(
												   double* rhou, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_left_periodic(
													double* rhou, 
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
				rhou[m + nobj] = rhou[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_left_neumann(
												   double* rhou, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhou_left[])(
												double* rhou, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhou_left_void,
NodeSpaceDiscr_ghost_rhou_left_wall,
NodeSpaceDiscr_ghost_rhou_left_inflow,
NodeSpaceDiscr_ghost_rhou_left_outflow,
NodeSpaceDiscr_ghost_rhou_left_periodic,
NodeSpaceDiscr_ghost_rhou_left_neumann,
NodeSpaceDiscr_ghost_rhou_left_neumann
};


/*------------------------------------------------------------------------------
 left rhov
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhov_left_void(
												double* rhov, 
												const NodeSpaceDiscr* node,
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhov_left_wall(
												double* rhov, 
												const NodeSpaceDiscr* node,
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_left_inflow(
												  double* rhov, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_left_outflow(
												   double* rhov, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_left_periodic(
													double* rhov, 
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
				rhov[m + nobj] = rhov[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_left_neumann(
												   double* rhov, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhov_left[])(
												double* rhov, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhov_left_void,
NodeSpaceDiscr_ghost_rhov_left_wall,
NodeSpaceDiscr_ghost_rhov_left_inflow,
NodeSpaceDiscr_ghost_rhov_left_outflow,
NodeSpaceDiscr_ghost_rhov_left_periodic,
NodeSpaceDiscr_ghost_rhov_left_neumann,
NodeSpaceDiscr_ghost_rhov_left_neumann
};


/*------------------------------------------------------------------------------
 left rhow
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhow_left_void(
												double* rhow, 
												const NodeSpaceDiscr* node,
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhow_left_wall(
												double* rhow, 
												const NodeSpaceDiscr* node,
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_left_inflow(
												  double* rhow, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_left_outflow(
												   double* rhow, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_left_periodic(
													double* rhow, 
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
				rhow[m + nobj] = rhow[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_left_neumann(
												   double* rhow, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhow_left[])(
												double* rhow, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhow_left_void,
NodeSpaceDiscr_ghost_rhow_left_wall,
NodeSpaceDiscr_ghost_rhow_left_inflow,
NodeSpaceDiscr_ghost_rhow_left_outflow,
NodeSpaceDiscr_ghost_rhow_left_periodic,
NodeSpaceDiscr_ghost_rhow_left_neumann,
NodeSpaceDiscr_ghost_rhow_left_neumann
};


/*------------------------------------------------------------------------------
 left rhoe
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoe_left_void(
												double* rhoe, 
												const NodeSpaceDiscr* node,
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhoe_left_wall(
												double* rhoe, 
												const NodeSpaceDiscr* node,
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_left_inflow(
												  double* rhoe, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_left_outflow(
												   double* rhoe, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_left_periodic(
													double* rhoe, 
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
				rhoe[m + nobj] = rhoe[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_left_neumann(
												   double* rhoe, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoe_left[])(
												double* rhoe, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhoe_left_void,
NodeSpaceDiscr_ghost_rhoe_left_wall,
NodeSpaceDiscr_ghost_rhoe_left_inflow,
NodeSpaceDiscr_ghost_rhoe_left_outflow,
NodeSpaceDiscr_ghost_rhoe_left_periodic,
NodeSpaceDiscr_ghost_rhoe_left_neumann,
NodeSpaceDiscr_ghost_rhoe_left_neumann
};


/*------------------------------------------------------------------------------
 left rhoY
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoY_left_void(
												double* rhoY, 
												const NodeSpaceDiscr* node,
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhoY_left_wall(
												double* rhoY, 
												const NodeSpaceDiscr* node,
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_left_inflow(
												  double* rhoY, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_left_outflow(
												   double* rhoY, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_left_periodic(
													double* rhoY, 
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
				rhoY[m + nobj] = rhoY[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_left_neumann(
												   double* rhoY, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoY_left[])(
												double* rhoY, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhoY_left_void,
NodeSpaceDiscr_ghost_rhoY_left_wall,
NodeSpaceDiscr_ghost_rhoY_left_inflow,
NodeSpaceDiscr_ghost_rhoY_left_outflow,
NodeSpaceDiscr_ghost_rhoY_left_periodic,
NodeSpaceDiscr_ghost_rhoY_left_neumann,
NodeSpaceDiscr_ghost_rhoY_left_neumann
};


/*------------------------------------------------------------------------------
 left rhoZ
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoZ_left_void(
												double** rhoZ, 
												const NodeSpaceDiscr* node,
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhoZ_left_wall(
												double** rhoZ, 
												const NodeSpaceDiscr* node,
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_left_inflow(
												  double** rhoZ, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_left_outflow(
												   double** rhoZ, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_left_periodic(
													double** rhoZ, 
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
				rhoZ[PRES][m + nobj] = rhoZ[PRES][m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_left_neumann(
												   double** rhoZ, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoZ_left[])(
												double** rhoZ, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhoZ_left_void,
NodeSpaceDiscr_ghost_rhoZ_left_wall,
NodeSpaceDiscr_ghost_rhoZ_left_inflow,
NodeSpaceDiscr_ghost_rhoZ_left_outflow,
NodeSpaceDiscr_ghost_rhoZ_left_periodic,
NodeSpaceDiscr_ghost_rhoZ_left_neumann,
NodeSpaceDiscr_ghost_rhoZ_left_neumann
};



/*------------------------------------------------------------------------------
 right rho
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rho_right_void(
												double* rho, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rho_right_wall(
												double* rho, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_right_inflow(
												  double* rho, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_right_outflow(
												   double* rho, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_right_periodic(
													double* rho, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rho[m + nobj] = rho[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_right_neumann(
												   double* rho, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rho_right[])(
												double* rho, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rho_right_void,
NodeSpaceDiscr_ghost_rho_right_wall,
NodeSpaceDiscr_ghost_rho_right_inflow,
NodeSpaceDiscr_ghost_rho_right_outflow,
NodeSpaceDiscr_ghost_rho_right_periodic,
NodeSpaceDiscr_ghost_rho_right_neumann,
NodeSpaceDiscr_ghost_rho_right_neumann
};


/*------------------------------------------------------------------------------
 right rhou
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhou_right_void(
												 double* rhou, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhou_right_wall(
												 double* rhou, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_right_inflow(
												   double* rhou, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_right_outflow(
													double* rhou, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_right_periodic(
													 double* rhou, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhou[m + nobj] = rhou[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_right_neumann(
													double* rhou, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhou_right[])(
												 double* rhou, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhou_right_void,
NodeSpaceDiscr_ghost_rhou_right_wall,
NodeSpaceDiscr_ghost_rhou_right_inflow,
NodeSpaceDiscr_ghost_rhou_right_outflow,
NodeSpaceDiscr_ghost_rhou_right_periodic,
NodeSpaceDiscr_ghost_rhou_right_neumann,
NodeSpaceDiscr_ghost_rhou_right_neumann
};


/*------------------------------------------------------------------------------
 right rhov
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhov_right_void(
												 double* rhov, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhov_right_wall(
												 double* rhov, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_right_inflow(
												   double* rhov, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_right_outflow(
													double* rhov, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_right_periodic(
													 double* rhov, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhov[m + nobj] = rhov[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_right_neumann(
													double* rhov, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhov_right[])(
												 double* rhov, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhov_right_void,
NodeSpaceDiscr_ghost_rhov_right_wall,
NodeSpaceDiscr_ghost_rhov_right_inflow,
NodeSpaceDiscr_ghost_rhov_right_outflow,
NodeSpaceDiscr_ghost_rhov_right_periodic,
NodeSpaceDiscr_ghost_rhov_right_neumann,
NodeSpaceDiscr_ghost_rhov_right_neumann
};


/*------------------------------------------------------------------------------
 right rhow
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhow_right_void(
												 double* rhow, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhow_right_wall(
												 double* rhow, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_right_inflow(
												   double* rhow, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_right_outflow(
													double* rhow, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_right_periodic(
													 double* rhow, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhow[m + nobj] = rhow[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_right_neumann(
													double* rhow, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhow_right[])(
												 double* rhow, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhow_right_void,
NodeSpaceDiscr_ghost_rhow_right_wall,
NodeSpaceDiscr_ghost_rhow_right_inflow,
NodeSpaceDiscr_ghost_rhow_right_outflow,
NodeSpaceDiscr_ghost_rhow_right_periodic,
NodeSpaceDiscr_ghost_rhow_right_neumann,
NodeSpaceDiscr_ghost_rhow_right_neumann
};


/*------------------------------------------------------------------------------
 right rhoe
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoe_right_void(
												 double* rhoe, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhoe_right_wall(
												 double* rhoe, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_right_inflow(
												   double* rhoe, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_right_outflow(
													double* rhoe, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_right_periodic(
													 double* rhoe, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhoe[m + nobj] = rhoe[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_right_neumann(
													double* rhoe, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoe_right[])(
												 double* rhoe, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhoe_right_void,
NodeSpaceDiscr_ghost_rhoe_right_wall,
NodeSpaceDiscr_ghost_rhoe_right_inflow,
NodeSpaceDiscr_ghost_rhoe_right_outflow,
NodeSpaceDiscr_ghost_rhoe_right_periodic,
NodeSpaceDiscr_ghost_rhoe_right_neumann,
NodeSpaceDiscr_ghost_rhoe_right_neumann
};


/*------------------------------------------------------------------------------
 right rhoY
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoY_right_void(
												 double* rhoY, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhoY_right_wall(
												 double* rhoY, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_right_inflow(
												   double* rhoY, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_right_outflow(
													double* rhoY, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_right_periodic(
													 double* rhoY, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhoY[m + nobj] = rhoY[m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_right_neumann(
													double* rhoY, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoY_right[])(
												 double* rhoY, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhoY_right_void,
NodeSpaceDiscr_ghost_rhoY_right_wall,
NodeSpaceDiscr_ghost_rhoY_right_inflow,
NodeSpaceDiscr_ghost_rhoY_right_outflow,
NodeSpaceDiscr_ghost_rhoY_right_periodic,
NodeSpaceDiscr_ghost_rhoY_right_neumann,
NodeSpaceDiscr_ghost_rhoY_right_neumann
};


/*------------------------------------------------------------------------------
 right rhoZ
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoZ_right_void(
												 double** rhoZ, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhoZ_right_wall(
												 double** rhoZ, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_right_inflow(
												   double** rhoZ, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {  
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_right_outflow(
													double** rhoZ, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_right_periodic(
													 double** rhoZ, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igx + 1;
	const int bobj = icx - igx;
	int i, j, k, l, m, nsrc, nobj;
	
	assert(ig <= igx);
	
	for(i = 0; i < ig; i++) {
		nsrc = bsrc + i;
		nobj = bobj + i;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(j = 0; j < icy; j++) {m = l + j * icx;
				rhoZ[PRES][m + nobj] = rhoZ[PRES][m + nsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_right_neumann(
													double** rhoZ, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoZ_right[])(
												 double** rhoZ, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhoZ_right_void,
NodeSpaceDiscr_ghost_rhoZ_right_wall,
NodeSpaceDiscr_ghost_rhoZ_right_inflow,
NodeSpaceDiscr_ghost_rhoZ_right_outflow,
NodeSpaceDiscr_ghost_rhoZ_right_periodic,
NodeSpaceDiscr_ghost_rhoZ_right_neumann,
NodeSpaceDiscr_ghost_rhoZ_right_neumann
};



/*------------------------------------------------------------------------------
 bottom rho
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rho_bottom_void(
												 double* rho, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rho_bottom_wall(
												 double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_bottom_inflow(
												   double* rho, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_bottom_outflow(
													double* rho, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_bottom_periodic(
													 double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_bottom_neumann(
													double* rho, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rho_bottom[])(
												 double* rho, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rho_bottom_void,
NodeSpaceDiscr_ghost_rho_bottom_wall,
NodeSpaceDiscr_ghost_rho_bottom_inflow,
NodeSpaceDiscr_ghost_rho_bottom_outflow,
NodeSpaceDiscr_ghost_rho_bottom_periodic,
NodeSpaceDiscr_ghost_rho_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhou
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhou_bottom_void(
												  double* rhou, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {}

static void NodeSpaceDiscr_ghost_rhou_bottom_wall(
												  double* rhou, 
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
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_bottom_inflow(
													double* rhou, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_bottom_outflow(
													 double* rhou, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_bottom_periodic(
													  double* rhou, 
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
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_bottom_neumann(
													 double* rhou, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhou_bottom[])(
												  double* rhou, 
												  const NodeSpaceDiscr* node, 
												  const int ig) = {
NodeSpaceDiscr_ghost_rhou_bottom_void,
NodeSpaceDiscr_ghost_rhou_bottom_wall,
NodeSpaceDiscr_ghost_rhou_bottom_inflow,
NodeSpaceDiscr_ghost_rhou_bottom_outflow,
NodeSpaceDiscr_ghost_rhou_bottom_periodic,
NodeSpaceDiscr_ghost_rhou_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhov
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhov_bottom_void(
												  double* rhov, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {}

static void NodeSpaceDiscr_ghost_rhov_bottom_wall(
												  double* rhov, 
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
				rhov[n + mobj] = - rhov[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_bottom_inflow(
													double* rhov, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_bottom_outflow(
													 double* rhov, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_bottom_periodic(
													  double* rhov, 
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
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_bottom_neumann(
													 double* rhov, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhov_bottom[])(
												  double* rhov, 
												  const NodeSpaceDiscr* node, 
												  const int ig) = {
NodeSpaceDiscr_ghost_rhov_bottom_void,
NodeSpaceDiscr_ghost_rhov_bottom_wall,
NodeSpaceDiscr_ghost_rhov_bottom_inflow,
NodeSpaceDiscr_ghost_rhov_bottom_outflow,
NodeSpaceDiscr_ghost_rhov_bottom_periodic,
NodeSpaceDiscr_ghost_rhov_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhow
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhow_bottom_void(
												  double* rhow, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {}

static void NodeSpaceDiscr_ghost_rhow_bottom_wall(
												  double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_bottom_inflow(
													double* rhow, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_bottom_outflow(
													 double* rhow, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_bottom_periodic(
													  double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_bottom_neumann(
													 double* rhow, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhow_bottom[])(
												  double* rhow, 
												  const NodeSpaceDiscr* node, 
												  const int ig) = {
NodeSpaceDiscr_ghost_rhow_bottom_void,
NodeSpaceDiscr_ghost_rhow_bottom_wall,
NodeSpaceDiscr_ghost_rhow_bottom_inflow,
NodeSpaceDiscr_ghost_rhow_bottom_outflow,
NodeSpaceDiscr_ghost_rhow_bottom_periodic,
NodeSpaceDiscr_ghost_rhow_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhoe
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoe_bottom_void(
												  double* rhoe, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {}

static void NodeSpaceDiscr_ghost_rhoe_bottom_wall(
												  double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_bottom_inflow(
													double* rhoe, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_bottom_outflow(
													 double* rhoe, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_bottom_periodic(
													  double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_bottom_neumann(
													 double* rhoe, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoe_bottom[])(
												  double* rhoe, 
												  const NodeSpaceDiscr* node, 
												  const int ig) = {
NodeSpaceDiscr_ghost_rhoe_bottom_void,
NodeSpaceDiscr_ghost_rhoe_bottom_wall,
NodeSpaceDiscr_ghost_rhoe_bottom_inflow,
NodeSpaceDiscr_ghost_rhoe_bottom_outflow,
NodeSpaceDiscr_ghost_rhoe_bottom_periodic,
NodeSpaceDiscr_ghost_rhoe_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhoY
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoY_bottom_void(
												  double* rhoY, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {}

static void NodeSpaceDiscr_ghost_rhoY_bottom_wall(
												  double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_bottom_inflow(
													double* rhoY, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_bottom_outflow(
													 double* rhoY, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_bottom_periodic(
													  double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_bottom_neumann(
													 double* rhoY, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoY_bottom[])(
												  double* rhoY, 
												  const NodeSpaceDiscr* node, 
												  const int ig) = {
NodeSpaceDiscr_ghost_rhoY_bottom_void,
NodeSpaceDiscr_ghost_rhoY_bottom_wall,
NodeSpaceDiscr_ghost_rhoY_bottom_inflow,
NodeSpaceDiscr_ghost_rhoY_bottom_outflow,
NodeSpaceDiscr_ghost_rhoY_bottom_periodic,
NodeSpaceDiscr_ghost_rhoY_bottom_neumann
};


/*------------------------------------------------------------------------------
 bottom rhoZ
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoZ_bottom_void(
												  double** rhoZ, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {}

static void NodeSpaceDiscr_ghost_rhoZ_bottom_wall(
												  double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_bottom_inflow(
													double** rhoZ, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_bottom_outflow(
													 double** rhoZ, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_bottom_periodic(
													  double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_bottom_neumann(
													 double** rhoZ, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoZ_bottom[])(
												  double** rhoZ, 
												  const NodeSpaceDiscr* node, 
												  const int ig) = {
NodeSpaceDiscr_ghost_rhoZ_bottom_void,
NodeSpaceDiscr_ghost_rhoZ_bottom_wall,
NodeSpaceDiscr_ghost_rhoZ_bottom_inflow,
NodeSpaceDiscr_ghost_rhoZ_bottom_outflow,
NodeSpaceDiscr_ghost_rhoZ_bottom_periodic,
NodeSpaceDiscr_ghost_rhoZ_bottom_neumann
};


/*------------------------------------------------------------------------------
 top rho
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rho_top_void(
											  double* rho, 
											  const NodeSpaceDiscr* node, 
											  const int ig) {}

static void NodeSpaceDiscr_ghost_rho_top_wall(
											  double* rho, 
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
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_top_inflow(
												double* rho, 
												const NodeSpaceDiscr* node, 
												const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_top_outflow(
												 double* rho, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_top_periodic(
												  double* rho, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rho[n + mobj] = rho[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_top_neumann(
												 double* rho, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rho_top[])(
											  double* rho, 
											  const NodeSpaceDiscr* node, 
											  const int ig) = {
NodeSpaceDiscr_ghost_rho_top_void,
NodeSpaceDiscr_ghost_rho_top_wall,
NodeSpaceDiscr_ghost_rho_top_inflow,
NodeSpaceDiscr_ghost_rho_top_outflow,
NodeSpaceDiscr_ghost_rho_top_periodic,
NodeSpaceDiscr_ghost_rho_top_neumann
};


/*------------------------------------------------------------------------------
 top rhou
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhou_top_void(
											   double* rhou, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhou_top_wall(
											   double* rhou, 
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
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_top_inflow(
												 double* rhou, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_top_outflow(
												  double* rhou, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_top_periodic(
												   double* rhou, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhou[n + mobj] = rhou[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_top_neumann(
												  double* rhou, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhou_top[])(
											   double* rhou, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhou_top_void,
NodeSpaceDiscr_ghost_rhou_top_wall,
NodeSpaceDiscr_ghost_rhou_top_inflow,
NodeSpaceDiscr_ghost_rhou_top_outflow,
NodeSpaceDiscr_ghost_rhou_top_periodic,
NodeSpaceDiscr_ghost_rhou_top_neumann
};


/*------------------------------------------------------------------------------
 top rhov
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhov_top_void(
											   double* rhov, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhov_top_wall(
											   double* rhov, 
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
				rhov[n + mobj] = - rhov[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_top_inflow(
												 double* rhov, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_top_outflow(
												  double* rhov, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_top_periodic(
												   double* rhov, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhov[n + mobj] = rhov[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_top_neumann(
												  double* rhov, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhov_top[])(
											   double* rhov, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhov_top_void,
NodeSpaceDiscr_ghost_rhov_top_wall,
NodeSpaceDiscr_ghost_rhov_top_inflow,
NodeSpaceDiscr_ghost_rhov_top_outflow,
NodeSpaceDiscr_ghost_rhov_top_periodic,
NodeSpaceDiscr_ghost_rhov_top_neumann
};


/*------------------------------------------------------------------------------
 top rhow
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhow_top_void(
											   double* rhow, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhow_top_wall(
											   double* rhow, 
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
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_top_inflow(
												 double* rhow, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_top_outflow(
												  double* rhow, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_top_periodic(
												   double* rhow, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhow[n + mobj] = rhow[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_top_neumann(
												  double* rhow, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhow_top[])(
											   double* rhow, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhow_top_void,
NodeSpaceDiscr_ghost_rhow_top_wall,
NodeSpaceDiscr_ghost_rhow_top_inflow,
NodeSpaceDiscr_ghost_rhow_top_outflow,
NodeSpaceDiscr_ghost_rhow_top_periodic,
NodeSpaceDiscr_ghost_rhow_top_neumann
};


/*------------------------------------------------------------------------------
 top rhoe
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoe_top_void(
											   double* rhoe, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhoe_top_wall(
											   double* rhoe, 
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
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_top_inflow(
												 double* rhoe, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_top_outflow(
												  double* rhoe, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_top_periodic(
												   double* rhoe, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoe[n + mobj] = rhoe[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_top_neumann(
												  double* rhoe, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoe_top[])(
											   double* rhoe, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhoe_top_void,
NodeSpaceDiscr_ghost_rhoe_top_wall,
NodeSpaceDiscr_ghost_rhoe_top_inflow,
NodeSpaceDiscr_ghost_rhoe_top_outflow,
NodeSpaceDiscr_ghost_rhoe_top_periodic,
NodeSpaceDiscr_ghost_rhoe_top_neumann
};


/*------------------------------------------------------------------------------
 top rhoY
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoY_top_void(
											   double* rhoY, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhoY_top_wall(
											   double* rhoY, 
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
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_top_inflow(
												 double* rhoY, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_top_outflow(
												  double* rhoY, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_top_periodic(
												   double* rhoY, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n  = l + i;
				rhoY[n + mobj] = rhoY[n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_top_neumann(
												  double* rhoY, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoY_top[])(
											   double* rhoY, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhoY_top_void,
NodeSpaceDiscr_ghost_rhoY_top_wall,
NodeSpaceDiscr_ghost_rhoY_top_inflow,
NodeSpaceDiscr_ghost_rhoY_top_outflow,
NodeSpaceDiscr_ghost_rhoY_top_periodic,
NodeSpaceDiscr_ghost_rhoY_top_neumann
};


/*------------------------------------------------------------------------------
 top rhoZ
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoZ_top_void(
											   double** rhoZ, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rhoZ_top_wall(
											   double** rhoZ, 
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
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_top_inflow(
												 double** rhoZ, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_top_outflow(
												  double** rhoZ, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_top_periodic(
												   double** rhoZ, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int icz = node->icz;
	const int bsrc = igy + 1;
	const int bobj = icy - igy;
	int i, j, k, l, n, msrc, mobj;
	
	assert(ig <= igy);
	
	for(j = 0; j < ig; j++) {
		msrc = (bsrc + j) * icx;
		mobj = (bobj + j) * icx;
		for(k = 0; k < icz; k++) {l = k * icx * icy;
			for(i = 0; i < icx; i++) {n = l + i;
				rhoZ[PRES][n + mobj] = rhoZ[PRES][n + msrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_top_neumann(
												  double** rhoZ, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoZ_top[])(
											   double** rhoZ, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rhoZ_top_void,
NodeSpaceDiscr_ghost_rhoZ_top_wall,
NodeSpaceDiscr_ghost_rhoZ_top_inflow,
NodeSpaceDiscr_ghost_rhoZ_top_outflow,
NodeSpaceDiscr_ghost_rhoZ_top_periodic,
NodeSpaceDiscr_ghost_rhoZ_top_neumann
};


/*------------------------------------------------------------------------------
 back rho
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rho_back_void(
											   double* rho, 
											   const NodeSpaceDiscr* node, 
											   const int ig) {}

static void NodeSpaceDiscr_ghost_rho_back_wall(
											   double* rho, 
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
				rho[n + lobj] = rho[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_back_inflow(
												 double* rho, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_back_outflow(
												  double* rho, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_back_periodic(
												   double* rho, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_back_neumann(
												  double* rho, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rho_back[])(
											   double* rho, 
											   const NodeSpaceDiscr* node, 
											   const int ig) = {
NodeSpaceDiscr_ghost_rho_back_void,
NodeSpaceDiscr_ghost_rho_back_wall,
NodeSpaceDiscr_ghost_rho_back_inflow,
NodeSpaceDiscr_ghost_rho_back_outflow,
NodeSpaceDiscr_ghost_rho_back_periodic,
NodeSpaceDiscr_ghost_rho_back_neumann
};


/*------------------------------------------------------------------------------
 back rhou
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhou_back_void(
												double* rhou, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhou_back_wall(
												double* rhou, 
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
				rhou[n + lobj] = rhou[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_back_inflow(
												  double* rhou, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_back_outflow(
												   double* rhou, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_back_periodic(
													double* rhou, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_back_neumann(
												   double* rhou, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhou_back[])(
												double* rhou, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhou_back_void,
NodeSpaceDiscr_ghost_rhou_back_wall,
NodeSpaceDiscr_ghost_rhou_back_inflow,
NodeSpaceDiscr_ghost_rhou_back_outflow,
NodeSpaceDiscr_ghost_rhou_back_periodic,
NodeSpaceDiscr_ghost_rhou_back_neumann
};


/*------------------------------------------------------------------------------
 back rhov
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhov_back_void(
												double* rhov, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhov_back_wall(
												double* rhov, 
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
				rhov[n + lobj] = rhov[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_back_inflow(
												  double* rhov, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_back_outflow(
												   double* rhov, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_back_periodic(
													double* rhov, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_back_neumann(
												   double* rhov, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhov_back[])(
												double* rhov, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhov_back_void,
NodeSpaceDiscr_ghost_rhov_back_wall,
NodeSpaceDiscr_ghost_rhov_back_inflow,
NodeSpaceDiscr_ghost_rhov_back_outflow,
NodeSpaceDiscr_ghost_rhov_back_periodic,
NodeSpaceDiscr_ghost_rhov_back_neumann
};

/*------------------------------------------------------------------------------
 back rhow
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhow_back_void(
												double* rhow, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhow_back_wall(
												double* rhow, 
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
				rhow[n + lobj] = - rhow[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_back_inflow(
												  double* rhow, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_back_outflow(
												   double* rhow, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_back_periodic(
													double* rhow, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_back_neumann(
												   double* rhow, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhow_back[])(
												double* rhow, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhow_back_void,
NodeSpaceDiscr_ghost_rhow_back_wall,
NodeSpaceDiscr_ghost_rhow_back_inflow,
NodeSpaceDiscr_ghost_rhow_back_outflow,
NodeSpaceDiscr_ghost_rhow_back_periodic,
NodeSpaceDiscr_ghost_rhow_back_neumann
};


/*------------------------------------------------------------------------------
 back rhoe
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoe_back_void(
												double* rhoe, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhoe_back_wall(
												double* rhoe, 
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
				rhoe[n + lobj] = rhoe[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_back_inflow(
												  double* rhoe, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_back_outflow(
												   double* rhoe, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_back_periodic(
													double* rhoe, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_back_neumann(
												   double* rhoe, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoe_back[])(
												double* rhoe, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhoe_back_void,
NodeSpaceDiscr_ghost_rhoe_back_wall,
NodeSpaceDiscr_ghost_rhoe_back_inflow,
NodeSpaceDiscr_ghost_rhoe_back_outflow,
NodeSpaceDiscr_ghost_rhoe_back_periodic,
NodeSpaceDiscr_ghost_rhoe_back_neumann
};


/*------------------------------------------------------------------------------
 back rhoY
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoY_back_void(
												double* rhoY, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhoY_back_wall(
												double* rhoY, 
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
				rhoY[n + lobj] = rhoY[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_back_inflow(
												  double* rhoY, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_back_outflow(
												   double* rhoY, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_back_periodic(
													double* rhoY, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_back_neumann(
												   double* rhoY, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoY_back[])(
												double* rhoY, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhoY_back_void,
NodeSpaceDiscr_ghost_rhoY_back_wall,
NodeSpaceDiscr_ghost_rhoY_back_inflow,
NodeSpaceDiscr_ghost_rhoY_back_outflow,
NodeSpaceDiscr_ghost_rhoY_back_periodic,
NodeSpaceDiscr_ghost_rhoY_back_neumann
};


/*------------------------------------------------------------------------------
 back rhoZ
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoZ_back_void(
												double** rhoZ, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rhoZ_back_wall(
												double** rhoZ, 
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
				rhoZ[PRES][n + lobj] = rhoZ[PRES][n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_back_inflow(
												  double** rhoZ, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_back_outflow(
												   double** rhoZ, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_back_periodic(
													double** rhoZ, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_back_neumann(
												   double** rhoZ, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoZ_back[])(
												double** rhoZ, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rhoZ_back_void,
NodeSpaceDiscr_ghost_rhoZ_back_wall,
NodeSpaceDiscr_ghost_rhoZ_back_inflow,
NodeSpaceDiscr_ghost_rhoZ_back_outflow,
NodeSpaceDiscr_ghost_rhoZ_back_periodic,
NodeSpaceDiscr_ghost_rhoZ_back_neumann
};


/*------------------------------------------------------------------------------
 front rho
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rho_front_void(
												double* rho, 
												const NodeSpaceDiscr* node, 
												const int ig) {}

static void NodeSpaceDiscr_ghost_rho_front_wall(
												double* rho, 
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
				rho[n + lobj] = rho[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rho_front_inflow(
												  double* rho, 
												  const NodeSpaceDiscr* node, 
												  const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_front_outflow(
												   double* rho, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_front_periodic(
													double* rho, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rho_front_neumann(
												   double* rho, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rho_front[])(
												double* rho, 
												const NodeSpaceDiscr* node, 
												const int ig) = {
NodeSpaceDiscr_ghost_rho_front_void,
NodeSpaceDiscr_ghost_rho_front_wall,
NodeSpaceDiscr_ghost_rho_front_inflow,
NodeSpaceDiscr_ghost_rho_front_outflow,
NodeSpaceDiscr_ghost_rho_front_periodic,
NodeSpaceDiscr_ghost_rho_front_neumann
};


/*------------------------------------------------------------------------------
 front rhou
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhou_front_void(
												 double* rhou, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhou_front_wall(
												 double* rhou, 
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
				rhou[n + lobj] = rhou[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhou_front_inflow(
												   double* rhou, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_front_outflow(
													double* rhou, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_front_periodic(
													 double* rhou, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhou_front_neumann(
													double* rhou, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhou_front[])(
												 double* rhou, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhou_front_void,
NodeSpaceDiscr_ghost_rhou_front_wall,
NodeSpaceDiscr_ghost_rhou_front_inflow,
NodeSpaceDiscr_ghost_rhou_front_outflow,
NodeSpaceDiscr_ghost_rhou_front_periodic,
NodeSpaceDiscr_ghost_rhou_front_neumann
};


/*------------------------------------------------------------------------------
 front rhov
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhov_front_void(
												 double* rhov, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhov_front_wall(
												 double* rhov, 
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
				rhov[n + lobj] = rhov[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhov_front_inflow(
												   double* rhov, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_front_outflow(
													double* rhov, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_front_periodic(
													 double* rhov, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhov_front_neumann(
													double* rhov, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhov_front[])(
												 double* rhov, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhov_front_void,
NodeSpaceDiscr_ghost_rhov_front_wall,
NodeSpaceDiscr_ghost_rhov_front_inflow,
NodeSpaceDiscr_ghost_rhov_front_outflow,
NodeSpaceDiscr_ghost_rhov_front_periodic,
NodeSpaceDiscr_ghost_rhov_front_neumann
};


/*------------------------------------------------------------------------------
 front rhow
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhow_front_void(
												 double* rhow, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhow_front_wall(
												 double* rhow, 
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
				rhow[n + lobj] = - rhow[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhow_front_inflow(
												   double* rhow, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_front_outflow(
													double* rhow, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_front_periodic(
													 double* rhow, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhow_front_neumann(
													double* rhow, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhow_front[])(
												 double* rhow, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhow_front_void,
NodeSpaceDiscr_ghost_rhow_front_wall,
NodeSpaceDiscr_ghost_rhow_front_inflow,
NodeSpaceDiscr_ghost_rhow_front_outflow,
NodeSpaceDiscr_ghost_rhow_front_periodic,
NodeSpaceDiscr_ghost_rhow_front_neumann
};


/*------------------------------------------------------------------------------
 front rhoe
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoe_front_void(
												 double* rhoe, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhoe_front_wall(
												 double* rhoe, 
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
				rhoe[n + lobj] = rhoe[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoe_front_inflow(
												   double* rhoe, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_front_outflow(
													double* rhoe, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_front_periodic(
													 double* rhoe, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoe_front_neumann(
													double* rhoe, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoe_front[])(
												 double* rhoe, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhoe_front_void,
NodeSpaceDiscr_ghost_rhoe_front_wall,
NodeSpaceDiscr_ghost_rhoe_front_inflow,
NodeSpaceDiscr_ghost_rhoe_front_outflow,
NodeSpaceDiscr_ghost_rhoe_front_periodic,
NodeSpaceDiscr_ghost_rhoe_front_neumann
};


/*------------------------------------------------------------------------------
 front rhoY
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoY_front_void(
												 double* rhoY, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhoY_front_wall(
												 double* rhoY, 
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
				rhoY[n + lobj] = rhoY[n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoY_front_inflow(
												   double* rhoY, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_front_outflow(
													double* rhoY, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_front_periodic(
													 double* rhoY, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoY_front_neumann(
													double* rhoY, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoY_front[])(
												 double* rhoY, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhoY_front_void,
NodeSpaceDiscr_ghost_rhoY_front_wall,
NodeSpaceDiscr_ghost_rhoY_front_inflow,
NodeSpaceDiscr_ghost_rhoY_front_outflow,
NodeSpaceDiscr_ghost_rhoY_front_periodic,
NodeSpaceDiscr_ghost_rhoY_front_neumann
};


/*------------------------------------------------------------------------------
 front rhoZ
 ------------------------------------------------------------------------------*/
static void NodeSpaceDiscr_ghost_rhoZ_front_void(
												 double** rhoZ, 
												 const NodeSpaceDiscr* node, 
												 const int ig) {}

static void NodeSpaceDiscr_ghost_rhoZ_front_wall(
												 double** rhoZ, 
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
				rhoZ[PRES][n + lobj] = rhoZ[PRES][n + lsrc];
			}
		}
	}
}

static void NodeSpaceDiscr_ghost_rhoZ_front_inflow(
												   double** rhoZ, 
												   const NodeSpaceDiscr* node, 
												   const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_front_outflow(
													double** rhoZ, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_front_periodic(
													 double** rhoZ, 
													 const NodeSpaceDiscr* node, 
													 const int ig) {
	ERROR("function not available");
}

static void NodeSpaceDiscr_ghost_rhoZ_front_neumann(
													double** rhoZ, 
													const NodeSpaceDiscr* node, 
													const int ig) {
	ERROR("function not available");
}

static void (*NodeSpaceDiscr_ghost_rhoZ_front[])(
												 double** rhoZ, 
												 const NodeSpaceDiscr* node, 
												 const int ig) = {
NodeSpaceDiscr_ghost_rhoZ_front_void,
NodeSpaceDiscr_ghost_rhoZ_front_wall,
NodeSpaceDiscr_ghost_rhoZ_front_inflow,
NodeSpaceDiscr_ghost_rhoZ_front_outflow,
NodeSpaceDiscr_ghost_rhoZ_front_periodic,
NodeSpaceDiscr_ghost_rhoZ_front_neumann
};



void NodeSpaceDiscr_ghost_rho(
							  double* rho, 
							  const NodeSpaceDiscr* node, 
							  const int ig) {
	
	(*NodeSpaceDiscr_ghost_rho_left[node->left])(rho, node, ig);
	(*NodeSpaceDiscr_ghost_rho_right[node->right])(rho, node, ig);
	(*NodeSpaceDiscr_ghost_rho_bottom[node->bottom])(rho, node, ig);
	(*NodeSpaceDiscr_ghost_rho_top[node->top])(rho, node, ig);
	(*NodeSpaceDiscr_ghost_rho_back[node->back])(rho, node, ig);
	(*NodeSpaceDiscr_ghost_rho_front[node->front])(rho, node, ig);
}


void NodeSpaceDiscr_ghost_rhou(
							   double* rhou, 
							   const NodeSpaceDiscr* node, 
							   const int ig) {
	
	(*NodeSpaceDiscr_ghost_rhou_left[node->left])(rhou, node, ig);
	(*NodeSpaceDiscr_ghost_rhou_right[node->right])(rhou, node, ig);
	(*NodeSpaceDiscr_ghost_rhou_bottom[node->bottom])(rhou, node, ig);
	(*NodeSpaceDiscr_ghost_rhou_top[node->top])(rhou, node, ig);
	(*NodeSpaceDiscr_ghost_rhou_back[node->back])(rhou, node, ig);
	(*NodeSpaceDiscr_ghost_rhou_front[node->front])(rhou, node, ig);
}


void NodeSpaceDiscr_ghost_rhov(
							   double* rhov, 
							   const NodeSpaceDiscr* node, 
							   const int ig) {
	
	(*NodeSpaceDiscr_ghost_rhov_left[node->left])(rhov, node, ig);
	(*NodeSpaceDiscr_ghost_rhov_right[node->right])(rhov, node, ig);
	(*NodeSpaceDiscr_ghost_rhov_bottom[node->bottom])(rhov, node, ig);
	(*NodeSpaceDiscr_ghost_rhov_top[node->top])(rhov, node, ig);
	(*NodeSpaceDiscr_ghost_rhov_back[node->back])(rhov, node, ig);
	(*NodeSpaceDiscr_ghost_rhov_front[node->front])(rhov, node, ig);
}


void NodeSpaceDiscr_ghost_rhow(
							   double* rhow, 
							   const NodeSpaceDiscr* node, 
							   const int ig) {
	
	(*NodeSpaceDiscr_ghost_rhow_left[node->left])(rhow, node, ig);
	(*NodeSpaceDiscr_ghost_rhow_right[node->right])(rhow, node, ig);
	(*NodeSpaceDiscr_ghost_rhow_bottom[node->bottom])(rhow, node, ig);
	(*NodeSpaceDiscr_ghost_rhow_top[node->top])(rhow, node, ig);
	(*NodeSpaceDiscr_ghost_rhow_back[node->back])(rhow, node, ig);
	(*NodeSpaceDiscr_ghost_rhow_front[node->front])(rhow, node, ig);
}


void NodeSpaceDiscr_ghost_rhoe(
							   double* rhoe, 
							   const NodeSpaceDiscr* node, 
							   const int ig) {
	
	(*NodeSpaceDiscr_ghost_rhoe_left[node->left])(rhoe, node, ig);
	(*NodeSpaceDiscr_ghost_rhoe_right[node->right])(rhoe, node, ig);
	(*NodeSpaceDiscr_ghost_rhoe_bottom[node->bottom])(rhoe, node, ig);
	(*NodeSpaceDiscr_ghost_rhoe_top[node->top])(rhoe, node, ig);
	(*NodeSpaceDiscr_ghost_rhoe_back[node->back])(rhoe, node, ig);
	(*NodeSpaceDiscr_ghost_rhoe_front[node->front])(rhoe, node, ig);
}


void NodeSpaceDiscr_ghost_rhoY(
							   double* rhoY, 
							   const NodeSpaceDiscr* node, 
							   const int ig) {
	
	(*NodeSpaceDiscr_ghost_rhoY_left[node->left])(rhoY, node, ig);
	(*NodeSpaceDiscr_ghost_rhoY_right[node->right])(rhoY, node, ig);
	(*NodeSpaceDiscr_ghost_rhoY_bottom[node->bottom])(rhoY, node, ig);
	(*NodeSpaceDiscr_ghost_rhoY_top[node->top])(rhoY, node, ig);
	(*NodeSpaceDiscr_ghost_rhoY_back[node->back])(rhoY, node, ig);
	(*NodeSpaceDiscr_ghost_rhoY_front[node->front])(rhoY, node, ig);
}


void NodeSpaceDiscr_ghost_rhoZ(
							   double** rhoZ, 
							   const NodeSpaceDiscr* node, 
							   const int ig) {
	
	(*NodeSpaceDiscr_ghost_rhoZ_left[node->left])(rhoZ, node, ig);
	(*NodeSpaceDiscr_ghost_rhoZ_right[node->right])(rhoZ, node, ig);
	(*NodeSpaceDiscr_ghost_rhoZ_bottom[node->bottom])(rhoZ, node, ig);
	(*NodeSpaceDiscr_ghost_rhoZ_top[node->top])(rhoZ, node, ig);
	(*NodeSpaceDiscr_ghost_rhoZ_back[node->back])(rhoZ, node, ig);
	(*NodeSpaceDiscr_ghost_rhoZ_front[node->front])(rhoZ, node, ig);
}


void NodeSpaceDiscr_ghost(
						  ConsVars* Sol, 
						  const NodeSpaceDiscr* node, 
						  const int ig) {
	
	NodeSpaceDiscr_ghost_rho(Sol->rho, node, ig);
	NodeSpaceDiscr_ghost_rhou(Sol->rhou, node, ig);
	NodeSpaceDiscr_ghost_rhov(Sol->rhov, node, ig);
	NodeSpaceDiscr_ghost_rhow(Sol->rhow, node, ig);
	NodeSpaceDiscr_ghost_rhoe(Sol->rhoe, node, ig);
	NodeSpaceDiscr_ghost_rhoY(Sol->rhoY, node, ig);
	NodeSpaceDiscr_ghost_rhoZ(Sol->rhoZ, node, ig);
}













/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: discretization.c,v $
 Revision 1.2  1998/03/07 09:56:45  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:33  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
