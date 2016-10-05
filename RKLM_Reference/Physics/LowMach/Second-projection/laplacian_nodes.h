/*******************************************************************************
 File:   laplacian_nodes.h
 Author: Nicola                        Rupert
 Date:   Thu Mar 12 08:20:38 CET 1998  March. 2004
 *******************************************************************************/
#ifndef LAPLACIAN_NODES_H
#define LAPLACIAN_NODES_H

#include "Common.h"
#include "kgrid.h"
#include "ProjectionType.h"
#include "mpv.h"

void precon_prepare(
                    const NodeSpaceDiscr* node,
                    const ElemSpaceDiscr* elem,
                    const double* hplus[3],
                    const double* wcenter,
                    const double* wgravb,
                    const int x_periodic,
                    const int y_periodic,
                    const int z_periodic);

void precon_apply(
                  double* vec_out,
                  const double* vec_in,
                  const NodeSpaceDiscr *node);

void precon_invert(
                   double* vec_out,
                   const double* vec_in,
                   const NodeSpaceDiscr *node);

void EnthalpyWeightedLap_Node_bilinear_p_scatter(
												 const NodeSpaceDiscr* node, 
												 const ElemSpaceDiscr* elem, 
												 const double* p,
												 const double* wplus[3],
												 const double* wcenter,
												 const double* wgrav,
												 const int x_periodic,
												 const int y_periodic,
												 const int z_periodic,
												 double* lap);

#endif /* LAPLACIAN_NODES_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/





