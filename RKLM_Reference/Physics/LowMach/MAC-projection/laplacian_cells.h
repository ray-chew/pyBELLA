/*******************************************************************************
 File:   laplacian_cells.h
 Author: Nicola                        Rupert
 Date:   Thu Mar 12 08:20:38 CET 1998  Feb. 2004
 *******************************************************************************/
#ifndef LAPLACIAN_CELLS_H
#define LAPLACIAN_CELLS_H

#include "kgrid.h"
#include "variable.h"
#include "mpv.h"

/*------------------------------------------------------------------------------
 nine points (2d) enthalpy weighted discrete Laplacian based on bilinear 
 pressure distributions within the rectangles spanned by cell centers
 ------------------------------------------------------------------------------*/

double precon_c_prepare(
                      const NodeSpaceDiscr* node,
                      const ElemSpaceDiscr* elem,
                      const double* hplus[3],
                      const double* hcenter);

void precon_c_apply(
                    double* vec_out,
                    const double* vec_in,
                    const ElemSpaceDiscr *elem);

void precon_c_invert(
                     double* vec_out,
                     const double* vec_in,
                     const ElemSpaceDiscr *elem);

void EnthalpyWeightedLap_bilinear_p(
									const ElemSpaceDiscr* elem, 
                                    const NodeSpaceDiscr* node,
									const double* p,
									const double* hplus[3],
									const double* hcenter,
									const MPV* mpv, 
									const double dt,
									double* lap);

#endif /* LAPLACIAN_CELLS_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/





