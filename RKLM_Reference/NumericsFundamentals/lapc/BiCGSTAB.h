/*******************************************************************************
 File:   BiCGSTAB.h
 Author: Nicola
 Date:   Sat Mar  7 19:26:55 CET 1998
 
 Notice: C version of Andreas Meister's C++ BiCGSTAB class
 *******************************************************************************/
#ifndef BiCGSTAB_H
#define BiCGSTAB_H

#include "Common.h"
#include "ProjectionType.h"
#include "BiCGSTAB.h"
#include "kgrid.h"
#include "variable.h"
#include "mpv.h"

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
typedef struct {
    
    int size;
    
    double* r_0;
    double* r_j;
    double* p_j;
    double* v_j;
    double* s_j;
    double* t_j;
    double* help_vec;
    double* help_vec2;
    
#ifdef SOLVER_CR2
    double* Lr_0;
#endif /* SOLVER_CR2 */
    
    double precision;
    double local_precision;
    int max_iterations;
    int actual_iterations;
    int output_period;
    
} BiCGSTABData;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
BiCGSTABData* BiCGSTABData_new(
                               const int size,
                               const double precision,
                               const double local_precision,
                               const int max_iterations,
                               const int output_period);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void BiCGSTABData_free(BiCGSTABData* var);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double SOLVER(
              BiCGSTABData* data,
              const ElemSpaceDiscr* elem,
              const NodeSpaceDiscr* node,
              const double* hplus[3],
              const double* hcenter,
              const double* hgrav,
              const ConsVars* Sol,
              const MPV* mpv,
              const double dt,
              const double theta,
              double* rhs,
              double* p2);
#endif



/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/



