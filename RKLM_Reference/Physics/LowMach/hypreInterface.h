#include "Common.h"

#ifdef SOLVER_2_HYPRE_RUPE

#ifndef HYPRE_INTERFACE_H
#define HYPRE_INTERFACE_H

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "krylov.h"
#include "enumerator.h"

void initialize_Hypre_Axb(const int nx,
                          const int ny,
                          const enum Boolean cellPoisson,
                          HYPRE_SStructMatrix *A,
                          HYPRE_SStructVector *x,
                          HYPRE_SStructVector *b );


void solve_with_Hypre(HYPRE_SStructMatrix A,
                      HYPRE_SStructVector x,
                      HYPRE_SStructVector b);


#endif

#endif /* SOLVER_2_HYPRE_RUPE */
