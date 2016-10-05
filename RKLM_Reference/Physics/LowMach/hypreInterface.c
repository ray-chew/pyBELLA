#include "Common.h"
#include "ProjectionType.h"
#include "userdata.h"
#include "enum_bdry.h"

#ifdef SOLVER_2_HYPRE_RUPE

/* #define HYPRE_PRINT */

#include "hypreInterface.h"

/* ======================================================================== */

void initialize_Hypre_Axb(
                          const int nx,
                          const int ny,
                          const enum Boolean cellPoisson,
                          HYPRE_SStructMatrix *A,
                          HYPRE_SStructVector *x,
                          HYPRE_SStructVector *b )
{
    extern User_Data ud;
    
    const int object_type = HYPRE_PARCSR;
    const int ndim = 2, nparts = 1, nvars = 1, part = 0;
    
    HYPRE_SStructGrid grid;
    HYPRE_SStructGridCreate(MPI_COMM_WORLD, ndim, nparts, &grid);
    
    if (cellPoisson == CORRECT)
    {
        int extents[2][2] = {{0,0}, {nx-1, ny-1}};
        HYPRE_SStructVariable vartypes[] = {HYPRE_SSTRUCT_VARIABLE_CELL};
        // create grid
        HYPRE_SStructGridSetExtents(grid, part, extents[0], extents[1]);
        HYPRE_SStructGridSetVariables(grid, part, nvars, vartypes);
        
#ifdef PERIODIC_X
        int map[] = {0,1};
        int ghostLeft[][2]  = { { -1   , 0 },{   -1, ny-1} };
        int ghostRight[][2] = { { nx   , 0 },{   nx, ny-1} };
        int innerLeft[][2]  = { { 0    , 0 },{    0, ny-1} };
        int innerRight[][2] = { { nx-1 , 0 },{ nx-1, ny-1} };
        int periodic[] = {nx,0};
        HYPRE_SStructGridSetPeriodic(grid, part, periodic);
#endif
    }
    else
    {
        int extents[2][2] = {{1,1}, {nx-1, ny-1}};
        HYPRE_SStructVariable vartypes[] = {HYPRE_SSTRUCT_VARIABLE_NODE};
        // create grid
        HYPRE_SStructGridSetExtents(grid, part, extents[0], extents[1]);
        HYPRE_SStructGridSetVariables(grid, part, nvars, vartypes);
        
        if ((ud.bdrytype_min[0] == PERIODIC) || (ud.bdrytype_min[1] == PERIODIC)) {
            int periodic[] = {(ud.bdrytype_min[0] == PERIODIC)*(nx-1), (ud.bdrytype_min[1] == PERIODIC)*(ny-1)};
            HYPRE_SStructGridSetPeriodic(grid, part, periodic);
        }

    }
    
    HYPRE_SStructGridAssemble(grid);
    
    // create stencil
    HYPRE_SStructStencil stencil;
    int size = 9;
    HYPRE_SStructStencilCreate(ndim, size, &stencil);
    int offsets[][2] = {{0,0}, {-1,0},{1,0},{0,-1},{0,1}, {-1,-1},{1,-1},{-1,1},{1,1}};
    for (int k = 0; k < size; k++)
        HYPRE_SStructStencilSetEntry(stencil, k, offsets[k], 0);
    
    // create graph
    HYPRE_SStructGraph graph;
    HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);
    if ( object_type != HYPRE_SSTRUCT )
    {
        HYPRE_SStructGraphSetObjectType(graph, object_type);
    }
    HYPRE_SStructGraphSetStencil(graph, part, 0, stencil);
    HYPRE_SStructGraphAssemble(graph);
    
    // create matrix
    HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, A);
    if ( object_type != HYPRE_SSTRUCT )
    {
        HYPRE_SStructMatrixSetObjectType(*A, object_type);
    }
    HYPRE_SStructMatrixInitialize(*A);
    
    // create  vectors
    HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, b);
    HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, x);
    if ( object_type != HYPRE_SSTRUCT )
    {
        HYPRE_SStructVectorSetObjectType(*b, object_type);
    }
    
    if ( object_type != HYPRE_SSTRUCT )
    {
        HYPRE_SStructVectorSetObjectType(*x, object_type);
    }
    HYPRE_SStructVectorInitialize(*b);
    HYPRE_SStructVectorInitialize(*x);
    
}

/* ======================================================================== */

void solve_with_Hypre(HYPRE_SStructMatrix A,
                      HYPRE_SStructVector x,
                      HYPRE_SStructVector b)
{
    extern User_Data ud;
    
    const int object_type = HYPRE_PARCSR;
    const int use_Boomer  = 1;
    
    HYPRE_SStructMatrixAssemble(A);
    HYPRE_SStructVectorAssemble(b);
    HYPRE_SStructVectorAssemble(x);
    
#ifdef HYPRE_PRINT
    char fname[100];
    sprintf(fname, "%s/Tests/HYPRE_A", ud.file_name);
    HYPRE_SStructMatrixPrint(fname, A, 0);
    sprintf(fname, "%s/Tests/HYPRE_b", ud.file_name);
    HYPRE_SStructVectorPrint(fname, b, 0);
#endif
    
    HYPRE_ParCSRMatrix   par_A;
    HYPRE_ParVector      par_b;
    HYPRE_ParVector      par_x;
    
    if (object_type == HYPRE_PARCSR)
    {
        HYPRE_SStructMatrixGetObject(A, (void **) &par_A);
        HYPRE_SStructVectorGetObject(b, (void **) &par_b);
        HYPRE_SStructVectorGetObject(x, (void **) &par_x);
    }
    
    HYPRE_Solver          par_solver;
    HYPRE_Solver          par_precond;

/* #define GMRES */
#ifdef GMRES
    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &par_solver);
    HYPRE_GMRESSetKDim(par_solver, 10);
    HYPRE_GMRESSetMaxIter(par_solver, 100);
    HYPRE_GMRESSetTol(par_solver, 1.0e-10);
    HYPRE_GMRESSetPrintLevel(par_solver, 3);
    HYPRE_GMRESSetLogging(par_solver, 1);
    
    if (use_Boomer)
    {
        /* use BoomerAMG as preconditioner */
        HYPRE_BoomerAMGCreate(&par_precond); 
        //HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
        //HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.8);
        HYPRE_BoomerAMGSetTol(par_precond, 0.0);
        HYPRE_BoomerAMGSetPrintLevel(par_precond, 0);
        HYPRE_BoomerAMGSetPrintFileName(par_precond, "sstruct.out.log");
        HYPRE_BoomerAMGSetMaxIter(par_precond, 5);
        HYPRE_BoomerAMGSetCycleRelaxType(par_precond, 3, 3);
        
        HYPRE_GMRESSetPrecond( par_solver,
                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                              par_precond);
    }
    
    HYPRE_GMRESSetup( par_solver, (HYPRE_Matrix) par_A,
                     (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);
    
    HYPRE_GMRESSolve( par_solver, (HYPRE_Matrix) par_A,
                     (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);
    
    //HYPRE_GMRESGetNumIterations( par_solver, &num_iterations);
    //HYPRE_GMRESGetFinalRelativeResidualNorm( par_solver, &final_res_norm);
    //printf("\n");
    //printf("Iterations = %d\n", num_iterations);
    //printf("Final Relative Residual Norm = %e\n", final_res_norm);
    //printf("\n");
    
    HYPRE_ParCSRGMRESDestroy(par_solver);
#else /* GMRES */
    HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &par_solver);
    // HYPRE_BiCGSTABSetKDim(par_solver, 10);
    HYPRE_BiCGSTABSetMaxIter(par_solver, 500);
    HYPRE_BiCGSTABSetTol(par_solver, 1.0e-6);
    HYPRE_BiCGSTABSetPrintLevel(par_solver, 3);
    HYPRE_BiCGSTABSetLogging(par_solver, 1);
    
    if (use_Boomer)
    {
        /* use BoomerAMG as preconditioner */
        HYPRE_BoomerAMGCreate(&par_precond); 
        //HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
        //HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.8);
        HYPRE_BoomerAMGSetTol(par_precond, 0.0);
        HYPRE_BoomerAMGSetPrintLevel(par_precond, 0);
        HYPRE_BoomerAMGSetPrintFileName(par_precond, "sstruct.out.log");
        HYPRE_BoomerAMGSetMaxIter(par_precond, 5);
        HYPRE_BoomerAMGSetCycleRelaxType(par_precond, 3, 3);
        
        HYPRE_BiCGSTABSetPrecond( par_solver,
                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                              par_precond);
    }
    
    HYPRE_BiCGSTABSetup( par_solver, (HYPRE_Matrix) par_A,
                     (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);
    
    HYPRE_BiCGSTABSolve( par_solver, (HYPRE_Matrix) par_A,
                     (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);
    
    //HYPRE_BiCGSTABGetNumIterations( par_solver, &num_iterations);
    //HYPRE_BiCGSTABGetFinalRelativeResidualNorm( par_solver, &final_res_norm);
    //printf("\n");
    //printf("Iterations = %d\n", num_iterations);
    //printf("Final Relative Residual Norm = %e\n", final_res_norm);
    //printf("\n");
    
    HYPRE_ParCSRBiCGSTABDestroy(par_solver);
#endif
    
    if (use_Boomer) {  
        HYPRE_BoomerAMGDestroy(par_precond);
    }
}

/* ======================================================================== */

void free_Hypre( HYPRE_IJMatrix hA,
                HYPRE_IJVector hx,
                HYPRE_IJVector hb )
{
    HYPRE_IJMatrixDestroy(hA);
    HYPRE_IJVectorDestroy(hb);
    HYPRE_IJVectorDestroy(hx);
}

#endif /* SOLVER_2_HYPRE_RUPE */