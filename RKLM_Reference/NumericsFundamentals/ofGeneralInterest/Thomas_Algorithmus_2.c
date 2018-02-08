#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "ThomasAlgorithmus.h"
#include "math_own.h"

/*********************************************************************
 tridiag(): solution of a tridiagonal system of equations a * x = b of size n
 a = tridiagonal matrix, x = solution vector, b = right hand side
 structure of a: [a(0,0), a(0,1), a(0,2), 
 a(2,1), a(2,2), a(2,3),
 ... a(i,i-1), a(i,i), a(i,i+1) ...
 a(n-1,n-3), a(n-1,n-2), a(n-1,n-1)]
 diagonal elements must be non zero!!!
 *********************************************************************/

void tridiag(double* a, double* x, double* b, int n)
{    
    double *l = (double*)malloc(n * sizeof(double));
    double *r = (double*)malloc(n * sizeof(double));
    double *m = (double*)malloc(n * sizeof(double));
    double *y = (double*)malloc(n * sizeof(double));
    
    int i, ii;
    
    /****************************/
    /* LR reduction of matrix a */
    /****************************/
    m[0]   = a[0];
    r[0]   = a[1];
    r[n-1] = a[2];
    l[1]   = a[3] / m[0];
    m[1]   = a[4] - l[1] * r[0];
    r[1]   = a[5] - l[1] * r[n-1];
    for(i = 2; i < n-1; i++)
    {
        ii = i * 3 + 1;
        l[i] = a[ii-1] / m[i-1];
        m[i] = a[ii] - l[i] * r[i-1];
        r[i] = a[ii+1];
    }
    ii = (n-1) * 3 + 1;
    l[0]   = a[ii-1] / m[n-3];
    l[n-1] = (a[ii] - l[0] * r[n-3]) / m[n-2];
    m[n-1] = a[ii+1] - l[n-1] * r[n-2];
    /*************************************************/
    /* intermediate vector from L y = b with y = R x */
    /*************************************************/
    y[0] = b[0];
    for(i = 1; i < n; i++)
        y[i] = b[i] - l[i] * y[i-1];
    y[n-1] -= l[0] * y[n-3];
    /*********************************/
    /* solution vector from R x = y; */
    /*********************************/
    x[n-1] = y[n-1] / m[n-1];
    for(i = n-2; i >= 0; i--)
        x[i] = (y[i] - r[i] * x[i+1]) / m[i];
    x[0] -= r[n-1] * x[2] / m[0];
    
    free (l);
    free (r);
    free (m);
    free (y);    
    
} /** tridiag() **/

/* ********************************************************************************** */

void Thomas_Algorithm(double* x,  double* rhs, double*  upper, double* diago, double* lower, int size)
{
    double *l_lower = (double*)malloc(size * sizeof(double));
    double *r_diago = (double*)malloc(size * sizeof(double));
    double *r_upper = (double*)malloc(size * sizeof(double));
    double *y       = (double*)malloc(size * sizeof(double));
    
    int i;
    /************************************************************************************/
    /* LR reduction of the tri-diagonal matrix   made up by (lower  ,  diago ,  upper ) */
    /* where L, R, are the tri-diagonal matrices made up by (l_lower,   1    ,    0   ) */
    /* and                                                  (  0    , r_diago, r_upper) */
    /************************************************************************************/
    r_diago[0]  = diago[0];
    r_upper[0]  = upper[0];
    
    for(i = 1; i < size; i++)
    {
        l_lower[i] = lower[i] / r_diago[i-1];
        r_diago[i] = diago[i] - l_lower[i] * r_upper[i-1];
        r_upper[i] = upper[i];
    }
    
    /*************************************************/
    /* intermediate vector from L y = b with y = R x */
    /*************************************************/
    y[0] = rhs[0];
    
    for(i = 1; i < size; i++)
    {
        y[i] = rhs[i] - l_lower[i] * y[i-1];
    }
    
    /*********************************/
    /* solution vector from R x = y; */
    /*********************************/
    x[size-1] = y[size-1] / r_diago[size-1];
    
    for(i = size-2; i >= 0; i--)
    {
        x[i] = (y[i] - r_upper[i] * x[i+1]) / r_diago[i];
    }
    
    free (l_lower);
    free (r_diago);
    free (r_upper);
    free (y      );
    
} /** Thomas_Algorithm() **/

/* ********************************************************************************** */

void LowerTripleTriangleSolve(double* x,  double* rhs,  double* diago, double* lower1, double* lower2, int size)
{
    /* solves a linear system with lower triangular structure */
    int i;
    
    x[0] = rhs[0] / diago[0];
    
    for(i=1; i<size; i++)
    {
        x[i] = (rhs[i] - lower1[i]*x[i-1]) / diago[i];
    }
    
    for(i=2; i<size; i++)
    {
        x[i] -= lower2[i]*x[i-2] / diago[i];
    }
}

/* ********************************************************************************** */
/* Templeton's Algorithm_4   from  (1975)                                             */
/* ********************************************************************************** */

void Algorithm_4(double* x, double* rhs, double lambda, int size)
{
    double *l_lower = (double*)malloc(size * sizeof(double));
    double *r_diago = (double*)malloc(size * sizeof(double));
    double *r_upper = (double*)malloc(size * sizeof(double));
    double *y       = (double*)malloc(size * sizeof(double));
    
    double *z = y;
    
    double alpha, sigma;
    int i;
    
    /************************************************************************************/
    /* preparation to get  x[0]  explicitly                                             */
    /************************************************************************************/
    x[0] = 0.0;
    
    assert(lambda > 2.0);
    alpha = 0.5* (-lambda + sqrt(lambda*lambda - 4));
    sigma = 1.0 / ( (alpha + pow(alpha,size-1)) + lambda * (1.0 + pow(alpha,size)));
    
    /* first row of the inverse matrix */
    for (int i=0; i<size; i++) {
        z[i] = sigma * (pow(alpha,i) + pow(alpha,size-i));
    }
    
    /* there is a more efficient version for this loop in Templeton (1975) */
    /* with first row of inverse matrix z[] determined, we have            */
    for (i = 0; i<size; i++) {
        x[0] += z[i]*rhs[i];
    }
    
    /************************************************************************************/
    /* Thomas algorithm for the rest case      lower = upper = 1.0; diago = lambda      */
    /************************************************************************************/
    
    /************************************************************************************/
    /* LR reduction of the tri-diagonal matrix   made up by (lower  ,  diago ,  upper ) */
    /* where L, R, are the tri-diagonal matrices made up by (l_lower,   1    ,    0   ) */
    /* and                                                  (  0    , r_diago, r_upper) */
    /************************************************************************************/
    r_diago[1]  = lambda;
    r_upper[1]  = 1.0;
    
    for(i = 2; i < size; i++)
    {
        l_lower[i] = 1.0 / r_diago[i-1];
        r_diago[i] = lambda - l_lower[i] * r_upper[i-1];
        r_upper[i] = 1.0;
    }
    
    /*************************************************/
    /* intermediate vector from L y = b with y = R x */
    /*************************************************/
    y[1] = rhs[1] - x[0];
    
    for(i = 2; i < size-1; i++)
    {
        y[i] = rhs[i] - l_lower[i] * y[i-1];
    }
    
    y[size-1] = (rhs[size-1]-x[0]) - l_lower[size-1] * y[size-2];
    
    /*********************************/
    /* solution vector from R x = y; */
    /*********************************/
    x[size-1] = y[size-1] / r_diago[size-1];
    
    for(i = size-2; i >= 1; i--)
    {
        x[i] = (y[i] - r_upper[i] * x[i+1]) / r_diago[i];
    }
    
    free (l_lower);
    free (r_diago);
    free (r_upper);
    free (y      );
    
#if 0
    /* test */
    {
        double delta = 0.0;
        double rmax  = 0.0;
        
        delta = MAX_own(delta, fabs(rhs[0] - (x[size-1] + lambda*x[0] + x[1])));
        rmax  = MAX_own(rmax, fabs(rhs[0]));
        
        for (int i=1; i<size-1; i++) {
            delta = MAX_own(delta, fabs(rhs[i] - (x[i-1] + lambda*x[i] + x[i+1])));
            rmax  = MAX_own(rmax, fabs(rhs[i]));
        }
        
        delta = MAX_own(delta, fabs(rhs[size-1] - (x[size-2] + lambda*x[size-1] + x[0])));
        rmax  = MAX_own(rmax, fabs(rhs[size-1]));

        printf("Algorithm 4:  delta = %e, rmax = %e\n", delta, rmax);
    }
#endif
    
} /** Algorithm_4() **/
