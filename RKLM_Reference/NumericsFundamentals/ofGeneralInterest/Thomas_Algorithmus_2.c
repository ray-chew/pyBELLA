#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "ThomasAlgorithmus.h"
#include "math_own.h"

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
