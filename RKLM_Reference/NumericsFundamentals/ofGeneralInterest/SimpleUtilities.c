#include <assert.h> 
#include "SimpleUtilities.h"

/* ================================================================================== */
double Invert_3x3_Matrix(double Tinv[3][3], double T[3][3])
{
	/* returns the determinant detT of T[][] and computes the inverse of it
     up to a factor of 1/detT                                              */
	double detT;
	int i, ii; 
	static int next[3] = {1,2,0}, prev[3] = {2,0,1};
	
	detT =   T[0][0] * (T[1][1]*T[2][2] - T[1][2]*T[2][1]) 
	- T[1][0] * (T[0][1]*T[2][2] - T[0][2]*T[2][1]) 
	+ T[2][0] * (T[0][1]*T[1][2] - T[0][2]*T[1][1]);
	
	for(i=0; i<3; i++)
    {
		for(ii=0; ii<3; ii++)
        {
			Tinv[ii][i] =  T[next[i]][next[ii]] * T[prev[i]][prev[ii]] 
			- T[next[i]][prev[ii]] * T[prev[i]][next[ii]];
        }
    }
	
	return(detT);
}

/* ================================================================================== */
void Multiply_3x3_Matrices(double AB[3][3], double A[3][3], double B[3][3])
{
	/* returns the determinant detT of T[][] and computes the inverse of it
     up to a factor of 1/detT                                              */
	int i, j, k; 
	
	for(i=0; i<3; i++)
    {
		for(j=0; j<3; j++)
        {
			for(k=0, AB[i][j] = 0; k<3; k++)
			{
				AB[i][j] += A[i][k] * B[k][j];
			}
        }
    }
}

/* ================================================================================== */

int power_of(int n, int base, int exponent)
{
    if(n==1)
    { 
        return(exponent);
    }
    else if(n%base == 1)
    {
        return(0);
    }
    else
    {
        return(power_of(n/base,base,exponent+1));
    }
}

/* ================================================================================== */

int integer_power_detailed(int result, int n, int exponent, int counter);

int integer_power(int n, int exponent)
{
    return(integer_power_detailed(n, n, exponent, 0));
}

/* ================================================================================== */

int integer_power_detailed(int result, int n, int exponent, int counter)
{
    assert(n>=0);
    assert(exponent>=0);
    
    if(exponent==0)
    { 
        return(1);
    }
    else if(exponent == 1)
    {
        return(result);
    }
    else
    {
        return(integer_power_detailed(n*result, n, exponent-1, counter+1));
    }
}


