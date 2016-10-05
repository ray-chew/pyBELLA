#ifndef SIMPLEUTILITIES_H
#define SIMPLEUTILITIES_H
double Invert_3x3_Matrix(double Tinv[3][3], double T[3][3]);
void Multiply_3x3_Matrices(double AB[3][3], double A[3][3], double B[3][3]);

int power_of_two(int n, int base, int exponent);
int integer_power(int n, int exponent);
#endif
