// אם ירצה השם

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define sqr(x) ((x)*(x))



static double 
l2norm(unsigned int dim, double *p1, double *p2);

double*
WAM(double *points, unsigned int obsCount, unsigned int dim);

double 
sumRow(double *matrix, unsigned int cols, unsigned int rowIndex);

double* 
DDG(double *WAM, unsigned int obsCount);

void 
DHalf(double *DDG, unsigned int obsCount);

double* 
prettyDDG(double* DDG, unsigned int obsCount);

double* // HEZY 
laplacian(double* WAM, double* DHalf, unsigned int obsCount);

void 
printMatrix(double* matrix, unsigned int rows, unsigned int cols);