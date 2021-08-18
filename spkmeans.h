// אם ירצה השם

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

#define sqr(x) ((x)*(x))
#define EPSILON 0.001


static double 
l2norm(uint32_t dim, double *p1, double *p2);

double*
WAM(double *points, uint32_t obsCount, uint32_t dim);

double 
sumRow(double *matrix, uint32_t cols, uint32_t rowIndex);

double* 
DDG(double *WAM, uint32_t obsCount);

void 
DHalf(double *DDG, uint32_t obsCount);

double* 
prettyDDG(double* DDG, uint32_t obsCount);

double* // HEZY 
laplacian(double* WAM, double* DHalf, uint32_t obsCount);

void 
multiplyMatrices(double *A, double *B, double *C,
                 uint32_t m, uint32_t k, uint32_t n);

void 
maxItem(double *matrix, int rows, int *row, int *col);

void 
Jacobi(double *matrix, uint32_t rows, double old);

int 
argmax(double *eigenArray, int count);

static void 
sumPoints(uint32_t dim, double *p1, double *p2);

static void
normalize(uint32_t dim, double *p, uint32_t factor);

static int 
closestCluster(uint32_t dim, double *point,
               double *centers, uint32_t clusterCount);

static void 
calcNewCenters(double *newCenters, uint32_t *count,
               double *datapoints, uint32_t dim,
               uint32_t datasetSize, uint32_t clusterCount,
               double *centers);

static double *
fit(double *centroids, double *datapoints,
    uint32_t datasetSize,uint32_t dim,
    uint32_t clusterCount, uint32_t MAX_ITER);

void 
printMatrix(double* matrix, uint32_t rows, uint32_t cols);