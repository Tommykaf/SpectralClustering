/* אם ירצה השם */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


#include "heap.h"
#include "parsematrix.h"

#include "matrix.h"


#define sqr(x) ((x)*(x))
#define EPSILON 1e-15
#define MAX_JACOBI_ITER 100
#define MAX_KMEANS_ITER 300


double 
l2norm(uint32_t dim, double *p1, double *p2);

void
WAM(matrix_t *matrix, matrix_t *res);

double 
sumRow(matrix_t *matrix, uint32_t rowIndex);

void 
DDG(matrix_t *WAM, double *res);

void 
DHalf(double *DDG, uint32_t obsCount);

void 
prettyDDG(double* DDG, uint32_t obsCount);

void /* HEZY */
laplacian(matrix_t* WAM, double* DHalf, matrix_t *LNorm);

void 
maxItem(matrix_t *matrix, uint32_t *row, uint32_t *col);

void 
Jacobi(matrix_t *input_matrix, double *V, double* eigenArray);

uint32_t 
argmax(double *eigenArray, uint32_t count);

void 
buildT(double *eigenArray, uint32_t count,
            uint32_t k, double* V, matrix_t* ret);

matrix_t*
prepareData(matrix_t *points, uint32_t k);

void 
sumPoints(uint32_t dim, double *p1, double *p2);

void
normalize(uint32_t dim, double *p, double factor);

uint32_t 
closestCluster(uint32_t dim, double *point, 
               double *centers, uint32_t clusterCount);

void 
calcNewCenters(double *newCenters, uint32_t *count,
               double *datapoints, uint32_t dim,
               uint32_t datasetSize, uint32_t clusterCount,
               double *centers);

double *
kmeansFit(double *centroids, double *datapoints,
    uint32_t datasetSize,uint32_t dim,
    uint32_t clusterCount);

void 
printMatrix(matrix_t* matrix);