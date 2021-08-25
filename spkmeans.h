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
#define EPSILON 0.001
#define MAX_JACOBI_ITER 100


static double 
l2norm(uint32_t dim, double *p1, double *p2);

static void
WAM(double *points, uint32_t obsCount, uint32_t dim, double *res);

double 
sumRow(double *matrix, uint32_t cols, uint32_t rowIndex);

static void 
DDG(double *WAM, uint32_t obsCount, double *res);

/* double 
    invsqrtQuake( double number );
*/
void 
DHalf(double *DDG, uint32_t obsCount);

static double* 
prettyDDG(double* DDG, uint32_t obsCount);

static void /* HEZY */
laplacian(double* WAM, double* DHalf, uint32_t obsCount, double *LNorm);

void 
maxItem(double *matrix, uint32_t rows, uint32_t *row, uint32_t *col);

static void 
Jacobi(double *matrix, uint32_t rows, double *V, double* eigenArray);

static uint32_t 
argmax(double *eigenArray, uint32_t count);

static void 
buildU(double *eigenArray, uint32_t count,
            uint32_t k, double* V, matrix* ret);

static matrix*
prepareData(double *points, uint32_t obsCount, uint32_t dim);

static void 
sumPoints(uint32_t dim, double *p1, double *p2);

static void
normalize(uint32_t dim, double *p, uint32_t factor);

static uint32_t 
closestCluster(uint32_t dim, double *point, 
               double *centers, uint32_t clusterCount);

static void 
calcNewCenters(double *newCenters, uint32_t *count,
               double *datapoints, uint32_t dim,
               uint32_t datasetSize, uint32_t clusterCount,
               double *centers);

static double *
kmeansFit(double *centroids, double *datapoints,
    uint32_t datasetSize,uint32_t dim,
    uint32_t clusterCount, uint32_t MAX_ITER);

void 
printMatrix(double* matrix, uint32_t rows, uint32_t cols);