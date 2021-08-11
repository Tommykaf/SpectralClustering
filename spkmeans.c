// בס"ד

#define sqr(x) ((x)*(x))

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static double l2norm(unsigned int dim, double *p1, double *p2) {
  double dist = 0;
  unsigned int i;
  
  for (i = 0; i < dim; i++){
    dist += sqr(p1[i] - p2[i]);
  }

  dist = sqrt(dist);
  return dist;
}

double* WAM(double *points, unsigned int obsCount, unsigned int dim) {
    double* WAM = (double*) malloc(obsCount*obsCount*sizeof(double));
    unsigned int i,j;

    if (WAM == NULL) {
        assert("malloc did an oopsie");
    }

    for (i = 0; i < obsCount; i++) {
        for (j = 0; j < i; j++) {
            WAM[i*obsCount + j] = exp(-0.5*l2norm(points + i*dim, points + j*dim, dim));    
        }
        WAM[i*obsCount + i] = 0;
    }
    
    return WAM;
}

double sumRow(double *matrix, unsigned int cols, unsigned int rowIndex){
    double sum = 0;
    unsigned int i;
    for (i = 0; i < cols; i++) {
        sum += matrix[rowIndex*cols + i];
    }
    return sum;
    
}

double* DDG(double *WAM, unsigned int obsCount) {
    // Diagonal matrix, so let's just set all memory to 0s
    double* DDG = (double*) calloc(obsCount*obsCount, sizeof(double));
    unsigned int i;

    if (DDG == NULL) {
        assert("malloc did an oopsie");
    }

    for (i = 0; i < obsCount; i++) {
        DDG[i*obsCount + i] = sumRow(WAM, obsCount, i);
    }
    
    return DDG;
}

void DHalf(double *DDG, unsigned int obsCount) {
    unsigned int i;
    for (i = 0; i < obsCount; i++) {
        DDG[i*obsCount + i] = 1/sqrt(DDG[i*obsCount + i]);
    }
}
