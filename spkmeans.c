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

double* WAM(double *points, unsigned int obs_count, unsigned int dim) {
    double* WAM = (double*) malloc(obs_count*obs_count*sizeof(double));
    unsigned int i=0,j=0;

    if (WAM == NULL) {
        assert("malloc did an oopsie");
    }

    for (i = 0; i < obs_count; i++) {
        for (j = 0; j < i; j++) {
            WAM[i*obs_count + j] = exp(-0.5 * l2norm(points+i*dim, points+j*dim, dim));    
        }
        WAM[i*obs_count + i] = 0;
    }
    
    return WAM;
}

