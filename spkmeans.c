// בס"ד

#include <spkmeans.h>

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
    double tmp;

    if (WAM == NULL) {
        assert("malloc did an oopsie");
    }

    for (i = 0; i < obsCount; i++) {
        for (j = 0; j < i; j++) {
            tmp = exp(-0.5*l2norm(points + i*dim, points + j*dim, dim));    
            WAM[i*obsCount + j] = tmp;
            WAM[j*obsCount + i] = tmp;
        }
        WAM[i*obsCount + i] = 0;
    }
    
    return WAM;
}

double sumRow(double *matrix, unsigned int cols, unsigned int rowIndex) {
    double sum = 0;
    unsigned int i;
    for (i = 0; i < cols; i++) {
        sum += matrix[rowIndex*cols + i];
    }
    return sum;
    
}

double* DDG(double *WAM, unsigned int obsCount) {
    double* DDG = (double*) malloc(obsCount * sizeof(double));
    unsigned int i;

    if (DDG == NULL) {
        assert("malloc did an oopsie");
    }

    for (i = 0; i < obsCount; i++) {
        DDG[i] = sumRow(WAM, obsCount, i);
    }
    
    return DDG;
}

void DHalf(double *DDG, unsigned int obsCount) {
    unsigned int i;
    for (i = 0; i < obsCount; i++) {
        DDG[i] = 1/sqrt(DDG[i]);
    }
}

double* prettyDDG(double* DDG, unsigned int obsCount) {
    double *pretty = calloc(obsCount*obsCount, sizeof(double));
    unsigned int i;

    for(i = 0; i < obsCount; i++)
    {
        pretty[i * obsCount + i] = DDG[i];
    }
    return pretty;
}

double* laplacian(double* WAM, double* DHalf, unsigned int obsCount){   
    double* LNorm = (double* ) malloc(obsCount*obsCount*sizeof(double));
    unsigned int i,j;
    double tmp;

    for(i = 0; i < obsCount; i++){
        for(j = 0; j < i; j++){
            tmp = -(DHalf[i] * WAM[i * obsCount + j] * DHalf[j]);
            LNorm[i * obsCount + j] = tmp; 
            LNorm[j * obsCount + i] = tmp;
        }
        LNorm[i * obsCount + i] = 1.0;
    }
    return LNorm;
}


void
multiplyMatrices(double* A, double* B, double* C, unsigned int m, unsigned int k, unsigned int n){
    // A = Mxk, B = KxN => C=MxN
    unsigned int i, j, l;
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            for (l = 0; l < k; l++){
                C[i * n + j] += A[i * k + l] * B[l * n + j];
            }
        }
    }
}

int 
isDiagonal(double* matrix, unsigned int rows){
    return 0;
}

void 
generateRot(double* matrix, unsigned int rows, double* p, double* pt){

}

void
Jacobi(double* matrix, unsigned int rows){
    double* V = (double*) calloc(rows * rows, sizeof(double));
    double *P = (double*) calloc(rows * rows, sizeof(double));
    double *Pt = (double*) calloc(rows * rows, sizeof(double));
    double *tmp= (double*) calloc(rows * rows, sizeof(double));
    unsigned int i;

    for (i = 0; i < rows; i++){
        V[i * rows + i] = 1.0;
    }

    while(!isDiagonal(matrix, rows)){
        generateRot(matrix, rows, P, Pt);
        multiplyMatrices(Pt, matrix,tmp, rows, rows, rows); // Pt * matrix => temp
        multiplyMatrices(tmp, P, matrix, rows, rows, rows); // pt* matrix * p => matrix
        multiplyMatrices(V, P, tmp, rows, rows, rows);      // p1 * p2 * ... = V * P_n+1 => tmp
        memcpy(V, tmp, rows * rows * sizeof(double));       // tmp => V
    }
    free(P);
    free(Pt);
}




void printMatrix(double* matrix, unsigned int rows, unsigned int cols){
    unsigned int i, j;
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            printf("%.4f", matrix[i*cols + j]);
            if (j < cols - 1)
                printf(",");
        }
        if (i + 1 < rows) // if not last row - tommy 2021 
            printf("\n");
    }
}