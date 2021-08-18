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
maxItem(double* matrix, int rows, int* row, int* col) {
    // This function places the bigget absolute value item i at indices[0] and his j at indices[1]
    // TODO make it not O(n^2)
    unsigned int i, j, currI = 1, currJ = 0;
    double tmp, currMax = abs(matrix[currJ + currI*rows]);
    for(i = 1; i < rows; i++){
        for(j = 0; j < i; j++){
            if ((tmp = abs(matrix[j + i*rows])) > currMax) {
                currI = i;
                currJ = j;
                currMax = tmp;
            }
        }
    }
    *row = currI;
    *col = currJ;
}

void
Jacobi(double* matrix, unsigned int rows, double old){
    double* V = (double*) calloc(rows * rows, sizeof(double));
    unsigned int i, j, r, count;
    int indices[2];
    double theta, s, t, c, tmp;
    double arj, ari, aii, ajj, aij;
    double diff;

    for (i = 0; i < rows; i++){
        V[i * rows + i] = 1.0;
    }

    while(diff < EPSILON && count++ < 100){
        diff = 0;
        maxItem(matrix, rows, &i, &j);
        aij = matrix[i*rows + j];
        aii = matrix[i*(rows+1)];
        ajj = matrix[j*(rows+1)];
        theta = (ajj - aii) / (2 * aij);
        t = (theta >= 0 ? 1 : -1) / (abs(theta) + sqrt(sqr(theta) + 1));
        c = 1 / sqrt(sqr(t) + 1);
        s = t * c;

        for(r = 0; r < rows; r++){
            if (r != j && r != i) {
                arj = matrix[r*rows + j];
                ari = matrix[r*rows + i];
                matrix[r*rows+i] = c*ari - s*arj;
                matrix[i*rows+r] = c*ari - s*arj;
                matrix[r*rows+j] = c*arj + s*ari;
                matrix[j*rows+r] = c*arj + s*ari;
                diff += 2 * (sqr(arj) + sqr(ari));
                diff -= 2 * (sqr(matrix[r*rows+i]) + sqr(matrix[r*rows+j]));
            }
            tmp = V[r*rows + i];
            V[r*rows + i] = tmp - s*(V[r*rows + j] + s*tmp/(1 + c));
            V[r*rows + j] = V[r*rows + j] + s*(tmp - s*V[r*rows + j]/(1 + c));
        }
        matrix[i*(rows + 1)] = sqr(c)*aii + sqr(s)*ajj - 2*c*s*aij;
        matrix[j*(rows + 1)] = sqr(s)*aii + sqr(c)*ajj + 2*c*s*aij;
        matrix[i*rows + j] = 0;

        diff -= (sqr(aii) + sqr(ajj)); 
        diff += sqr(matrix[i*(rows + 1)]) + sqr(matrix[j*(rows + 1)]);
    }
}

int
argmax(double* eigenArray,int count) {
    unsigned int i, k;
    double delta, tmp;
    for(i = 0; i < count / 2; i++){
        if ((tmp = abs(eigenArray[i] - eigenArray[i + 1])) > delta){
            delta = tmp;
            k = i;
        }
    }
    return k;
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