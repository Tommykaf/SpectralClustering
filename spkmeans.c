/* בס"ד */

#include "spkmeans.h"

double l2norm(uint32_t dim, double *p1, double *p2)
{
  double dist = 0;
  uint32_t i;

  for (i = 0; i < dim; i++)
  {
    dist += sqr(p1[i] - p2[i]);
  }

  dist = sqrt(dist);
  return dist;
}

double vectorLength(uint32_t dim, double *v)
{
  double dist = 0;
  uint32_t i;

  for (i = 0; i < dim; i++)
  {
    dist += sqr(v[i]);
  }

  dist = sqrt(dist);
  return dist;
}

/*
Calculates the WAM of the graph of the given <points>, places it in res
args:
  matrix - observation matrix
  res - WAM of the original point matrix - size rows x rows
*/
void WAM(matrix_t *matrix, matrix_t *res)
{
  uint32_t i, j;
  double tmp;
  uint32_t obsCount = matrix->rows;
  uint32_t dim = matrix->cols;
  double *points = matrix->values;
  
  res->cols = obsCount;
  res->rows = obsCount;

  res->values = (double *) safeMalloc(res->rows * res->cols * sizeof(double));

  for (i = 0; i < obsCount; i++)
  {
    for (j = 0; j < i; j++)
    {
      tmp = l2norm(dim, &points[i * dim], &points[j * dim]);
      tmp = exp(-0.5 * l2norm(dim, &points[i * dim], &points[j * dim]));
      res->values[i * obsCount + j] = tmp;
      res->values[j * obsCount + i] = tmp;
    }
    res->values[i * obsCount + i] = 0;
  }
}

/*
Sum <cols> elements in a given <rowIndex> from <matrix> and returns the sum
args:
  matrix - the matrix to sum from
  rowIndex - the row index to sum
*/
double sumRow(matrix_t *matrix, uint32_t rowIndex)
{
  double sum = 0;
  uint32_t i;
  for (i = 0; i < matrix->cols; i++)
  {
    sum += matrix->values[rowIndex * matrix->cols + i];
  }
  return sum;
}

/*
Calculates the diagonal of the DDG of the given WAM, places it in res
args:
  WAM - a WAM matrix
  res - double array representing diagonal of the DDG
*/
void DDG(matrix_t *WAM, double *res)
{
  uint32_t i;

  for (i = 0; i < WAM->rows; i++)
  {
    res[i] = sumRow(WAM, i);
  }
}

/*
Modifies DDG to be DHalf
args:
  obsCount - number of observations
  DDG - double array of obsCount length
*/
void DHalf(double *DDG, uint32_t obsCount)
{
  uint32_t i;
  for (i = 0; i < obsCount; i++)
  {
    DDG[i] = 1 / sqrt(DDG[i]);
  }
}

/*
creates a nice matrix represantion of a DDG/DHALF
args:
  DDG - the diagonal values array
  obsCount - number of values in the array
*/
void prettyDDG(double *DDG, uint32_t obsCount)
{
  matrix_t *pretty = (matrix_t *) safeMalloc(sizeof(matrix_t));

  uint32_t i;

  pretty->cols = obsCount;
  pretty->rows = obsCount;

  pretty->values = (double *) safeCalloc(obsCount * obsCount, sizeof(double));

  for (i = 0; i < pretty->rows; i++)
  {
    pretty->values[i * pretty->cols + i] = DDG[i];
  }
  
  printMatrix(pretty->values, pretty->rows, pretty->cols);

  free(pretty->values);
  free(pretty);
}


/*
calculates the laplacian from WAM and DDG/DHALF, places it into LNorm
args:
  DHalf - the Dhalf diagonal values array
  WAM - the WAM matrix
  obsCount - number of values in the array
  LNorm - The returned laplacian matrix
*/
void laplacian(matrix_t *WAM, double *DHalf, matrix_t *LNorm)
{
  uint32_t i, j;
  double tmp;
  LNorm->rows = WAM->rows;
  LNorm->cols = WAM->cols;

  LNorm->values = (double *) safeMalloc(sqr(LNorm->rows) * sizeof(double));

  for (i = 0; i < LNorm->rows; i++)
  {
    for (j = 0; j < i; j++)
    {
      tmp = -(DHalf[i] * WAM->values[i * WAM->cols + j] * DHalf[j]);
      LNorm->values[i * LNorm->cols + j] = tmp;
      LNorm->values[j * LNorm->cols + i] = tmp;
    }
    LNorm->values[i * LNorm->cols + i] = 1.0;
  }
}


/* 
Find the max item in the upper triangle Matrix and return his indices
args:
  matrix - the square matrix of dim rows to search in
  row - the max item row index - return value 0
  col - the max item col index - return value 2
*/
void maxItem(matrix_t *matrix, uint32_t *row, uint32_t *col)
{
  uint32_t i, j, currI = 1, currJ = 0;
  double tmp, currMax = fabs(matrix->values[currJ + currI * matrix->rows]);
  for (i = 1; i < matrix->rows; i++)
  {
    for (j = 0; j < i; j++)
    {
      if ((tmp = fabs(matrix->values[j + i * matrix->rows])) > currMax)
      {
        currI = i;
        currJ = j;
        currMax = tmp;
      }
    }
  }
  *row = currI;
  *col = currJ;
}

/* Assumes matrix is symetric, and V is full of zeros */
void Jacobi(matrix_t *input_matrix, double *V, double *eigenArray)
{
  uint32_t i, j, r, count = 0;
  double theta, s, t, c, tmp;
  double arj, ari, aii, ajj, aij;
  double diff = 2 * EPSILON;

  uint32_t rows = input_matrix->rows;
  double *matrix = input_matrix->values;

  for (i = 0; i < rows; i++)
  {
    V[i * rows + i] = 1.0;
  }

  while (diff > EPSILON && count++ < MAX_JACOBI_ITER)
  {
    diff = 0;
    maxItem(input_matrix, &j, &i);
    aij = matrix[i * rows + j]; /* The Max off diag */
    aii = matrix[i * (rows + 1)];
    ajj = matrix[j * (rows + 1)];
    theta = (ajj - aii) / (2 * aij);
    t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(sqr(theta) + 1));
    c = 1 / sqrt(sqr(t) + 1);
    s = t * c;

    for (r = 0; r < rows; r++)
    {
      if (r != j && r != i)
      {
        arj = matrix[r * rows + j];
        ari = matrix[r * rows + i];
        matrix[r * rows + i] = c * ari - s * arj;
        matrix[i * rows + r] = c * ari - s * arj;
        matrix[r * rows + j] = c * arj + s * ari;
        matrix[j * rows + r] = c * arj + s * ari;
        diff += 2 * (sqr(arj) + sqr(ari));
        diff -= 2 * (sqr(matrix[r * rows + i]) + sqr(matrix[r * rows + j]));
      }
      /* REGULAR
      tmp = V[r * rows + i];
      V[r * rows + i] = tmp - s * (V[r * rows + j] + s * tmp / (1 + c));
      V[r * rows + j] = V[r * rows + j] + s * (tmp - s * V[r * rows + j] / (1 + c));
      "Transpose"  - calculating on the side
      */
      tmp = V[i * rows + r];
      V[i * rows + r] = tmp - s * (V[j * rows + r] + s * tmp / (1 + c));
      V[j * rows + r] = V[j * rows + r] + s * (tmp - s * V[j * rows + r] / (1 + c));
    }
    matrix[i * (rows + 1)] = sqr(c) * aii + sqr(s) * ajj - 2 * c * s * aij;
    matrix[j * (rows + 1)] = sqr(s) * aii + sqr(c) * ajj + 2 * c * s * aij;
    matrix[i * rows + j] = 0;
    matrix[j * rows + i] = 0;

    diff += 2 * sqr(aij);
  }
  for (r = 0; r < rows; r++)
  {
    eigenArray[r] = matrix[r * (rows + 1)];
  }
}

/* Eigengap hurestic*/
uint32_t argmax(double *eigenArray, uint32_t count, uint32_t *indices)
{
  uint32_t i, k = 0;
  double tmp, delta = -1.0;

  for (i = 0; i < count / 2; i++)
  {
    if ((tmp = fabs(eigenArray[indices[i]] - eigenArray[indices[i + 1]])) > delta)
    {
      delta = tmp;
      k = i+1;
    }
  }
  return k;
}

/* build Eigenvector matrix, ret must be initialized */
void buildT(double *eigenArray, uint32_t count, uint32_t k, double* V, matrix_t *ret){
  uint32_t i, j;
  uint32_t *indices = (uint32_t *) safeMalloc(sizeof(uint32_t) * k);
  double len;


  heapSort(eigenArray, count, indices, k);
  /* indices now contains the k first eigenValues indices */
  ret->rows = count;
  ret->cols = k;

  ret->values = (double *) safeMalloc(sizeof(double) * k * count); /* n by k */

  for (i = 0; i < k; i++)
  {
    for (j = 0; j < count; j++)
    {
      ret->values[j * k + i] = V[indices[i] * count + j];
      /* V was calculted transposed => V[i,j] is now V[j,i]
      equivlent to ret[j*k + i] = V[j*count + indices[i]] if V wasn't rotated
      should save us some cache misses $P */
    }
  }
  for (i = 0; i < count; i++)
  {
    len = vectorLength(ret->cols, &ret->values[i*ret->cols]);
    normalize(ret->cols, &(ret->values[i*ret->cols]), len);
  }
  free(indices);
}

matrix_t *prepareData(matrix_t *points, uint32_t k){
  uint32_t obsCount = points->rows;
  matrix_t *wam = (matrix_t *) safeMalloc(sizeof(matrix_t));
  matrix_t *LNorm = (matrix_t *) safeMalloc(sizeof(matrix_t));
  matrix_t *DATA = (matrix_t *) safeMalloc(sizeof(matrix_t));
  double *D = (double *) safeMalloc(obsCount * sizeof(double));
  double *V = (double *) safeCalloc(sqr(obsCount), sizeof(double));
  double *eigenArray = (double *) safeMalloc(obsCount * sizeof(double));
  uint32_t *indices = (uint32_t*) safeCalloc(obsCount, sizeof(uint32_t));
  
  /* wam */
  WAM(points, wam);
  /* ddg + dhalf */
  DDG(wam, D);
  DHalf(D, obsCount);
  /*lnorm */
  laplacian(wam, D, LNorm);
  /* jacob */
  Jacobi(LNorm, V, eigenArray);
  /* build T */
  heapSort(eigenArray, obsCount, indices, obsCount);
  k = k != 0 ? k : argmax(eigenArray, obsCount, indices); /* if k is 0 we do the hueristic */
  buildT(eigenArray, obsCount, k, V, DATA);

  free(wam->values);
  free(wam);
  free(LNorm->values);
  free(LNorm);

  free(D);
  free(eigenArray);
  free(indices);
  free(V);

  return DATA;
}

void sumPoints(uint32_t dim, double *p1, double *p2)
{
  uint32_t i;
  for (i = 0; i < dim; i++)
  {
    p1[i] += p2[i];
  }
}

void normalize(uint32_t dim, double *p, double factor)
{
  uint32_t i;
  for (i = 0; i < dim; i++)
  {
    p[i] /= factor;
  }
}

uint32_t closestCluster(uint32_t dim, double *point, double *centers, uint32_t clusterCount)
{
  double min = 0, dist = 0;
  uint32_t minIndex = 0, firstRun = 1;
  uint32_t i;
  for (i = 0; i < clusterCount; i++)
  {
    dist = l2norm(dim, point, &(centers[i * dim]));
    if ((dist < min) || firstRun)
    {
      min = dist;
      minIndex = i;
      firstRun = 0;
    }
  }
  return minIndex;
}

void calcNewCenters(double *newCenters, uint32_t *count, double *datapoints, uint32_t dim,
                           uint32_t datasetSize, uint32_t clusterCount, double *centers)
{
  uint32_t cluster;
  uint32_t i;
  uint32_t j;

  memset(newCenters, 0, sizeof(double) * dim * clusterCount);
  memset(count, 0, sizeof(uint32_t) * clusterCount);

  for (j = 0; j < datasetSize; j++)
  {
    cluster = closestCluster(dim, &(datapoints[j * dim]), centers, clusterCount);
    sumPoints(dim, &(newCenters[cluster * dim]), &(datapoints[j * dim]));
    count[cluster]++;
  }

  for (i = 0; i < clusterCount; i++)
  {
    normalize(dim, &(newCenters[i * dim]), (double) count[i]);
  }
}

void kmeansFit(double *centroids, double *datapoints, uint32_t datasetSize,
                         uint32_t dim, uint32_t clusterCount)
{
  uint32_t i;

  double *newCentroids = (double *) safeCalloc(clusterCount, dim * sizeof(double));
  uint32_t *count = (uint32_t *) safeCalloc(clusterCount, sizeof(int));
  

  for (i = 0; i < MAX_KMEANS_ITER; i++)
  {
    calcNewCenters(newCentroids, count, datapoints, dim, datasetSize, clusterCount, centroids);
    if (memcmp(centroids, newCentroids, clusterCount * dim * sizeof(double)) == 0)
    {
      break;
    }

    memcpy(centroids, newCentroids, sizeof(double) * dim * clusterCount);
  }

  free(newCentroids);
  free(count);
}

/* Handles the -0 situation */
double zerod(double num)
{
  return -0.00005 < num && num < 0 ? 0 : num;
}

void printMatrix(double* values, uint32_t rows, uint32_t cols)
{
  uint32_t i, j;
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      printf("%.4f", zerod(values[i * cols + j]));
      if (j < cols - 1)
        printf(",");
    }
    if (i + 1 < rows)
      /* if not last row - tommy 2021 */
      printf("\n");
  }
}

int main(int argc, char* argv[]) {
  uint32_t K;
  uint32_t *indices;
  matrix_t *initial_input, *wam, *lNorm, *data;
  double *eigenArray, *V, *D, *centroids;

  assert(argc == 4);
  K = atoi(argv[1]);
  initial_input = (matrix_t *) safeCalloc(1, sizeof(matrix_t));

  assert(strcmp(argv[2], "jacobi") == 0 || 
        strcmp(argv[2], "wam") == 0||
        strcmp(argv[2], "ddg") == 0||
        strcmp(argv[2], "lnorm") == 0||
        strcmp(argv[2], "spk") == 0 
  );

  parseFile(argv[3], initial_input);

  if (strcmp(argv[2], "jacobi") == 0) {
    V = (double *) safeCalloc(sqr(initial_input->rows), sizeof(double));
    eigenArray = (double *) safeMalloc(initial_input->rows * sizeof(double));
    Jacobi(initial_input, V, eigenArray);

    printMatrix(eigenArray, 1, initial_input->rows);
    printf("\n"); /* printMatrix doesnt print another \n at the end*/
    printMatrix(V, initial_input->rows, initial_input->rows);
    
    free(V);
    free(eigenArray);

  } else {
    wam = (matrix_t *) safeCalloc(1, sizeof(matrix_t));
    WAM(initial_input, wam);
    
    if (strcmp(argv[2], "wam") == 0) {  
      printMatrix(wam->values, wam->rows, wam->cols);
    } else {
      D = (double *) safeMalloc(wam->rows * sizeof(double));
      DDG(wam, D);

      if (strcmp(argv[2], "ddg") == 0) {
        prettyDDG(D, wam->rows);
      } else {
        lNorm = (matrix_t *) safeCalloc(1, sizeof(matrix_t));
        
        DHalf(D, wam->rows);
        laplacian(wam, D, lNorm);
        if (strcmp(argv[2], "lnorm") == 0) {
          printMatrix(lNorm->values, lNorm->rows, lNorm->cols);
        } 
        else {
          if (strcmp(argv[2], "spk") == 0) {
            V = (double *) safeCalloc(sqr(initial_input->rows), sizeof(double));
            eigenArray = (double*) safeMalloc(initial_input->rows * sizeof(double));
            indices = (uint32_t*) safeCalloc(initial_input->rows, sizeof(uint32_t));

            Jacobi(lNorm, V, eigenArray);
            heapSort(eigenArray, lNorm->rows, indices, lNorm->rows);
            K = K != 0 ? K : argmax(eigenArray, lNorm->rows, indices);
            
            data = (matrix_t *) safeCalloc(1, sizeof(matrix_t));
            buildT(eigenArray, lNorm->rows, K, V, data);

          
            centroids = (double *) safeMalloc(sqr(data->cols) * sizeof(double));
            memcpy(centroids, data->values, sqr(data->cols) * sizeof(double));

            kmeansFit(centroids, data->values, data->rows, data->cols, data->cols);
            printMatrix(centroids, data->cols, data->cols);

            free(V);
            free(eigenArray);
            free(indices);
            free(data->values);
            free(data);
            free(centroids);
            }
          }
        free(lNorm->values);
        free(lNorm);
        }
      free(D);
      }
    free(wam->values);
    free(wam);
    }
    free(initial_input->values);
    free(initial_input);
    return 0;
}