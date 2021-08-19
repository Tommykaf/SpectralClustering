/* בס"ד */

#include "spkmeans.h"

static double l2norm(uint32_t dim, double *p1, double *p2)
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

static void WAM(double *points, uint32_t obsCount, uint32_t dim, double *res)
{
  uint32_t i, j;
  double tmp;

  for (i = 0; i < obsCount; i++)
  {
    for (j = 0; j < i; j++)
    {
      tmp = exp(-0.5 * l2norm(dim, points + i * dim, points + j * dim));
      res[i * obsCount + j] = tmp;
      res[j * obsCount + i] = tmp;
    }
    res[i * obsCount + i] = 0;
  }
}

double sumRow(double *matrix, uint32_t cols, uint32_t rowIndex)
{
  double sum = 0;
  uint32_t i;
  for (i = 0; i < cols; i++)
  {
    sum += matrix[rowIndex * cols + i];
  }
  return sum;
}

static void DDG(double *WAM, uint32_t obsCount, double *res)
{
  uint32_t i;

  for (i = 0; i < obsCount; i++)
  {
    res[i] = sumRow(WAM, obsCount, i);
  }
}

/* double invsqrtQuake( double number )
   {
       double y = number;
       double x2 = y * 0.5;
       int64_t i = *(int64_t *) &y;
       //The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
       i = 0x5fe6eb50c7b537a9 - (i >> 1);
       y = *(double *) &i;
       y = y * (1.5 - (x2 * y * y));   // 1st iteration
       //      y  = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
       return y;
   }
*/

void DHalf(double *DDG, uint32_t obsCount)
{
  uint32_t i;
  for (i = 0; i < obsCount; i++)
  {
    /*DDG[i] = 1 / invsqrtQuake(DDG[i]);*/
    DDG[i] = 1 / sqrt(DDG[i]);
  }
}

static double *prettyDDG(double *DDG, uint32_t obsCount)
{
  double *pretty = calloc(obsCount * obsCount, sizeof(double));
  uint32_t i;

  for (i = 0; i < obsCount; i++)
  {
    pretty[i * obsCount + i] = DDG[i];
  }
  return pretty;
}

static void laplacian(double *WAM, double *DHalf, uint32_t obsCount, double *LNorm)
{
  uint32_t i, j;
  double tmp;

  for (i = 0; i < obsCount; i++)
  {
    for (j = 0; j < i; j++)
    {
      tmp = -(DHalf[i] * WAM[i * obsCount + j] * DHalf[j]);
      LNorm[i * obsCount + j] = tmp;
      LNorm[j * obsCount + i] = tmp;
    }
    LNorm[i * obsCount + i] = 1.0;
  }
}


void maxItem(double *matrix, uint32_t rows, uint32_t *row, uint32_t *col)
{
  /* This function places the bigget absolute value item i at indices[0] and his j at indices[1]
  TODO make it not O(n^2) */
  uint32_t i, j, currI = 1, currJ = 0;
  double tmp, currMax = fabs(matrix[currJ + currI * rows]);
  for (i = 1; i < rows; i++)
  {
    for (j = 0; j < i; j++)
    {
      if ((tmp = fabs(matrix[j + i * rows])) > currMax)
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
static void Jacobi(double *matrix, uint32_t rows, double *V)
{
  uint32_t i, j, r, count;
  double theta, s, t, c, tmp;
  double arj, ari, aii, ajj, aij;
  double diff;

  for (i = 0; i < rows; i++)
  {
    V[i * rows + i] = 1.0;
  }

  while (diff < EPSILON && count++ < MAX_JACOBI_ITER)
  {
    diff = 0;
    maxItem(matrix, rows, &i, &j);
    aij = matrix[i * rows + j];
    aii = matrix[i * (rows + 1)];
    ajj = matrix[j * (rows + 1)];
    theta = (ajj - aii) / (2 * aij);
    t = (theta >= 0 ? 1 : -1) / (fabs(theta) + sqrt(sqr(theta) + 1));
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

    diff -= (sqr(aii) + sqr(ajj));
    diff += sqr(matrix[i * (rows + 1)]) + sqr(matrix[j * (rows + 1)]);
  }
}

uint32_t argmax(double *eigenArray, uint32_t count)
{
  uint32_t i, k;
  double delta, tmp;
  for (i = 0; i < count / 2; i++)
  {
    if ((tmp = fabs(eigenArray[i] - eigenArray[i + 1])) > delta)
    {
      delta = tmp;
      k = i;
    }
  }
  return k;
}

/* build Eigenvector matrix, ret can't be initialized */
void buildU(double *eigenArray, uint32_t count, double* V, double* ret){
  uint32_t k = argmax(eigenArray, count);
  uint32_t i, j;
  uint32_t *indices = malloc(sizeof(uint32_t) * k);
  ret = malloc(sizeof(double) * k * count); /* n by k */
  
  assert(indices != NULL);
  assert(ret != NULL);
  
  heapSort(eigenArray, count, indices ,k);
  /* indices now contains the k first eigenValues indices */

  for(i = 0; i < k; i++){
    for(j = 0; j < count; j++){
      ret[j * k + i] = V[indices[i] * count + j];
      /* V was calculted transposed => V[i,j] is now V[j,i]
      equivlent to ret[j*k + i] = V[j*count + indices[i]] if V wasn't rotated
      should save us some cache misses $P */
    }
  }
}


static void sumPoints(uint32_t dim, double *p1, double *p2)
{
  uint32_t i;
  for (i = 0; i < dim; i++)
  {
    p1[i] += p2[i];
  }
}

static void normalize(uint32_t dim, double *p, uint32_t factor)
{
  uint32_t i;
  for (i = 0; i < dim; i++)
  {
    p[i] /= factor;
  }
}

static uint32_t closestCluster(uint32_t dim, double *point, double *centers, uint32_t clusterCount)
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

static void calcNewCenters(double *newCenters, uint32_t *count, double *datapoints, uint32_t dim,
               uint32_t datasetSize, uint32_t clusterCount, double *centers)
{
  uint32_t cluster;
  uint32_t i;
  uint32_t j;

  memset(newCenters, 0, sizeof(double) * dim * clusterCount);
  memset(count, 0, sizeof(int) * clusterCount);

  for (j = 0; j < datasetSize; j++)
  {
    cluster = closestCluster(dim, &(datapoints[j * dim]), centers, clusterCount);
    sumPoints(dim, &(newCenters[cluster * dim]), &(datapoints[j * dim]));
    count[cluster]++;
  }

  for (i = 0; i < clusterCount; i++)
  {
    normalize(dim, &(newCenters[i * dim]), count[i]);
  }
}

static double *fit(double *centroids, double *datapoints, uint32_t datasetSize,
    uint32_t dim, uint32_t clusterCount, uint32_t MAX_ITER)
{
  uint32_t i;

  double *newCentroids = calloc(clusterCount, dim * sizeof(double));
  uint32_t *count = calloc(clusterCount, sizeof(int));

  for (i = 0; i < MAX_ITER; i++)
  {
    calcNewCenters(newCentroids, count, datapoints, dim, datasetSize, clusterCount, centroids);
    if (memcmp(centroids, newCentroids, clusterCount * dim * sizeof(double)))
    {
      break;
    }

    memcpy(centroids, newCentroids, sizeof(double) * dim * clusterCount);
  }

  free(newCentroids);
  free(count);

  return centroids;
}

static void printMatrix(double *matrix, uint32_t rows, uint32_t cols)
{
  uint32_t i, j;
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      printf("%.4f", matrix[i * cols + j]);
      if (j < cols - 1)
        printf(",");
    }
    if (i + 1 < rows) 
    /* if not last row - tommy 2021 */
      printf("\n");
  }
}