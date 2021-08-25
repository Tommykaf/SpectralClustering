#include "parsematrix.h"

void addPointToDataset(matrix *dataset, char line[], unsigned int *datasetMaxLen)
{
  char *endptr;
  unsigned int i = 0;

  if (dataset->rows == (*datasetMaxLen) - 1 && DATASET_MAX_LEN - 1 > dataset->rows)
  {
    (*datasetMaxLen) = MIN(2 * (*datasetMaxLen), DATASET_MAX_LEN);
    dataset->values = (double *)realloc(dataset->values, (*datasetMaxLen) * sizeof(double));
    if (dataset->values == NULL)
    {
      assert("Realloc failed :(");
    }
  }

  endptr = line;
  for (i = 0; i < (*dataset).cols; i++)
  {
    dataset->values[dataset->rows * dataset->cols + i] = strtod(endptr, &endptr);
    endptr++;
  }

  (dataset->rows)++;
}

void shrinkDataset(matrix *dataset)
{
  dataset->values = (double *)realloc(dataset->values, dataset->rows * dataset->cols * sizeof(double));
  if ((*dataset).values == NULL)
  {
    assert("Realloc failed :(");
  }
}

/*
  Assumes that *dataset is allocated (not necessarily initialized) but (*dataset).values isn't allocated
  Also, in_file should be an open file on read mode
*/
void parseFile(FILE *in_file, matrix *dataset)
{

  unsigned int datasetMaxLen = DATASET_INIT_LEN;
  int firstline = 1;
  char line[LINE_MAX_LEN];
  unsigned int i = 0;

  dataset->cols = 0;
  dataset->rows = 0;

  while (EOF != fscanf(in_file, "%s\n", line) && (*line) != EOF)
  {
    if (firstline)
    {
      firstline = 0;
      for (i = 0; i < sizeof(line); i++)
      {
        if (',' == line[i])
        {
          (*dataset).cols++;
        }
      }
      dataset->values = calloc(DATASET_INIT_LEN, dataset->cols * sizeof(double));
    }
    addPointToDataset(dataset, line, &datasetMaxLen);
  }

  shrinkDataset(dataset);
}