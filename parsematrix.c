#include "parsematrix.h"

void addPointToDataset(matrix_t *dataset, char line[])
{
  char *endptr;
  uint32_t datasetMaxLen = DATASET_INIT_LEN;
  uint32_t i = 0;

  if (dataset->rows == datasetMaxLen - 1 && DATASET_MAX_LEN - 1 > dataset->rows)
  {
    datasetMaxLen = MIN(2 * datasetMaxLen, DATASET_MAX_LEN);
    dataset->values = (double *)realloc(dataset->values, datasetMaxLen * sizeof(double));
    assert(dataset->values != NULL);
    {
      
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

void shrinkDataset(matrix_t *dataset)
{
  dataset->values = (double *)realloc(dataset->values, dataset->rows * dataset->cols * sizeof(double));
  assert(dataset->values != NULL);
}

/*
  Assumes that *dataset is allocated (not necessarily initialized) but (*dataset).values isn't allocated
  Also, in_file should be an open file on read mode
*/
void parseFile(char *in_file, matrix_t *dataset)
{
  FILE *fp = fopen(in_file, "r");
  uint32_t firstline = 1;
  char line[LINE_MAX_LEN];
  uint32_t i = 0;

  dataset->cols = 1;
  dataset->rows = 0;

  while (EOF != fscanf(fp, "%s\n", line) && (*line) != EOF)
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
      assert(dataset->values != NULL);
    }
    addPointToDataset(dataset, line);
  }

  shrinkDataset(dataset);
  fclose(fp);
}