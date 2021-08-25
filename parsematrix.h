#define DATASET_INIT_LEN (20)
#define DATASET_MAX_LEN (50)
#define LINE_MAX_LEN (200) /* TODO: Change according to instructions */

#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "matrix.h"

void addPointToDataset(matrix* dataset, char line[], unsigned int * datasetMaxLen);
void shrinkDataset(matrix* dataset);
void parseFile(FILE *in_file, matrix* dataset);