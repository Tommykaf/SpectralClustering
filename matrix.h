#include <stdint.h>

#ifndef MATRIX_H
#define MATRIX_H

typedef struct _matrix {
  double* values;
  uint32_t rows;
  uint32_t cols;
} matrix;

#endif