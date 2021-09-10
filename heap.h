#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "safelib.h"

#ifndef HEAP_H
#define HEAP_H

typedef struct _heapElement
{
  uint32_t index; /* in original list */
  double value;   /* the value to sort by */
} element;

uint32_t leftChild(uint32_t);
uint32_t rightChild(uint32_t);

void swap(element *, element *);

void heapify(element *, uint32_t, uint32_t);

element 
extractMin(element *heap, uint32_t len);

void heapSort(double *vals, uint32_t len, uint32_t *ret, uint32_t k);

#endif