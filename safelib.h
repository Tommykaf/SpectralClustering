#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#ifndef SAFELIB
#define SAFELIB

void *safeCalloc(size_t count, size_t size);
void *safeMalloc(size_t size);

#endif