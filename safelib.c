#include "safelib.h"

void *safeCalloc(size_t count, size_t size)
{
    void *ptr = calloc(count, size);
    if (ptr == NULL){
        printf("An error has Occured");
        assert(ptr != NULL);
    }
    return ptr;
}

void *safeMalloc(size_t size){
    void *ptr = malloc(size);
    if (ptr == NULL){
        printf("An error has Occured");
        assert(ptr != NULL);
    }
    return ptr;
}

