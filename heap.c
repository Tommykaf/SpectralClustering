#include "heap.h"


uint32_t leftChild(uint32_t i){
    return 1 + (i << 1);
}

uint32_t rightChild(uint32_t i){
    return (1 + i) << 1;
}

void swap(element* x, element *y){
    element tmp = *x;
    *x = *y;
    *y = tmp;
}

void heapify(element *arr, uint32_t arrSize, uint32_t index){
    uint32_t minIndex = index;
    uint32_t rIndex = rightChild(index);
    uint32_t lIndex = leftChild(index);
    
    if (lIndex < arrSize && arr[lIndex].value < arr[minIndex].value){
        minIndex = lIndex;
    }

    if (rIndex < arrSize && arr[rIndex].value < arr[minIndex].value){
        minIndex = rIndex;
    }

    if (minIndex != index){
        swap(&arr[index], &arr[minIndex]);
        heapify(arr, arrSize, minIndex);
    }

}

element extractMin(element* heap, uint32_t len){
    element min = heap[0];
    swap(&heap[0], &heap[len - 1]);
    heapify(heap, len - 1, 0);
    return min;
}

void heapSort(double* vals, uint32_t len, uint32_t *ret, uint32_t k){
    element *heap = (element*) calloc(len, sizeof(element));
    uint32_t i;

    assert(heap != NULL);

    for(i = 0; i < len; i++){
        heap[i].index = i;
        heap[i].value = vals[i];
    }

    for(i = (len >> 1); i > 0; i--){
        heapify(heap, len, i - 1);
    }    
    
    for(i = 0; i < k; i++){
        ret[i] = extractMin(heap, len--).index;
    }

    free(heap);
}
