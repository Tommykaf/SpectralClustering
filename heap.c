#include <heap.h>




uint32_t leftChild(uint32_t i){
    return 1 + i << 1;
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
    uint32_t maxIndex = index;
    uint32_t rIndex = rightChild(index);
    uint32_t lIndex = leftChild(index);
    
    if (lIndex < arrSize && arr[lIndex].value > arr[maxIndex].value){
        maxIndex = lIndex;
    }

    if (rIndex < arrSize && arr[rIndex].value > arr[maxIndex].value){
        maxIndex = rIndex;
    }

    if (maxIndex != index){
        swap(&arr[index], &arr[maxIndex]);
    }

    heapify(arr, arrSize, maxIndex);
}


void heapSort(double* vals, uint32_t len, uint32_t *ret, uint32_t k){
    element *heap = (element*) calloc(len, sizeof(element));
    uint32_t i;
    double val;
    
    assert(heap != NULL);

    for(i = 0; i < len; i++){
        val = vals[i];
        heap[i].index = i;
        heap[i].value = val;
    }

    for(i = (len >> 1) - 1; i >= 0; i--){
        heapify(heap, len, i);
    }    
    
    for(i = len - 1; i >= 0; i--){
        swap(&heap[i], &heap[0]);
        heapify(heap, len, i);
    }

    for(i = 0; i < k; i++){
        ret[i] = heap[i].index;
    }

    free(heap);
}
