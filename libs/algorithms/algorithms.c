#include "algorithms.h"

void inputArray(int *a, size_t size) {
    for (size_t i = 0; i < size; i++) {
        scanf("%d", &a[i]);
    }
}

void outputArray(int *a, size_t size) {
    for (size_t i = 0; i < size; i++) {
        printf("%d ", a[i]);
    }
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

bool isOrdered(const int *a, size_t size) {
    for (size_t i = 1; i < size - 1; i++) {
        if(a[i - 1] > a[i])
            return false;
    }
    return true;
}