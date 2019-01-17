#include <stdio.h>
#include <malloc.h>
#include "lib/utils.h"

void initParameters(int *L, int *M, double *W, int dim, int n_data, const double *data) {
    double max = 0, min = 0;

    for (int i = 0; i < n_data; ++i) {
// create getElementmethod
    }
}

int main() {
    const int dim = 4;
    const int n_data = 20;
    double *data = generateDataSet(dim, n_data);

    int *L = (int*)malloc(sizeof(int));
    int *M = (int*)malloc(sizeof(int));
    double *W = (double*)malloc(sizeof(double));

//    initParameters(L, M, W, dim, n_data, data);
//    double ***hashTables = generateHashTables(*L, *M, dim);

    printDataSet(dim,n_data, data);
}