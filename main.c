#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lib/utils.h"

//TODO: more research about number L and M, how many are needed

void initParameters(int *L, int *M, double *W, int dim, int n_data, const double *data) {
    //comeback and pick this up later
    *M = (int)floor(dim/2.0);

    *L = *M;

    double **buff = (double**) malloc(dim * sizeof(double *));

    for (int i = 0; i < dim; ++i) {
        buff[i] = (double *) malloc(2 * sizeof(double));
        buff[i][0] = 0;
        buff[i][1] = RAND_MAX;
    }

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        for (int j = 0; j < dim; ++j) {
            if (ele[j] > buff[j][0]) {
                buff[j][0] = ele[j];
            }
            if (ele[j] < buff[j][1]) {
                buff[j][1] = ele[j];
            }
        }
        free(ele);
    }

    double maxDistance = 0;
    for (int i = 0; i < dim; ++i) {
        if (maxDistance < buff[i][0] - buff[i][1]) {
            maxDistance = buff[i][0] - buff[i][1];
        }
    }

    *W = maxDistance / 4;

    for (int i = 0; i < dim; ++i) {
        free(buff[i]);
    }
    free(buff);
}

int main() {
    const int dim = 8;
    const int n_data = 100;
    double *data = generateDataSet(dim, n_data);

    srand((unsigned int) time(NULL));

    int *L = (int *) malloc(sizeof(int));
    int *M = (int *) malloc(sizeof(int));
    double *W = (double *) malloc(sizeof(double));

    initParameters(L, M, W, dim, n_data, data);

    printf("L - %d, M - %d, W - %f, dim - %d \n", *L, *M, *W, dim);
    double ***hashTables = generateHashTables(*L, *M, dim);

    printHashTables(dim, *L, *M, hashTables);
//    printDataSet(dim, n_data, data);


    //free pointer variables
    for (int i = 0; i < *L; ++i) {
        for (int j = 0; j < *M; ++j) {
            free(hashTables[i][j]);
        }
        free(hashTables[i]);
    }
    free(hashTables);
    free(L); free(M); free(W);free(data);

}