//
// Created by huyen on 1/15/19.
//

#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"

double *generateDataSet(int dim, int n_data) {
    double *data = (double *) malloc(sizeof(double) * dim * n_data);

    srand((unsigned int) time(NULL));

    for (int i = 0; i < dim * n_data; ++i) {
        data[i] = (double) rand() / RAND_MAX;
    }

    return data;
}

double *newUnitVector(int dim) {
    double *unitVector = generateDataSet(dim, 1);

    double vectorLength = 0;
    for (int i = 0; i < dim; ++i) {
        vectorLength += unitVector[i] * unitVector[i];
    }

    vectorLength = sqrt(vectorLength);

    for (int i = 0; i < dim; ++i) {
        unitVector[i] = unitVector[i] / vectorLength;
    }

    return unitVector;
}

double **generateHashTable(int m, int dim) {
    double **h = (double **) malloc(m * sizeof(double));

    for (int i = 0; i < m; ++i) {
        h[i] = newUnitVector(dim);
    }

    return h;
}

double *getElementAtIndex(int idx, int dim, int n_data, const double *data) {
    double *ele = (double *) malloc(dim * sizeof(double));

    if (idx < n_data)
        for (int i = 0; i < dim; ++i) {
            ele[i] = data[idx * dim + i];
        }
    return ele;
}

double ***generateHashTables(int l, int m, int dim) {
    double ***hashTables = (double ***) malloc(l * sizeof(double **));

    for (int i = 0; i < l; ++i) {
        hashTables[i] = generateHashTable(m, dim);
    }
}

void printDataSet(int dim, int n_data, const double *data) {
    int counter = 0;
    for (int i = 0; i < dim * n_data; ++i) {
        if (i % 4 == 0) {
            printf("%d ---\n", counter++);
        }
        printf("%f \n", data[i]);
    }
}

