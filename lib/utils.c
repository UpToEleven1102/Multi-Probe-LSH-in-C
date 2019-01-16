//
// Created by huyen on 1/15/19.
//

#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"

float* generateDataSet(int dim, int n_data){
    float* data = (float*)malloc(sizeof(float)*dim*n_data);

    srand((unsigned int) time(NULL));

    for (int i = 0; i < dim * n_data; ++i) {
        data[i] = (float)rand()/RAND_MAX;
    }

    return data;
}

void printDataSet(int dim, int n_data, float *data) {
    int counter = 0;
    for (int i = 0; i < dim * n_data; ++i) {
        if (i %4 == 0) {
            printf("%d ---\n", counter++);
        }
        printf("%f \n", data[i]);
    }
}

