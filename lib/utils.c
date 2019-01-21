//
// Created by huyen on 1/15/19.
//

#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"


double distanceOfTwoPoints(int dim, double *point1, double *point2){
    double distance = 0;

    for (int i = 0; i < dim; ++i) {
        distance += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }

    return sqrt(distance);
}

double innerProduct(const double *vector1, const double *vector2, int dim) {
    double product = 0;
    for (int i = 0; i < dim; ++i) {
        product += vector1[i] * vector2[i];
    }
    return product;
}

double *generateDataSet(int dim, int n_data) {
    double *data = (double *) malloc(sizeof(double) * dim * n_data);


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

    return hashTables;
}

void printDataSet(int dim, int n_data, const double *data) {
    int counter = 0;
    for (int i = 0; i < dim * n_data; ++i) {
        if (i % dim == 0) {
            printf("%d ---\n", counter++);
        }
        printf("%f \n", data[i]);
    }
}

double calculateHashValue(int dim, double w, double *ele, double *hashFunc) {
    double hashValue = innerProduct(ele, hashFunc, dim) / w;
    return floor(hashValue);
}

void printHashTables(int dim, int l, int m, double ***tables) {
    for (int i = 0; i < l; ++i) {
        printf("table %d \n", i);
        for (int j = 0; j < m; ++j) {
            printf("function # %d \n", j);
            printDataSet(dim, 1, tables[i][j]);
        }
    }
}

void printHashBuckets(int dim, int l, int m, HashBucket *buckets) {
    HashBucket *ite = buckets;
    int counter = 0;

    printf("-- parameters: dim %d l %d m %d \n", dim, l, m);

    while (ite != NULL) {
        printf("Bucket %d \n Hash Values: -- \n", counter++);

        for (int i = 0; i < l; ++i) {
            printf("Table %d -- \n", i);
            for (int j = 0; j < m; ++j) {
                printf("e %f \n", ite->hashValues[i][j]);
            }
        }


        printf("Elements: \n");

        LinkedList *listIte = ite->head;

        while(listIte != NULL) {
            printDataSet(dim, 1, listIte->data);
            listIte = listIte->next;
        }

        free(listIte);
        ite = ite->next;
    }

    free(ite);
}
