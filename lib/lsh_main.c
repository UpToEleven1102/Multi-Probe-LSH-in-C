//
// Created by huyen on 3/29/19.
//

#include <stdlib.h>
#include <stdio.h>
#include <values.h>
#include <math.h>
#include "lsh_main.h"
#include "utils.h"
#include "lsh.h"

double *generateDataSet(int dim, int n_data) {
    double *data = (double *) malloc(sizeof(double) * dim * n_data);

    for (int i = 0; i < dim * n_data; ++i) {
        data[i] = (double) rand() / RAND_MAX;
    }

    return data;
}

double z1;

//what mean and deviation that want?

//phase = 0 or 1 equivalent to z0 or z1
//Box-Mulller transformation
double gaussian_rand(double mean, double stdDev, int phase) {
    const double epsilon = 2.22507e-308;
    const double two_pi = 2 * M_PI;

    printf("%f \n", epsilon);

    if (phase == 1)
        return z1 * stdDev + mean;

    double u1, u2;
    do {
        u1 = (rand() + 1.) / (RAND_MAX + 2.);
        u2 = rand() / (RAND_MAX + 1.);
    } while (u1 <= epsilon);

    double z0;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

    return z0 * stdDev + mean;
}

double *newUnitVector(int dim, double *mean, double *stdDev) {
//    double *unitVector = generateDataSet(dim, 1);
//

    double *unitVector = (double *) malloc(dim * sizeof(double));

    int phase = 0;

    for (int i = 0; i < dim; ++i) {
        unitVector[i] = gaussian_rand(mean[i], stdDev[i], phase);
        phase = 1 - phase;
    }

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

double ***generateHashTables(int l, int m, int dim, double *mean, double *stdDev) {
    double ***hashTables = (double ***) malloc(l * sizeof(double **));

    for (int i = 0; i < l; ++i) {
        hashTables[i] = (double **) malloc(m * sizeof(double *));
        for (int j = 0; j < m; ++j) {
            hashTables[i][j] = newUnitVector(dim, mean, stdDev);
        }
    }

    return hashTables;
}

void initParameters(int *L, int *M, double *W, double *mean, double *stdDev, int dim, int n_data, const double *data,
                    double *centroid) {
    //comeback and pick this up later
    *M = (int) floor(dim / 2.0);

    *L = 1;
    double **buff = (double **) malloc(dim * sizeof(double *));

    for (int i = 0; i < dim; ++i) {
        buff[i] = (double *) malloc(3 * sizeof(double));
        buff[i][0] = 0;
        buff[i][1] = RAND_MAX;
        buff[i][2] = 0;
        stdDev[i] = 0;
    }

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);

        for (int j = 0; j < dim; ++j) {
            centroid[j] += ele[j];

            if (ele[j] > buff[j][0]) {
                buff[j][0] = ele[j];
            }
            if (ele[j] < buff[j][1]) {
                buff[j][1] = ele[j];
            }
            buff[j][2] += ele[j];
        }
        free(ele);
    }

    double maxDistance = 0;
    for (int i = 0; i < dim; ++i) {
        centroid[i] /= n_data;
        mean[i] = buff[i][2] / n_data;
        if (maxDistance < buff[i][0] - buff[i][1]) {
            maxDistance = buff[i][0] - buff[i][1];
        }
    }

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);

        for (int j = 0; j < dim; ++j) {
            stdDev[j] = (ele[j] - mean[j]) * (ele[j] - mean[j]);
        }

        free(ele);
    }

    double _mean = 0;

    for (int i = 0; i < dim; ++i) {
        stdDev[i] = sqrt(stdDev[i] / n_data);
        _mean += mean[i];
    }

    _mean = _mean / dim;

    *W = _mean * 1.5;

    for (int i = 0; i < dim; ++i) {
        free(buff[i]);
    }
    free(buff);
}



int LSH_main(int dim, int n_data, double *data,
             double ***hashTables, HashBucket *buckets, double *centroid, double *result,
             double *datum) {
    int *L = (int *) malloc(sizeof(int));
    int *M = (int *) malloc(sizeof(int));
    double *W = (double *) malloc(sizeof(double));
    const int NUM_PERTURBATION_VECTORS = 100;
    centroid = (double *) calloc(dim, sizeof(double));

    double *mean = (double *) malloc(dim * sizeof(double));
    double *stdDev = (double *) malloc(dim * sizeof(double));

    initParameters(L, M, W, mean, stdDev, dim, n_data, data, centroid);
//    printf("L - %d, M - %d, W - %f, dim - %d \n", *L, *M, *W, dim);

    hashTables = generateHashTables(*L, *M, dim, mean, stdDev);
    printHashTables(dim, *L, *M, hashTables);

    *buckets = *LSH(dim, n_data, *L, *M, *W, hashTables, data, NULL, centroid);

    printf("hash buckets: \n");

    int numBuckets = printHashBuckets(dim, *L, *M, buckets);

    printf("number of buckets: %d \n Enter to continue \n", numBuckets);

    getchar();

    printf("Query point: \n");
    printDataSet(dim, 1, datum);

    result = LSH_probing(dim, *L, *M, *W, hashTables, buckets, datum, NUM_PERTURBATION_VECTORS, centroid);

    //second chunk comes in
//    buckets = LSH(dim, n_data, *L, *M, *W, hashTables, dataSets[1], buckets, centroid);
//    numBuckets = printHashBuckets(dim, *L, *M, buckets);

//    printf("number of buckets: %d \n Enter to continue \n", numBuckets);

    printf("Result: \n");
    printDataSet(dim, 1, result);
    getchar();

//    generatePerturbationVectors(dim, *M, *W,
//                                5, query, hashTables[0]);

//    printf("Distance of query to data points in data set: \n");

//verify distance
    int closestIdx = 0;
    double closestDistance = MAXDOUBLE;



    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        double distance = distanceOfTwoPoints(dim, datum, ele);
        if (distance < closestDistance) {
            closestDistance = distance;
            closestIdx = i;
            for (int j = 0; j < dim; ++j) {
                result[j] = ele[j];
            }
        }
        printf("data %d: %f \n", i, distance);
        free(ele);
    }

    if (result == NULL) {
        printf("result = NULL");
    } else {
        printf("Closest data point: \n");
        printDataSet(dim, 1, result);
        closestDistance = distanceOfTwoPoints(dim, datum, result);
    }

    printf("Closest idx: %d - distance: %f \n", closestIdx, closestDistance);

    HashBucket *ite = buckets;
    while (ite != NULL) {
        HashBucket *temp = ite;
        ite = ite->next;
        for (
                int i = 0;
                i < *
                        L;
                ++i) {
            free(temp
                         ->hashValues[i]);
        }
        free(temp
                     ->hashValues);
        LinkedList *listIte = temp->head;

        while (listIte != NULL) {
            LinkedList *tempListIte = listIte;
            listIte = listIte->next;
            free(tempListIte
                         ->data);
            free(tempListIte);
        }

        free(temp);
    }

    for (
            int i = 0;
            i < *
                    L;
            ++i) {
        for (
                int j = 0;
                j < *
                        M;
                ++j) {
            free(hashTables[i][j]);
        }
        free(hashTables[i]);
    }
    free(hashTables);
    free(L);
    free(M);
    free(W);
}
