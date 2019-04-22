//
// Created by huyen on 4/22/19.
//

#include <values.h>
#include <stdio.h>
#include <malloc.h>
#include "lsh_search.h"
#include "utils.h"

double search(int dim, HashBucket *bucket, double *query, double minDistance, double *result_ptr) {
    double distance = MAXDOUBLE;
    double *result = (double *) malloc(dim * sizeof(double));

    LinkedList *data = bucket->head;
    while (data->next) {
        double currentDistance = distanceOfTwoPoints(dim, data->data, query);
        if (distance > currentDistance) {
            distance = currentDistance;
            for (int i = 0; i < dim; ++i) {
                result[i] = data->data[i];
            }
        }
        data = data->next;
    }

    if (distance < minDistance) {
        for (int i = 0; i < dim; ++i) {
            result_ptr[i] = result[i];
        }
        return distance;
    }

    return minDistance;
}


int **probing(int numOfVectors, int dim, int l, int m, double w,
              double *query, double **hashFuncs) {
    int **perturbationVectors = (int **) malloc(numOfVectors * sizeof(int *));

    printf("perturbation vectors: \n");
    for (int i = 0; i < numOfVectors; ++i) {
        printf("%d -- \n", i);
        for (int j = 0; j < m; ++j) {
            printf("%d \n", perturbationVectors[i][j]);
        }
    }

    getchar();
    return perturbationVectors;
}

int LSH_search(int dim, int l, int m, double w, double ***hashTables,
                   HashBucket *buckets, double *query, double *centroid, double *distanceB4Probing, double *result) {
    double distance = MAXDOUBLE;
    int **hashVal = calculateHashValues(dim, l, m, w, centroid, hashTables, query);

    HashBucket *ite = buckets;

    while (ite != NULL) {
        if (compareHashValues(l, m, hashVal, ite->hashValues)) {
            for (int i = 0; i < l; ++i) {
                free(hashVal[i]);
            }
            free(hashVal);
            distance = search(dim, ite, query, MAXDOUBLE, result);
            break;
        }
        ite = ite->next;
    }

    *distanceB4Probing = distance;

    printf("Closest distance b4 probing: %f \n", distance);

    printDataSet(dim, 1, result);

    getchar();
    return 0;
}

int
_LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query, int NUM_VECTORS,
            double *centroid, double *distanceB4Probing, double *result) {
    int **hashVal = calculateHashValues(dim, l, m, w, centroid, hashTables, query);

//    int **perturVectors = probing(NUM_VECTORS, dim, l, m, w, query, hashTables[0]);

    LSH_search(dim, l, m, w, hashTables, buckets, query, centroid,
                        distanceB4Probing, result);

    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }

//    printf("Closest distance: %f", distanceOfTwoPoints(dim, result, query));
    free(hashVal);
//    getchar();
    return 0;
}