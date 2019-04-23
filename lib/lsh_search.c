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

int
_LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, int num_buckets, double *query,
            double *centroid, double *distanceB4Probing, double *result) {
    double distance = MAXDOUBLE;
    int **hashVal = (int**)malloc(l * sizeof(int*));

    for (int i = 0; i < l; ++i) {
        hashVal[i] = (int*)calloc(m,sizeof(int));
    }

    calculateHashValues(dim, l, m, w, centroid, hashTables, query, hashVal);

    HashBucket *ite = buckets;

    while (ite != NULL) {
        if (compareHashValues(l, m, hashVal, ite->hashValues)) {
            distance = search(dim, ite, query, MAXDOUBLE, result);
            break;
        }
        ite = ite->next;
    }

    *distanceB4Probing = distance;

    printf("Closest distance b4 probing: %f \n", distance);

    printDataSet(dim, 1, result);
    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }
    free(hashVal);
    return 0;
}