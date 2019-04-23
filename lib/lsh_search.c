//
// Created by huyen on 4/22/19.
//

#include <values.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
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

double calculateDistanceToBucket(int dim, int l, int m, double w, int **hashVal, int **bucketHashVal, double *query, double ***hashTables, double *centroid) {
    double distance = 0;

    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < m; ++j) {
            int steps = (hashVal[i][j] - bucketHashVal[i][j]);
            if (abs(steps) > 0){
                distance += (abs(steps) - 1) * w;
            }

            distance = distance + distanceToBoundary(dim, w, query, hashTables[i][j], centroid, steps); // + distance to the boundary
        }
    }

    return distance;
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

    BucketHashVal *bucketHashVal = NULL;

    bool found = false;
    BucketHashVal *newBucketHashVal;
    while (ite != NULL) {
        newBucketHashVal = (BucketHashVal*)malloc(sizeof(BucketHashVal));
        newBucketHashVal->next = NULL;
        newBucketHashVal->prev = NULL;
        newBucketHashVal->score = calculateDistanceToBucket(dim, l, m, w, hashVal, ite->hashValues, query, hashTables, centroid);
        newBucketHashVal->value = (int **)malloc(l * sizeof(int*));

        for (int i = 0; i < l; ++i) {
            newBucketHashVal->value[i] = (int*)calloc(m, sizeof(int));
            for (int j = 0; j < m; ++j) {
                newBucketHashVal->value[i][j] = ite->hashValues[i][j];
            }
        }

        if (bucketHashVal == NULL) {
            bucketHashVal = newBucketHashVal;
        } else {
            BucketHashVal *iteHashVal = bucketHashVal;


            while(iteHashVal != NULL) {
                if (iteHashVal->next == NULL) {
                    iteHashVal->next = newBucketHashVal;
                    newBucketHashVal->prev = iteHashVal;
                    break;
                }

                if (iteHashVal->score > newBucketHashVal->score) {
                    if (iteHashVal->prev == NULL) {
                        iteHashVal->prev = newBucketHashVal;
                        newBucketHashVal->next = iteHashVal;
                        bucketHashVal = newBucketHashVal;
                        break;
                    }

                    iteHashVal->prev->next = newBucketHashVal;
                    newBucketHashVal->prev = iteHashVal ->prev;
                    newBucketHashVal->next = iteHashVal;
                    iteHashVal->prev = newBucketHashVal;
                    break;
                }
                iteHashVal = iteHashVal->next;
            }


        }

        if (!found && compareHashValues(l, m, hashVal, ite->hashValues)) {
            printf("Found %f \n", newBucketHashVal->score);
            distance = search(dim, ite, query, MAXDOUBLE, result);
            found = true;
        }
        ite = ite->next;
    }

    *distanceB4Probing = distance;
    printf("Distance b4 probing: %f \n", *distanceB4Probing);


    //probe buckets
    BucketHashVal * iteHashVal = bucketHashVal;
    int counter = 0;

    while(iteHashVal != NULL) {
        printf("%d - score - %f \n", ++counter, iteHashVal->score);
        for (int i = 0; i < l; ++i) {
            free(iteHashVal->value[i]);
        }
        free(iteHashVal->value);
        if(iteHashVal->prev!=NULL) {
            free(iteHashVal->prev);
        }
        iteHashVal = iteHashVal->next;
    }

    getchar();

    printDataSet(dim, 1, result);
    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }
    free(hashVal);
    return 0;
}