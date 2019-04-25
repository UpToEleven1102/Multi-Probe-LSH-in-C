//
// Created by huyen on 4/22/19.
//

#include <values.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "lsh_search.h"
#include "utils.h"

double search(int dim, HashBucket *bucket, double *query, double minDistance, double *result_ptr, int* counter) {
    double distance = MAXDOUBLE;
    double *result = (double *) malloc(dim * sizeof(double));

    LinkedList *ite = bucket->head;
    while (ite != NULL) {
        *counter = *counter +1;
        double currentDistance = distanceOfTwoPoints(dim, ite->data, query);
        if (distance > currentDistance) {
            distance = currentDistance;
            for (int i = 0; i < dim; ++i) {
                result[i] = ite->data[i];
            }
        }
        ite = ite->next;
    }

    if (distance < minDistance) {
        for (int i = 0; i < dim; ++i) {
            result_ptr[i] = result[i];
        }
        free(result);
        return distance;
    }
    free(result);
    return minDistance;
}

int
_LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, int num_buckets, double *query,
            double *centroid, double *distanceB4Probing, double *result, int *num_checked_data_points) {
    double distance = MAXDOUBLE;
    *num_checked_data_points = 0;
    int num_probing_buckets = num_buckets * 5 / 100;
    printf("%d ----\n", num_probing_buckets);

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
            distance = search(dim, ite, query, MAXDOUBLE, result, num_checked_data_points);
            found = true;
        }
        ite = ite->next;
    }

    *distanceB4Probing = distance;
    printf("Distance b4 probing: %f, checked: %d \n", *distanceB4Probing, *num_checked_data_points);


    int *num_checked = (int*)malloc(sizeof(int));

    *num_checked = *num_checked_data_points;

    //probe buckets
    BucketHashVal *iteHashVal = bucketHashVal->next;
    int counter = 0;
    while(iteHashVal != NULL) {
        if (counter++<num_probing_buckets) {
            ite = buckets;
            while(ite != NULL) {
                if (compareHashValues(l, m, iteHashVal->value, ite->hashValues)) {
                    distance = search(dim, ite, query, distance, result, num_checked);
                    break;
                }
                ite = ite->next;
            }
        }

//        printf("%d - score - %f \n", ++counter, iteHashVal->score);
        for (int i = 0; i < l; ++i) {
            free(iteHashVal->value[i]);
        }
        free(iteHashVal->value);
        if(iteHashVal->prev!=NULL) {
            free(iteHashVal->prev);
        }
        iteHashVal = iteHashVal->next;
    }

    printf("distance after probing: %f, checked : %d \n", distance, *num_checked);

//    getchar();

    printDataSet(dim, 1, result);
    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }
    free(hashVal);
    return *num_checked;
}