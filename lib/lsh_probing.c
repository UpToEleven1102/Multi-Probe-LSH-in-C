//
// Created by huyen on 1/21/19.
//

#include <stdio.h>
#include <stdlib.h>
#include "lsh_probing.h"


double *lshProbing(int dim, int n_data, int l, int m, double w, double*** hashTables, HashBucket *buckets, double *query, double *data){
    double **queryHashValue = calculateHashValues(dim, l, m, w, hashTables, query);

    HashBucket *ite = buckets;
    LinkedList *currentBucketHead = NULL;

    while(ite!=NULL) {
        if(compareHashValues(l, m, queryHashValue, ite->hashValues)) {

            currentBucketHead = ite->head;
            break;
        }
        ite = ite->next;
    }

    if(currentBucketHead == NULL) {
        return NULL;
    }

    double *result_ptr;
    double shortestDistance = RAND_MAX;

    while(currentBucketHead != NULL) {
        double distance = distanceOfTwoPoints(dim, currentBucketHead->data, query);
        if (distance < shortestDistance) {
            shortestDistance = distance;
            result_ptr = currentBucketHead->data;
        }
        currentBucketHead = currentBucketHead->next;
    }

    printf("shortest distance: %f \n", shortestDistance);

    return result_ptr;
}