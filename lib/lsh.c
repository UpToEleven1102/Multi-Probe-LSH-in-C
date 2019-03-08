//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include <values.h>
#include "lsh.h"
#include "utils.h"

HashBucket *hashBuckets = NULL;

int insert(int dim, int l, int m, double w, double ***hashTables, double *ele) {
    int **hashValues = calculateHashValues(dim, l, m, w, hashTables, ele);

    if (hashBuckets == NULL) {
        hashBuckets = (HashBucket *) malloc(sizeof(HashBucket));
        hashBuckets->hashValues = hashValues;
        LinkedList *head = (LinkedList *) malloc(sizeof(LinkedList));
        head->data = ele;
        head->next = NULL;
        hashBuckets->head = head;
        hashBuckets->next = NULL;
        return 0;
    }

    HashBucket *ite = hashBuckets;

    while (ite != NULL) {
        if (compareHashValues(l, m, hashValues, ite->hashValues)) {
            LinkedList *node = (LinkedList *) malloc(sizeof(LinkedList));
            node->data = ele;
            node->next = ite->head;
            ite->head = node;
            return 0;
        }
        ite = ite->next;
    }

    HashBucket *bucket = (HashBucket *) malloc(sizeof(HashBucket));
    bucket->hashValues = hashValues;
    LinkedList *head = (LinkedList *) malloc(sizeof(LinkedList));
    head->data = ele;
    head->next = NULL;
    bucket->head = head;
    bucket->next = hashBuckets;

    hashBuckets = bucket;

    return 0;
}

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets) {
    hashBuckets = buckets;
    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        insert(dim, l, m, w, hashTables, ele);
    }

    return hashBuckets;
}

double search(int dim, HashBucket *bucket, double *query, double *result) {
    double distance = MAXDOUBLE;

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
    printf("Closest distance: %f", distance);

    return distance;
}

int probing() {

}

double *LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query) {
    double *result = (double *) malloc(sizeof(double) * dim);
    int **hashVal = calculateHashValues(dim, l, m, w, hashTables, query);
    double distance = MAXDOUBLE;

    HashBucket *ite = buckets;

    while (ite->next) {
        if (compareHashValues(l, m, hashVal, ite->hashValues)) {
            //do probing here
            double localDistance = search(dim, ite, query, result);

            distance = distance > localDistance? localDistance: distance;
            break;
        }
        ite = ite->next;
    }

    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }

    printf("Closest distance: %f", distance);

    free(hashVal);
    return result;
}



