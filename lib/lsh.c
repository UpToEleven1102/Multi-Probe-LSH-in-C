//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include <values.h>
#include "lsh.h"
#include "utils.h"

HashBucket *hashBuckets = NULL;

int insert(int dim, int l, int m, double w, double ***hashTables, double *ele, double *centroid, int *num_hash_buckets) {
    int **hashValues = calculateHashValues(dim, l, m, w, centroid, hashTables, ele);

    if (hashBuckets == NULL) {
        *num_hash_buckets = *num_hash_buckets + 1;
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

    *num_hash_buckets = *num_hash_buckets + 1;

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

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets,
                double *centroid, int* num_hash_buckets) {
    *num_hash_buckets = 0;
    hashBuckets = buckets;
    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        insert(dim, l, m, w, hashTables, ele, centroid, num_hash_buckets);
    }

    return hashBuckets;
}
