//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include "lsh.h"
#include "utils.h"

HashBucket *hashBuckets = NULL;

int insert(int dim, int l, int m, double w, double ***hashTables, double *ele, const double *centroid) {
    double **hashValues = calculateHashValues(dim, l, m, w, centroid, hashTables, ele);

    if (hashBuckets == NULL) {
        hashBuckets = (HashBucket*)malloc(sizeof(HashBucket));
        hashBuckets->hashValues = hashValues;
        LinkedList *head = (LinkedList*)malloc(sizeof(LinkedList));
        head->data = ele;
        head->next = NULL;
        hashBuckets->head = head;
        return 0;
    }

    HashBucket *ite = hashBuckets;

    while (ite != NULL) {
        if (compareHashValues(l, m, hashValues, ite-> hashValues)) {
            LinkedList *node = (LinkedList*)malloc(sizeof(LinkedList));
            node->data = ele;
            node->next = ite->head;
            ite->head = node;
            return 0;
        }
        ite = ite->next;
    }

    HashBucket *bucket = (HashBucket*)malloc(sizeof(HashBucket));
    bucket->hashValues = hashValues;
    LinkedList *head = (LinkedList*)malloc(sizeof(LinkedList));
    head->data = ele;
    head->next = NULL;
    bucket->head = head;
    bucket->next = hashBuckets;

    hashBuckets=bucket;

    return 0;
}

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data) {
    const double *centroid = calculateCentroid(dim, n_data, data);

    printf("Centroid: \n");
    printDataSet(dim, 1, centroid);

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        insert(dim, l, m, w, hashTables, ele, centroid);
    }

    return hashBuckets;
}
