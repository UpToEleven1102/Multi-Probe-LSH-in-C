//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include "lsh.h"
#include "utils.h"

HashBucket *hashBuckets = NULL;

double **calculateHashValues(int dim, int l, int m, double w, double ***hashTables, double *ele) {
    double **hashValues = (double **) malloc(l * sizeof(double *));
    for (int i = 0; i < l; ++i) {
        hashValues[i] = (double *) malloc(m * sizeof(double));
        for (int j = 0; j < m; ++j) {
            hashValues[i][j] = calculateHashValue(dim, w, ele, hashTables[i][j]);
        }
    }

    return hashValues;
}

bool compareHashValues(int l, int m, double **hashValue1, double **hashValue2) {
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < m; ++j) {
            if (hashValue1[i][j] != hashValue2[i][j])
                return false;
        }
    }

    return true;
}

int insert(int dim, int l, int m, double w, double ***hashTables, double *ele) {
    double **hashValues = calculateHashValues(dim, l, m, w, hashTables, ele);

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
    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        insert(dim, l, m, w, hashTables, ele);
    }

    return hashBuckets;
}
