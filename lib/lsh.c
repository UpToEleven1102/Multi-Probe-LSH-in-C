//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include <values.h>
#include "lsh.h"
#include "utils.h"
#include "data_structure.h"

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

double **probing(int numOfVectors, int dim, int l, int m, double w, double *query, double **hashFuncs) {
    double **perturbationVectors = (double **)malloc(numOfVectors * sizeof(double*));
    //find 2M array
    struct pairZ twoM[2*m];

    int counter = 0;
    for (int i = 0; i < m; ++i) {
        struct pairZ z1 = { .x = distanceToBoundary(dim, w, query, hashFuncs[i], -1), .i = i, .r = -1};
        twoM[counter++] = z1;
        struct pairZ z2 = { .x = distanceToBoundary(dim, w, query, hashFuncs[i], 1), .i = i, .r = 1};
        twoM[counter++] = z2;
    }

    for (int i = 0; i < 2 * m; ++i) {
        printf("%f - %d \n", twoM[i].x, twoM[i].i);
    }

    for (int i = 0; i < 2 * m; ++i) {
        for (int j = i+1; j < 2 * m - 1; ++j) {
            struct pairZ temp = twoM[i];
            if (twoM[j].x < temp.x) {
                twoM[i] = twoM[j];
                twoM[j] = temp;
            }
        }
    }

    //generate Heap
//    struct PerturVector *a0 = &(struct PerturVector){.data = 0, .score = 0, .length = 1, .next = NULL, .prev = NULL, .calculateScore = calculateScoreA, .isValid = isValidA, .shift = shiftA, .expand = expandA};
//    a0 ->calculateScore(a0, twoM);
//    struct Heap *heap = &(struct Heap){.head = a0, .add = addHeapLinkedList, .remove = removeHeapLinkedList, .extractMin = extractMinHeapLinkedList};

    //generate perturbation vectors
//    for (int i = 0; i < numOfVectors; ++i) {
//        perturbationVectors[i] = (double*) malloc(dim * sizeof(double));
//        struct PerturVector *minA;
//        do {
//            minA = heap.extractMin(&heap);
//            struct PerturVector *shifted = minA->shift(minA);
//            shifted->calculateScore(shifted, twoM);
//            heap.add(&heap, shifted);
//            struct PerturVector *expanded = minA->expand(minA);
//            expanded->calculateScore(expanded, twoM);
//            heap.add(&heap, expanded);
//        } while(!minA->isValid(minA, 2*m));
//
//        //extract minA to perturbationVectors[i]
//    }

    return perturbationVectors;
}

double *LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query) {
    double *result = (double *) malloc(sizeof(double) * dim);
    int **hashVal = calculateHashValues(dim, l, m, w, hashTables, query);
    double distance = MAXDOUBLE;
    const int NUM_VECTORS = 50;

    HashBucket *ite = buckets;

    while (ite->next) {
        if (compareHashValues(l, m, hashVal, ite->hashValues)) {
            //do probing here

            double **perturVectors = probing(NUM_VECTORS, dim, l, m, w, query, hashTables[0]);

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



