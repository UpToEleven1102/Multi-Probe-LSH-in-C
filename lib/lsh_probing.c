//
// Created by huyen on 1/21/19.
//

#include <stdio.h>
#include <stdlib.h>
#include "lsh_probing.h"
#include "data_structure.h"


void shift(int length, double *set) {
    set[length - 1]++;
}

double *expand(int length, double *set) {
    double *newSet = (double *) malloc((length + 1) * sizeof(double));
    for (int i = 0; i < length; ++i) {
        newSet[i] = set[i];
    }

    free(set);
    return newSet;
}

bool validA() {
    return true;
}

int **generatePerturbationVectors(int dim, int m, double w, int t, double *query, double **hashTable) {
    int counter = 0;
    int **perturbationSets = (int **) malloc(t * sizeof(int *));

    struct Z *zs = (struct Z*)malloc(2*m*sizeof(struct Z));

    for (int i = 0; i < m; ++i) {
        zs[i].x = distanceToBoundary(dim, w, query, hashTable[i], -1);
        zs[i].i = i;
        zs[i].r = -1;
        zs[i+m].x = distanceToBoundary(dim, w, query, hashTable[i], 1);
        zs[i+m].i = i;
        zs[i+m].r = 1;
    }

    for (int i = 0; i < 2*m; ++i) {
        for (int j = i+1; j < 2 * m; ++j) {
            if (zs[i].x > zs[j].x) {
                struct Z temp = zs[i];
                zs[i] = zs[j];
                zs[j] = temp;
            }
        }
    }

    //Pi j is the pair of i and r at index j

    //todo: write struct for A

    struct Heap heap = {NULL, addHeapLinkedList, removeHeapLinkedList};

    //initialize first perturbation set
    struct A *head = (struct A*)malloc(sizeof(struct A));
    head->data = (double*)malloc(sizeof(double));
    head->next=NULL;
    head->length = 1;
    head->data[0] = 1;

    heap.add(&heap, head);



    for (int i = 0; i < t; ++i) {
        while(1) {
            if(validA())
                break;
        }
    }

    return perturbationSets;
}


double *
lshProbing(int dim, int n_data, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query,
           double *data) {
    double **queryHashValue = calculateHashValues(dim, l, m, w, hashTables, query);

    int **perturbationVector = (int **) malloc(3 * sizeof(int *));

    for (int i = 0; i < 3; ++i) {
        perturbationVector[i] = (int *) malloc(m * sizeof(double));
        switch (i) {
            case 0:
                for (int j = 0; j < m; ++j)
                    perturbationVector[i][j] = -1;
                break;
            case 1:
                for (int j = 0; j < m; ++j)
                    perturbationVector[i][j] = 0;
                break;
            case 2:
                for (int j = 0; j < m; ++j)
                    perturbationVector[i][j] = 1;
            default:
                //do nothing
                break;
        }
    }

    HashBucket *ite = buckets;
    LinkedList *currentBucketHead = NULL;

    while (ite != NULL) {
        if (compareHashValues(l, m, queryHashValue, ite->hashValues)) {

            currentBucketHead = ite->head;
            break;
        }
        ite = ite->next;
    }

    if (currentBucketHead == NULL) {
        return NULL;
    }

    double *result_ptr;
    double shortestDistance = RAND_MAX;

    while (currentBucketHead != NULL) {
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