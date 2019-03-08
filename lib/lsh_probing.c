////
//// Created by huyen on 1/21/19.
////
//
//#include <stdio.h>
//#include <stdlib.h>
//#include "lsh_probing.h"
//#include "data_structure.h"
//
//int **generatePerturbationVectors(int dim, int m, double w, int t, double *query, double **hashTable) {
//    int counter = 0;
//    int **perturbationSets = (int **) malloc(t * sizeof(int *));
//
//    struct Z *zs = (struct Z *) malloc(2 * m * sizeof(struct Z));
//
//    for (int i = 0; i < m; ++i) {
//        zs[i].x = distanceToBoundary(dim, w, query, hashTable[i], -1);
//        zs[i].i = i;
//        zs[i].r = -1;
//        zs[i + m].x = distanceToBoundary(dim, w, query, hashTable[i], 1);
//        zs[i + m].i = i;
//        zs[i + m].r = 1;
//    }
//
//    for (int i = 0; i < 2 * m; ++i) {
//        for (int j = i + 1; j < 2 * m; ++j) {
//            if (zs[i].x > zs[j].x) {
//                struct Z temp = zs[i];
//                zs[i] = zs[j];
//                zs[j] = temp;
//            }
//        }
//    }
//
//    //Pi j is the pair of i and r at index j
//
//    //todo: write struct for PerturVector
//
//    struct Heap heap = {NULL, addHeapLinkedList, removeHeapLinkedList, extractMinHeapLinkedList};
//
//    //initialize first perturbation set
//
//
//    struct PerturVector head = {(int[]) {1}, 0, 1, NULL, calculateScoreA, isValidA, shiftA, expandA};
////    struct PerturVector *head = (struct PerturVector*)malloc(sizeof(struct PerturVector));
////    head->data = (double*)malloc(sizeof(double));
////    head->next=NULL;
////    head->length = 1;
////    head->data[0] = 1;
////    head->calculateScore = calculateScoreA;
////    head->isValid = isValidA;
//
//    heap.add(&heap, &head);
//
//    while(counter < t) {
//        struct PerturVector *topNode = heap.extractMin(&heap);
//
//        //get top node
//        // generate 2 new PerturVector sets
//        //output the valid top node.
//        //there is a stoping rule.......?????
//
//        struct PerturVector *shifted = topNode->shift(topNode);
//        struct PerturVector *expanded = topNode->expand(topNode);
//
//        heap.add(&heap, shifted);
//        heap.add(&heap, expanded);
//
//        if (topNode->isValid(topNode, 2 * m)) {
//            perturbationSets[counter] = topNode->data;
//            counter++;
//        }
//    }
//
//    return perturbationSets;
//}
//
//double *
//lshProbing(int dim, int n_data, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query,
//           double *data) {
//    int **queryHashValue = calculateHashValues(dim, l, m, w, hashTables, query);
//
//    int **perturbationVector = (int **) malloc(3 * sizeof(int *));
//
//    for (int i = 0; i < 3; ++i) {
//        perturbationVector[i] = (int *) malloc(m * sizeof(double));
//        switch (i) {
//            case 0:
//                for (int j = 0; j < m; ++j)
//                    perturbationVector[i][j] = -1;
//                break;
//            case 1:
//                for (int j = 0; j < m; ++j)
//                    perturbationVector[i][j] = 0;
//                break;
//            case 2:
//                for (int j = 0; j < m; ++j)
//                    perturbationVector[i][j] = 1;
//            default:
//                //do nothing
//                break;
//        }
//    }
//
//    HashBucket *ite = buckets;
//    LinkedList *currentBucketHead = NULL;
//
//    while (ite != NULL) {
//        if (compareHashValues(l, m, queryHashValue, ite->hashValues)) {
//
//            currentBucketHead = ite->head;
//            break;
//        }
//        ite = ite->next;
//    }
//
//    if (currentBucketHead == NULL) {
//        return NULL;
//    }
//
//    double *result_ptr = NULL;
//    double shortestDistance = RAND_MAX;
//
//    while (currentBucketHead != NULL) {
//        double distance = distanceOfTwoPoints(dim, currentBucketHead->data, query);
//        if (distance < shortestDistance) {
//            shortestDistance = distance;
//            result_ptr = currentBucketHead->data;
//        }
//        currentBucketHead = currentBucketHead->next;
//    }
//
//    printf("shortest distance: %f \n", shortestDistance);
//
//    return result_ptr;
//}