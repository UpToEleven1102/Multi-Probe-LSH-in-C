#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lib/utils.h"
#include "lib/lsh.h"
#include "lib/lsh_probing.h"

//TODO: more research about number L and M, how many are needed
//TODO: implement b
// how large we should choose value of L, paper
// hash functions are unit vectors?
// hashValue formula   <a.v -b> / w. what is b??

//b = centroid * hash functions??????? data dependent

//W is dependent to the number of data points

//slides: intuitions wrong: 1-step buckets are not better than 2-step buckets, paper: 1step buckets are not neccessary better than buckets that are 2 steps away

void initParameters(int *L, int *M, double *W, int dim, int n_data, const double *data) {
    //comeback and pick this up later
    *M = (int) floor(dim / 2.0);

    *L = *M;

    double **buff = (double **) malloc(dim * sizeof(double *));

    for (int i = 0; i < dim; ++i) {
        buff[i] = (double *) malloc(2 * sizeof(double));
        buff[i][0] = 0;
        buff[i][1] = RAND_MAX;
    }

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        for (int j = 0; j < dim; ++j) {
            if (ele[j] > buff[j][0]) {
                buff[j][0] = ele[j];
            }
            if (ele[j] < buff[j][1]) {
                buff[j][1] = ele[j];
            }
        }
        free(ele);
    }

    double maxDistance = 0;
    for (int i = 0; i < dim; ++i) {
        if (maxDistance < buff[i][0] - buff[i][1]) {
            maxDistance = buff[i][0] - buff[i][1];
        }
    }

    *W = maxDistance / 2;

    for (int i = 0; i < dim; ++i) {
        free(buff[i]);
    }
    free(buff);
}

int main() {
    const int dim = 8;
    const int n_data = 1000;
    double *data = generateDataSet(dim, n_data);

//    srand((unsigned int) time(NULL));
    srand(3);
    int *L = (int *) malloc(sizeof(int));
    int *M = (int *) malloc(sizeof(int));
    double *W = (double *) malloc(sizeof(double));

    initParameters(L, M, W, dim, n_data, data);

    printf("L - %d, M - %d, W - %f, dim - %d \n", *L, *M, *W, dim);
    double ***hashTables = generateHashTables(*L, *M, dim);

    HashBucket *buckets = LSH(dim, n_data, *L, *M, *W, hashTables, data);

    printf("hash buckets: \n");
    printHashBuckets(dim, *L, *M, buckets);

    double *query = generateDataSet(dim, 1);


    //start lsh_probing
    double *result = lshProbing(dim, n_data, *L, *M, *W, hashTables, buckets, query, data);

    printf("Query point: \n");
    printDataSet(dim, 1, query);

//    printf("Distance of query to data points in data set: \n");

    int closestIdx = 0;
    double closestDistance = RAND_MAX;

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        double distance = distanceOfTwoPoints(dim, query, ele);
        if (distance < closestDistance) {
            closestDistance = distance;
            closestIdx = i;
        }
//        printf("data %d: %f \n", i, distance);
        free(ele);
    }

    if (result == NULL) {
        printf("result = NULL");
    } else {
        printf("Closest data point: \n");
        printDataSet(dim, 1, result);
    }

    printf("Closest idx: %d - distance: %f \n", closestIdx, closestDistance);

    //verify variables
//    printHashTables(dim, *L, *M, hashTables);
//    printDataSet(dim, n_data, data);


    //free pointer variables
    HashBucket *ite = buckets;
    while (ite != NULL) {
        HashBucket *temp = ite;
        ite = ite->next;
        for (int i = 0; i < *L; ++i) {
            free(temp->hashValues[i]);
        }
        free(temp->hashValues);
        LinkedList *listIte = temp->head;

        while (listIte != NULL) {
            LinkedList *tempListIte = listIte;
            listIte = listIte->next;
            free(tempListIte->data);
            free(tempListIte);
        }

        free(temp);
    }

    for (int i = 0; i < *L; ++i) {
        for (int j = 0; j < *M; ++j) {
            free(hashTables[i][j]);
        }
        free(hashTables[i]);
    }
    free(hashTables);
    free(L);
    free(M);
    free(W);
    free(data);
}