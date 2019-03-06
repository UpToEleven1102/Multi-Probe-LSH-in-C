#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
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

//note: get the lowest score of the heap or get the top element

//construct perturbation vectors in each hashtable??

//know when and what table to apply the perturbation vector to


double *generateDataSet(int dim, int n_data) {
    double *data = (double *) malloc(sizeof(double) * dim * n_data);

    for (int i = 0; i < dim * n_data; ++i) {
        data[i] = (double) rand() / RAND_MAX;
    }

    return data;
}


double *newUnitVector(int dim) {
    double *unitVector = generateDataSet(dim, 1);

    double vectorLength = 0;
    for (int i = 0; i < dim; ++i) {
        vectorLength += unitVector[i] * unitVector[i];
    }

    vectorLength = sqrt(vectorLength);

    for (int i = 0; i < dim; ++i) {
        unitVector[i] = unitVector[i] / vectorLength;
    }

    return unitVector;
}


double **generateHashTable(int m, int dim) {
    double **h = (double **) malloc(m * sizeof(double));

    for (int i = 0; i < m; ++i) {
        h[i] = newUnitVector(dim);
    }

    return h;
}

double ***generateHashTables(int l, int m, int dim) {
    double ***hashTables = (double ***) malloc(l * sizeof(double **));

    for (int i = 0; i < l; ++i) {
        hashTables[i] = generateHashTable(m, dim);
    }

    return hashTables;
}

void initParameters(int *L, int *M, double *W, int dim, int n_data, const double *data) {
    //comeback and pick this up later
    *M = (int) floor(dim / 2.0);

//    *L = *M;
    *L = 1;
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
    srand(0);
    const int dim = 29;
    const int n_data = 1000;
    double *data = (double *) malloc(dim * n_data * sizeof(double));

    double *data2 = (double *) malloc(dim * n_data * sizeof(double));

    double *data3 = (double *) malloc(dim * n_data * sizeof(double));

    double *query, *result;
    query = (double *) malloc(dim * sizeof(double));

    FILE *file = fopen("../data_sets/HIGGS.csv", "rb");

    char line[1024];
    int counter = 0;

    //data 1
    for (int i = 0; (fscanf(file, "%s", line) == 1); ++i) {
        const char *tok;

        for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",")) {
            data[counter] = strtof(tok, NULL);
            counter++;
        }

        if (i == n_data)
            break;
    }

    //data 2
    counter = 0;
    for (int i = 0; (fscanf(file, "%s", line) == 1); ++i) {
        const char *tok;

        for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",")) {
            data2[counter] = strtof(tok, NULL);
            counter++;
        }

        if (i == n_data)
            break;
    }

    counter = 0;
    for (int i = 0; (fscanf(file, "%s", line) == 1); ++i) {
        const char *tok;

        for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",")) {
            data3[counter] = strtof(tok, NULL);
            counter++;
        }

        if (i == n_data)
            break;
    }

    //get a datapoint as a query
    counter = 0;
    fscanf(file, "%s", line);
    const char *tok;

    for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",")) {
        query[counter] = strtof(tok, NULL);
        counter++;
    }

    fclose(file);
//    printDataSet(dim, n_data, data);

    int *L = (int *) malloc(sizeof(int));
    int *M = (int *) malloc(sizeof(int));
    double *W = (double *) malloc(sizeof(double));

    initParameters(L, M, W, dim, n_data, data
    );
//    printf("L - %d, M - %d, W - %f, dim - %d \n", *L, *M, *W, dim);

    double ***hashTables = generateHashTables(*L, *M, dim);
    printHashTables(dim, *L, *M, hashTables);

    getchar();

    HashBucket *buckets = LSH(dim, n_data, *L, *M, *W, hashTables, data, NULL);

    printf("hash buckets: \n");

    int numBuckets = printHashBuckets(dim, *L, *M, buckets);

    printf("number of buckets: %d \n", numBuckets);
    getchar();

// TODO: classify other data sets
//    buckets = LSH(dim, n_data, *L, *M, *W, hashTables, data2, buckets);
//
//    numBuckets = printHashBuckets(dim, *L, *M, buckets);
//
//    printf("number of buckets: %d \n", numBuckets);
//    getchar();
//
//    buckets = LSH(dim, n_data, *L, *M, *W, hashTables, data3, buckets);
//
//    numBuckets = printHashBuckets(dim, *L, *M, buckets);
//
//    printf("number of buckets: %d \n", numBuckets);
//    getchar();

    printf("Query point: \n");
    printDataSet(dim, 1, query);
    getchar();

//    //start lsh_probing
//    result = lshProbing(dim, n_data, *L, *M, *W, hashTables, buckets, query, data);

    result = LSH_search(dim, *L, *M, *W, hashTables, buckets, query);

    printf("Result: \n");
    printDataSet(dim, 1, result);

//    generatePerturbationVectors(dim, *M, *W,
//                                5, query, hashTables[0]);

////    printf("Distance of query to data points in data set: \n");
//
//    int closestIdx = 0;
//    double closestDistance = RAND_MAX;
//
//    for (int i = 0; i < n_data; ++i) {
//        double *ele = getElementAtIndex(i, dim, n_data, data);
//        double distance = distanceOfTwoPoints(dim, query, ele);
//        if (distance < closestDistance) {
//            closestDistance = distance;
//            closestIdx = i;
//        }
////        printf("data %d: %f \n", i, distance);
//        free(ele);
//    }
//
//    if (result == NULL) {
//        printf("result = NULL");
//    } else {
//        printf("Closest data point: \n");
//        printDataSet(dim, 1, result);
//    }
//
//    printf("Closest idx: %d - distance: %f \n", closestIdx, closestDistance);
//
//    //verify variables
////    printHashTables(dim, *L, *M, hashTables);
////    printDataSet(dim, n_data, data);
//
//
//free pointer variables
    HashBucket *ite = buckets;
    while (ite != NULL) {
        HashBucket *temp = ite;
        ite = ite->next;
        for (
                int i = 0;
                i < *
                        L;
                ++i) {
            free(temp
                         ->hashValues[i]);
        }
        free(temp
                     ->hashValues);
        LinkedList *listIte = temp->head;

        while (listIte != NULL) {
            LinkedList *tempListIte = listIte;
            listIte = listIte->next;
            free(tempListIte
                         ->data);
            free(tempListIte);
        }

        free(temp);
    }

    for (
            int i = 0;
            i < *
                    L;
            ++i) {
        for (
                int j = 0;
                j < *
                        M;
                ++j) {
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