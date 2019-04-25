//
// Created by huyen on 1/15/19.
//

#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"


double distanceOfTwoPoints(int dim, const double *point1, const double *point2) {
    double distance = 0;

    for (int i = 0; i < dim; ++i) {
        distance += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }

    return sqrt(distance);
}

double innerProduct(const double *vector1, const double *vector2, int dim) {
    double product = 0;
    for (int i = 0; i < dim; ++i) {
        product += vector1[i] * vector2[i];
    }
    return product;
}

//bool isEqualArrays(int dim, const double *arr1, const double *arr2) {
//    for (int i = 0; i < dim; ++i) {
//        if (arr1[i] != arr2[i])
//            return false;
//    }
//
//    return true;
//}

bool compareHashValues(int l, int m, int **hashValue1, int **hashValue2) {
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < m; ++j) {
            if (hashValue1[i][j] != hashValue2[i][j])
                return false;
        }
    }

    return true;
}


double *getElementAtIndex(int idx, int dim, int n_data, const double *data) {
    double *ele = (double *) malloc(dim * sizeof(double));

    if (idx < n_data)
        for (int i = 0; i < dim; ++i) {
            ele[i] = data[idx * dim + i];
        }
    return ele;
}

void printDataSet(int dim, int n_data, const double *data) {
    int counter = 0;
    for (int i = 0; i < dim * n_data; ++i) {
        if (i % dim == 0) {
            printf("%d ---\n", counter++);
        }
        printf("%f \n", data[i]);
    }
}

double distanceToBoundary(int dim, double w, double *query, double *hashFunc, double *centroid, int r) {
    if (r == 0) return 0;

    double b = innerProduct(hashFunc, centroid, dim);

    double distanceToTheLeft = (innerProduct(query, hashFunc, dim) + b) - calculateHashValue(dim, w, query, hashFunc, centroid) * w;

    if (r < 0)
        return distanceToTheLeft;
    if (r > 0)
        return w - distanceToTheLeft;
}


int calculateHashValue(int dim, double w, double *ele, double *hashFunc, double *centroid) {
    double b = innerProduct(hashFunc, centroid, dim);

    return (int) floor((innerProduct(ele, hashFunc, dim) - b)/ w);
}

double calculateDistanceToBucket(int dim, int l, int m, double w, int **hashVal, int **bucketHashVal, double *query, double ***hashTables, double *centroid) {
    double distance = 0;

    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < m; ++j) {
            int steps = (hashVal[i][j] - bucketHashVal[i][j]);
            if (abs(steps) > 0){
                distance += (abs(steps) - 1) * w;
            }

            distance = distance + distanceToBoundary(dim, w, query, hashTables[i][j], centroid, steps); // + distance to the boundary
        }
    }

    return distance;
}


//b = centroid * hash functions??????? data dependent
int calculateHashValues(int dim, int l, int m, double w, double *centroid, double ***hashTables, double *ele, int **hashValues) {
    double b;

    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < m; ++j) {
            b = innerProduct(hashTables[i][j], centroid, dim);
            hashValues[i][j] = (int) floor((innerProduct(ele, hashTables[i][j], dim) - b) / w);
//            hashValues[i][j] = (int) ((innerProduct(ele, hashTables[i][j], dim)) / w);
        }
    }
    return 0;
}

void printHashTables(int dim, int l, int m, double ***tables) {
    for (int i = 0; i < l; ++i) {
        printf("table %d \n", i);
        for (int j = 0; j < m; ++j) {
            printf("function # %d \n", j);
            printDataSet(dim, 1, tables[i][j]);
        }
    }
}

int printHashBuckets(int dim, int l, int m, HashBucket *buckets) {
    HashBucket *ite = buckets;
    int counter = 0;

    printf("-- parameters: dim %d l %d m %d \n", dim, l, m);

    while (ite != NULL) {
        printf("Bucket %d \n Hash Values: -- \n", counter++);

        for (int i = 0; i < l; ++i) {
            printf("Table %d -- \n", i);
            for (int j = 0; j < m; ++j) {
                printf("%d \n", ite->hashValues[i][j]);
            }
        }

        printf("Elements: \n");

        LinkedList *listIte = ite->head;
        int counter2 = 0;
        while (listIte != NULL) {
            printf("%d ----------\n", counter2++);
            printDataSet(dim, 1, listIte->data);
            listIte = listIte->next;
        }

        free(listIte);
        ite = ite->next;
    }

    free(ite);

    return counter;
}

int printHashValues(int l, int m, int **hashValue) {
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < m; ++j) {
            printf("%d \n", hashValue[i][j]);
        }
    }
}
