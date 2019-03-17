//
// Created by huyen on 1/15/19.
//

#ifndef LSH_PROBING_UTILS_H
#define LSH_PROBING_UTILS_H

typedef int bool;
#define true 1;
#define  false 0;

struct LinkedList {
    double *data;
    struct LinkedList *next;
};

struct HashBucket {
    struct LinkedList *head;
    int **hashValues;
    struct HashBucket *next;
};

typedef struct LinkedList LinkedList;
typedef struct HashBucket HashBucket;

void printDataSet(int dim, int n_data, const double *data);

bool compareHashValues(int l, int m, int **hashValue1, int **hashValue2);

//bool isEqualArrays(int dim, const double *arr1, const double *arr2);

double distanceOfTwoPoints(int dim, double *point1, double *point2);

double *getElementAtIndex(int idx, int dim, int n_data, const double *data);

double scorePerturbationVector(int dim, int m, double w, double *query, double **hashTable, int *vector);

double distanceToBoundary(int dim, double w, double *query, double *hashFunc, int r);

int calculateHashValue(int dim, double w, double *ele, double *hashFunc);

int **calculateHashValues(int dim, int l, int m, double w, double *centroid, double ***hashTables, double *ele);

void printHashTables(int dim, int l, int m, double ***tables);

int printHashBuckets(int dim, int l, int m, HashBucket *buckets);

int printHashValues(int l, int m, double **hashValue);

#endif //LSH_PROBING_UTILS_H
