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
    double **hashValues;
    struct HashBucket *next;
};

typedef struct LinkedList LinkedList;
typedef struct HashBucket HashBucket;


double *calculateCentroid(int dim, int n_data, double *data);

double* generateDataSet(int dim, int n_data);
void printDataSet(int dim, int n_data, const double *data);
double ***generateHashTables(int l, int m, int dim);
bool compareHashValues(int l, int m, double **hashValue1, double **hashValue2);
double distanceOfTwoPoints(int dim, double *point1, double *point2);

double *getElementAtIndex(int idx, int dim, int n_data, const double *data);

double calculateHashValue(int dim, double w, double *ele, const double *centroid, double *hashFunc);

double **calculateHashValues(int dim, int l, int m, double w, const double *centroid, double ***hashTables, double *ele);

void printHashTables(int dim, int l, int m, double ***tables);

void printHashBuckets(int dim, int l, int m, HashBucket *buckets);

#endif //LSH_PROBING_UTILS_H
