//
// Created by huyen on 1/17/19.
//

#ifndef LSH_PROBING_LSH_H
#define LSH_PROBING_LSH_H

#include "utils.h"

struct HeapEle {
    int *data;
    double score;
    int length;
    struct HeapEle *next;
    struct HeapEle *prev;
};

struct HeapTreeNode {
    int *data;
    double score;
    int length;
    struct HeapTreeNode *left;
    struct HeapTreeNode *right;
    struct HeapTreeNode *parent;
};

//double **calculateHashValues(int dim, int l, int m, double w, double ***hashTables, double *ele);

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets, double *centroid, int *num_hash_buckets);

#endif //LSH_PROBING_LSH_H
