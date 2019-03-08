//
// Created by huyen on 1/21/19.
//

#ifndef LSH_PROBING_LSH_PROBING_H
#define LSH_PROBING_LSH_PROBING_H

#include "utils.h"

int **generatePerturbationVectors(int dim, int m, double w, int t, double *query, double **hashTable);

double *lshProbing(int dim, int n_data, int l, int m, double w, double*** hashTables, HashBucket *buckets, double *query, double *data);

#endif //LSH_PROBING_LSH_PROBING_H
