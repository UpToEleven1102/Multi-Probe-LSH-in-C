//
// Created by huyen on 1/17/19.
//

#ifndef LSH_PROBING_LSH_H
#define LSH_PROBING_LSH_H

#include "utils.h"

//double **calculateHashValues(int dim, int l, int m, double w, double ***hashTables, double *ele);

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets);

#endif //LSH_PROBING_LSH_H
