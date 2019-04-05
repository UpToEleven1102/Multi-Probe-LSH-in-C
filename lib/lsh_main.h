//
// Created by huyen on 3/29/19.
//

#ifndef LSH_PROBING_LSH_MAIN_H
#define LSH_PROBING_LSH_MAIN_H

#include "utils.h"

int LSH_main(int dim, int n_data, double *data,
             double ***hashTables, HashBucket *buckets, double *centroid, double *result,
             double *datum, FILE *file);
#endif //LSH_PROBING_LSH_MAIN_H
