//
// Created by huyen on 3/29/19.
//

#ifndef LSH_PROBING_LSH_MAIN_H
#define LSH_PROBING_LSH_MAIN_H

#include "utils.h"
#include <stdio.h>

int LSH_main(int dim, int n_data, double *data, HashBucket *buckets, int num_queries, double **queries, FILE *file);
#endif //LSH_PROBING_LSH_MAIN_H
