//
// Created by huyen on 1/17/19.
//

#ifndef LSH_PROBING_LSH_H
#define LSH_PROBING_LSH_H

#include "utils.h"

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data);

#endif //LSH_PROBING_LSH_H
