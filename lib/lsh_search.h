//
// Created by huyen on 4/22/19.
//

#ifndef LSH_PROBING_LSH_SEARCH_H
#define LSH_PROBING_LSH_SEARCH_H

#include "utils.h"

int _LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query, int NUM_VECTORS,
            double *centroid, double *distanceB4Probing, double *result);
#endif //LSH_PROBING_LSH_SEARCH_H
