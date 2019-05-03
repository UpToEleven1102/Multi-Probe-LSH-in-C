//
// Created by huyen on 4/22/19.
//

#ifndef LSH_PROBING_LSH_SEARCH_H
#define LSH_PROBING_LSH_SEARCH_H

#include "utils.h"

struct BucketHashVal {
    int **value;
    double score;
    struct BucketHashVal *next;
    struct BucketHashVal *prev;
};

typedef struct BucketHashVal BucketHashVal;

int _LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, int num_hash_buckets,
                double *query, double *centroid, double *distanceB4Probing, double *result, int *num_checked_data_points);

#endif //LSH_PROBING_LSH_SEARCH_H
