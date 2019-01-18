//
// Created by huyen on 1/17/19.
//

#import "utils.h"

#ifndef LSH_PROBING_LSH_H
#define LSH_PROBING_LSH_H

double **calculateHashValues(int dim, int l, int m, double w,  double ***hashTables, double *ele);

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data);

#endif //LSH_PROBING_LSH_H
