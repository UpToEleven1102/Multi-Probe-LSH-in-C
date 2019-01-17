//
// Created by huyen on 1/15/19.
//

#ifndef LSH_PROBING_UTILS_H
#define LSH_PROBING_UTILS_H

double* generateDataSet(int dim, int n_data);
void printDataSet(int dim, int n_data, double *data);
double ***generateHashTables(int l, int m, int dim);

#endif //LSH_PROBING_UTILS_H
