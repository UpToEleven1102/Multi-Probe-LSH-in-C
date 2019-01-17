//
// Created by huyen on 1/15/19.
//

#ifndef LSH_PROBING_UTILS_H
#define LSH_PROBING_UTILS_H

double* generateDataSet(int dim, int n_data);
void printDataSet(int dim, int n_data, const double *data);
double ***generateHashTables(int l, int m, int dim);

double *getElementAtIndex(int idx, int dim, int n_data, const double *data);

void printHashTables(int dim, int l, int m, double ***tables);

#endif //LSH_PROBING_UTILS_H
