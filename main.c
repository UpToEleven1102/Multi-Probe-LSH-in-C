#include <stdio.h>
#include "lib/utils.h"

void initParameters(int *L, int *M, float *W) {

}

int main() {
    const int dim = 4;
    const int n_data = 20;
    const double *data = generateDataSet(dim, n_data);

    int *L, *M;
    float *W;

    initParameters(L, M, W);
    double ***hashTables = generateHashTables(*L, *M, dim);

//    a(L, M, W, )

    printDataSet(dim,n_data, data);
}