#include <stdio.h>
#include "lib/utils.h"

int main() {
    const int dim = 4;
    const int n_data = 20;

    float *data = generateDataSet(dim, n_data);

    printDataSet(dim,n_data, data);
}