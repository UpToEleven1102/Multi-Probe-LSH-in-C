#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <values.h>
#include "lib/utils.h"
#include "lib/lsh.h"
#include "lib/lsh_main.h"

//TODO: more research about number L and M, how many are needed
//TODO: implement b
// how large we should choose value of L, paper
// hash functions are unit vectors?
// hashValue formula   <a.v -b> / w. what is b??

//b = centroid * hash functions??????? data dependent

//W determine how wide the slot is - hi(q) = hash function,  fi(q) = hi(q) * w (without getting floor)

//W is dependent to the number of data points

//slides: intuitions wrong: 1-step buckets are not better than 2-step buckets, paper: 1step buckets are not neccessary better than buckets that are 2 steps away

//note: get the lowest score of the heap or get the top element

//construct perturbation vectors in each hashtable??

//know when and what table to apply the perturbation vector to

//what if there is no collision


double **readCSVFile(int n_data, int dim, int num_data_sets, double *query) {
    double **dataSets = (double **) malloc(num_data_sets * sizeof(double *));
    for (int i = 0; i < num_data_sets; ++i) {
        dataSets[i] = (double *) malloc(dim * n_data * sizeof(double));
    }

    FILE *file = fopen("../data_sets/HIGGS.csv", "rb");

    char line[1024];
    int counter;

    for (int j = 0; j < num_data_sets; ++j) {
        counter = 0;
        for (int i = 0; (fscanf(file, "%s", line) == 1); ++i) {
            const char *tok;

            for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",")) {
                dataSets[j][counter] = strtof(tok, NULL);
                counter++;
            }

            if (i == n_data)
                break;
        }
    }

    //get a data point as a query
    counter = 0;
    fscanf(file, "%s", line);
    const char *tok;

    for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",")) {
        query[counter] = strtof(tok, NULL);
        counter++;
    }

    fclose(file);

    return dataSets;
}

double **readBinaryFile(int n_data, int dim, int num_data_sets, double *query) {
    double **dataSets = (double **) malloc(num_data_sets * sizeof(double *));
    for (int i = 0; i < num_data_sets; ++i) {
        dataSets[i] = (double *) calloc(dim * n_data, sizeof(double));
    }

    FILE *file = fopen("../data_sets/tr_HIGGS.dat", "rb");

    char line[1024];
    int counter = 0;

    //data 1
    if (file != NULL) {
        for (int i = 0; i < num_data_sets; ++i) {
            fread(dataSets[i], sizeof(double), dim * n_data, file);
            if (feof(file))
                break;
        }

        fread(query, sizeof(double), dim, file);
        fclose(file);
    } else {
        printf("Failed to open file. ");
    }

    return dataSets;
}

int main() {
    srand(1);
    const int dim = 29;
    const int n_data = 1000000;
    const int NUM_DATA_SETS = 3;

    double *query, *result, *centroid, ***hashTables;
    query = (double *) malloc(dim * sizeof(double));
    result = (double *) malloc(dim * sizeof(double));

    double **dataSets = readBinaryFile(n_data, dim, NUM_DATA_SETS, query);
    double *data = dataSets[1];

//    printDataSet(dim, n_data, data);
    HashBucket *buckets = malloc(sizeof(HashBucket));

    FILE *oF = fopen("output.txt", "a");


    LSH_main(dim, n_data, data, hashTables, buckets, centroid, result, query, oF);

//verify distance
    int closestIdx = 0;
    double closestDistance = MAXDOUBLE;

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        double distance = distanceOfTwoPoints(dim, query, ele);
        if (distance < closestDistance) {
            closestDistance = distance;
            closestIdx = i;
            for (int j = 0; j < dim; ++j) {
                result[j] = ele[j];
            }
        }
//        printf("data %d: %f \n", i, distance);
        free(ele);
    }

    if (result == NULL) {
        printf("result = NULL");
    } else {
        printf("Closest data point: \n");
//        printDataSet(dim, 1, result);
        closestDistance = distanceOfTwoPoints(dim, query, result);
    }

    printf("Closest idx: %d - distance: %f \n", closestIdx, closestDistance);
    fprintf(oF, " exact closest distance: %f \n", closestDistance);
    fclose(oF);

    //verify variables
//    printHashTables(dim, *L, *M, hashTables);
//    printDataSet(dim, n_data, data);


//free pointer variables
    for (int i = 0; i < NUM_DATA_SETS; ++i) {
        free(dataSets[i]);
    }
    free(dataSets);
}