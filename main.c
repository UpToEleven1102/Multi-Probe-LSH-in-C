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



//TODO why so many data points get into 1 bucket

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

int readBinaryFile(int n_data, int dim, int num_data_sets, double **dataSets, int num_queries, double **queries) {
    FILE *file = fopen("../data_sets/tlc_nyc2016_norm_41M_dim16.dat", "rb");
//    FILE *file = fopen("../data_sets/tr_HIGGS.dat", "rb");
//    FILE *file = fopen("../data_sets/heterogeneity_activity_norm.dat", "rb");

    char line[1024];
    int counter = 0;

    //data 1
    if (file != NULL) {
        for (int i = 0; i < num_data_sets; ++i) {
            fread(dataSets[i], sizeof(double), dim * n_data, file);
            if (feof(file))
                break;
        }

        for (int i = 0; i < num_queries; ++i) {
            fread(queries[i], sizeof(double), dim, file);
            if (feof(file))
                break;
        }
        fclose(file);
    } else {
        printf("Failed to open file. ");
    }

    return 0;
}

int main() {
    srand(1);
    const int dim = 16;
//    const int dim = 28;
    const int num_data_points = 1000000;
    const int num_queries = num_data_points / 100;
    const int n_data = num_data_points - num_queries;
    const int NUM_DATA_SETS = 1;

    double **dataSets = (double **) malloc(NUM_DATA_SETS * sizeof(double *));
    for (int i = 0; i < NUM_DATA_SETS; ++i) {
        dataSets[i] = (double *) calloc(dim * n_data, sizeof(double));
    }

    double **queries = (double **) malloc(num_queries * sizeof(double *));
    for (int i = 0; i < num_queries; ++i) {
        queries[i] = (double *) calloc(dim, sizeof(double));
    }

    //read data
    readBinaryFile(n_data, dim,
            NUM_DATA_SETS, dataSets, num_queries, queries);//outputs


    double *data = dataSets[0];

//    printDataSet(dim, n_data, data);
    HashBucket *buckets = malloc(sizeof(HashBucket));


    FILE *oF = fopen("tlc.txt", "a");


    LSH_main(dim, n_data, data, buckets, num_queries, queries, oF);


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