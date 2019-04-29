//
// Created by huyen on 3/29/19.
//

#include <stdlib.h>
#include <stdio.h>
#include <values.h>
#include <math.h>
#include "lsh_main.h"
#include "utils.h"
#include "lsh.h"
#include "lsh_search.h"
#include <time.h>

double z1;

//what mean and deviation that want?

//phase = 0 or 1 equivalent to z0 or z1
//Box-Mulller transformation
double gaussian_rand(double mean, double stdDev, int phase) {
    const double epsilon = 2.22507e-308;
    const double two_pi = 2 * M_PI;

    if (phase == 1)
        return z1;

    double u1, u2;
    do {
        u1 = (rand() + 1.) / (RAND_MAX + 2.);
        u2 = rand() / (RAND_MAX + 1.);
    } while (u1 <= epsilon);

    double z0;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

    return z0;
}

double *newUnitVector(int dim, double *mean, double *stdDev) {
//    double *unitVector = generateDataSet(dim, 1);
//

    double *unitVector = (double *) malloc(dim * sizeof(double));

    int phase = 0;

    for (int i = 0; i < dim; ++i) {
        unitVector[i] = gaussian_rand(mean[i], stdDev[i], phase);
        phase = 1 - phase;
    }

    double vectorLength = 0;
    for (int i = 0; i < dim; ++i) {
        vectorLength += unitVector[i] * unitVector[i];
    }

    vectorLength = sqrt(vectorLength);

    for (int i = 0; i < dim; ++i) {
        unitVector[i] = unitVector[i] / vectorLength;
    }

    return unitVector;
}

double ***generateHashTables(int l, int m, int dim, double *mean, double *stdDev) {
    double ***hashTables = (double ***) malloc(l * sizeof(double **));

    for (int i = 0; i < l; ++i) {
        hashTables[i] = (double **) malloc(m * sizeof(double *));
        for (int j = 0; j < m; ++j) {
            hashTables[i][j] = newUnitVector(dim, mean, stdDev);
        }
    }

    return hashTables;
}

//probe paremeters, write a function
void initParameters(int *l_ptr, int *m_ptr, double *w_ptr, double *mean, double *stdDev, int dim, int n_data,
                    const double *data,
                    double *centroid, double *dataSpread) {
    //comeback and pick this up later
//    *m_ptr = (int) floor(dim / 4.0);

    *m_ptr = 4;

//    *l_ptr = *m_ptr * 3;

    *l_ptr = 3;

    double **buff = (double **) malloc(dim * sizeof(double *));

    for (int i = 0; i < dim; ++i) {
        buff[i] = (double *) malloc(3 * sizeof(double));
        buff[i][0] = -RAND_MAX;
        buff[i][1] = RAND_MAX;
        buff[i][2] = 0;
        stdDev[i] = 0;
    }

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);

        for (int j = 0; j < dim; ++j) {
            centroid[j] += ele[j];

            if (ele[j] > buff[j][0]) {
                buff[j][0] = ele[j];
            }
            if (ele[j] < buff[j][1]) {
                buff[j][1] = ele[j];
            }
            buff[j][2] += ele[j];
        }
        free(ele);
    }

    double distance = 0;
    for (int i = 0; i < dim; ++i) {
        printf("%f, %f, %f \n", buff[i][0], buff[i][1], buff[i][2] );
        centroid[i] /= n_data;
        mean[i] = buff[i][2] / n_data;
        distance += (buff[i][0] - buff[i][1]) * (buff[i][0] - buff[i][1]);

    }

//    getchar();

    distance = sqrt(distance / dim);

    printf("distance: %f \n", distance);

//    getchar();

    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);

        for (int j = 0; j < dim; ++j) {
            stdDev[j] = (ele[j] - mean[j]) * (ele[j] - mean[j]);
        }

        free(ele);
    }

    double _mean = 0;

    for (int i = 0; i < dim; ++i) {
        stdDev[i] = sqrt(stdDev[i] / n_data);
        _mean += mean[i];
    }

    _mean = _mean / dim;

    //stop changing this line, you're dumb-ass
    *dataSpread = distance;


//    *w_ptr = _mean;

    *w_ptr = .6 * distance;

    for (int i = 0; i < dim; ++i) {
        free(buff[i]);
    }
    free(buff);
}

int LSH_main(int dim, int n_data, double *data, HashBucket *buckets, int num_queries, double **queries,
             FILE *file) {
    clock_t start, end;
    double generateBucketsTime, searchTime;
    double *distanceB4Probing = (double *) malloc(sizeof(double));
    double *dataSpread = (double *) malloc(sizeof(double));

    int *l_ptr = (int *) malloc(sizeof(int));
    int *m_ptr = (int *) malloc(sizeof(int));
    double *w_ptr = (double *) malloc(sizeof(double));

    double *result = (double *) calloc(dim, sizeof(double));

    double *centroid = (double *) calloc(dim, sizeof(double));

    double *mean = (double *) malloc(dim * sizeof(double));
    double *stdDev = (double *) malloc(dim * sizeof(double));

    initParameters(l_ptr, m_ptr, w_ptr, mean, stdDev, dim, n_data, data, centroid, dataSpread);
//    printf("l_ptr - %d, m_ptr - %d, w_ptr - %f, dim - %d \n", *l_ptr, *m_ptr, *w_ptr, dim);

    double ***hashTables = generateHashTables(*l_ptr, *m_ptr, dim, mean, stdDev);
//    printHashTables(dim, *l_ptr, *m_ptr, hashTables);
    start = clock();

    int *num_hash_buckets = (int *) malloc(sizeof(int));

    buckets = LSH(dim, n_data, *l_ptr, *m_ptr, *w_ptr, hashTables, data, NULL, centroid, num_hash_buckets);

    end = clock();
    generateBucketsTime = end - start;

//    printf("hash buckets: \n");
//
//    int numBuckets = printHashBuckets(dim, *l_ptr, *m_ptr, buckets);


    printf("number of buckets: %d \n Enter to continue \n", *num_hash_buckets);

//    for (int i = 0; i < num_queries; ++i) {
        start = clock();
        int *num_checked_data_points = (int *) calloc(1, sizeof(int));

        int num_check_after_probing = _LSH_search(dim, *l_ptr, *m_ptr, *w_ptr, hashTables, buckets, *num_hash_buckets,
                                                  queries[0], centroid,
                                                  distanceB4Probing, result, num_checked_data_points);
        end = clock();

        searchTime = end - start;


        fprintf(file,
                "n_data: %d, dim: %d, data spread: %f, w: %f, l: %d, m: %d, num bucket: %d, generate buckets time: %f, search time: %f, distance b4 probing: %f, after probing: %f, num data checked b4 probing: %d, after probing: %d \n",
                n_data, dim, *dataSpread, *w_ptr, *l_ptr, *m_ptr, *num_hash_buckets, generateBucketsTime, searchTime,
                *distanceB4Probing,
                distanceOfTwoPoints(dim, result, queries[0]), *num_checked_data_points, num_check_after_probing);

//    getchar();

        //verify distance
        int closestIdx = 0;
        double closestDistance = MAXDOUBLE;

        for (int i = 0; i < n_data; ++i) {
            double *ele = getElementAtIndex(i, dim, n_data, data);
            double distance = distanceOfTwoPoints(dim, queries[0], ele);
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
            closestDistance = distanceOfTwoPoints(dim, queries[0], result);
        }

        printf("Closest idx: %d - distance: %f \n", closestIdx, closestDistance);
        fprintf(file, " exact closest distance: %f \n", closestDistance);
//    }

//verify distance

    HashBucket *ite = buckets;
    while (ite != NULL) {
        HashBucket *temp = ite;
        ite = ite->next;
        for (
                int i = 0;
                i < *
                        l_ptr;
                ++i) {
            free(temp->hashValues[i]);
        }
        free(temp->hashValues);
        LinkedList *listIte = temp->head;

        while (listIte != NULL) {
            LinkedList *tempListIte = listIte;
            listIte = listIte->next;
            free(tempListIte->data);
            free(tempListIte);
        }

        free(temp);
    }

    for (
            int i = 0;
            i < *l_ptr;
            ++i) {
        for (int j = 0; j < *m_ptr; ++j) {
            free(hashTables[i][j]);
        }
        free(hashTables[i]);
    }
    free(hashTables);
    free(l_ptr);
    free(m_ptr);
    free(w_ptr);
}
