//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include <values.h>
#include "lsh.h"
#include "utils.h"
#include "data_structure.h"

HashBucket *hashBuckets = NULL;

int insert(int dim, int l, int m, double w, double ***hashTables, double *ele, double *centroid, int *num_hash_buckets) {
    int **hashValues = calculateHashValues(dim, l, m, w, centroid, hashTables, ele);

    if (hashBuckets == NULL) {
        *num_hash_buckets = *num_hash_buckets + 1;
        hashBuckets = (HashBucket *) malloc(sizeof(HashBucket));
        hashBuckets->hashValues = hashValues;
        LinkedList *head = (LinkedList *) malloc(sizeof(LinkedList));
        head->data = ele;
        head->next = NULL;
        hashBuckets->head = head;
        hashBuckets->next = NULL;
        return 0;
    }

    HashBucket *ite = hashBuckets;

    while (ite != NULL) {
        if (compareHashValues(l, m, hashValues, ite->hashValues)) {
            LinkedList *node = (LinkedList *) malloc(sizeof(LinkedList));
            node->data = ele;
            node->next = ite->head;
            ite->head = node;
            return 0;
        }
        ite = ite->next;
    }

    *num_hash_buckets = *num_hash_buckets + 1;

    HashBucket *bucket = (HashBucket *) malloc(sizeof(HashBucket));
    bucket->hashValues = hashValues;
    LinkedList *head = (LinkedList *) malloc(sizeof(LinkedList));
    head->data = ele;
    head->next = NULL;
    bucket->head = head;
    bucket->next = hashBuckets;

    hashBuckets = bucket;

    return 0;
}

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets,
                double *centroid, int* num_hash_buckets) {
    *num_hash_buckets = 0;
    hashBuckets = buckets;
    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        insert(dim, l, m, w, hashTables, ele, centroid, num_hash_buckets);
    }

    return hashBuckets;
}

double search(int dim, HashBucket *bucket, double *query, double minDistance, double *result_ptr) {
    double distance = MAXDOUBLE;
    double *result = (double *) malloc(dim * sizeof(double));

    LinkedList *data = bucket->head;
    while (data->next) {
        double currentDistance = distanceOfTwoPoints(dim, data->data, query);
        if (distance > currentDistance) {
            distance = currentDistance;
            for (int i = 0; i < dim; ++i) {
                result[i] = data->data[i];
            }
        }
        data = data->next;
    }

    if (distance < minDistance) {
        for (int i = 0; i < dim; ++i) {
            result_ptr[i] = result[i];
        }
        return distance;
    }

    return minDistance;
}


int **probing(int numOfVectors, int dim, int l, int m, double w,
              double *query, double **hashFuncs) {
    int **perturbationVectors = (int **) malloc(numOfVectors * sizeof(int *));

    printf("perturbation vectors: \n");
    for (int i = 0; i < numOfVectors; ++i) {
        printf("%d -- \n", i);
        for (int j = 0; j < m; ++j) {
            printf("%d \n", perturbationVectors[i][j]);
        }
    }

    getchar();
    return perturbationVectors;
}

double *LSH_search(int dim, int l, int m, double w, double ***hashTables,
                   HashBucket *buckets, double *query,
                   int **perturVectors, int num_vectors, double *centroid, double *distanceB4Probing) {
    double distance = MAXDOUBLE;
    double *result = (double *) malloc(dim * sizeof(double));

    int **hashVal = calculateHashValues(dim, l, m, w, centroid, hashTables, query);

    //convert perturbation vectors
    int ***probingHashVals = (int ***) malloc(num_vectors * sizeof(int **));
    for (int i = 0; i < num_vectors; ++i) {
//        printf("pertur vector %d \n", i);
        probingHashVals[i] = (int **) malloc(l * sizeof(int *));
        for (int j = 0; j < l; ++j) {
            probingHashVals[i][j] = (int *) malloc(m * sizeof(int));
            for (int k = 0; k < m; ++k) {
                probingHashVals[i][j][k] = perturVectors[i][k] + hashVal[j][k];
//                printf("%d \n", probingHashVals[i][j][k]);
            }
        }
    }

//    getchar();

    HashBucket *ite = buckets;

    while (ite != NULL) {
        if (compareHashValues(l, m, hashVal, ite->hashValues)) {
            for (int i = 0; i < l; ++i) {
                free(hashVal[i]);
            }
            free(hashVal);
            distance = search(dim, ite, query, MAXDOUBLE, result);
            break;
        }
        ite = ite->next;
    }

    *distanceB4Probing = distance;

    printf("Closest distance b4 probing: %f \n", distance);

    for (int i = 0; i < num_vectors; ++i) {
        ite = buckets;
        while (ite != NULL) {
            if (compareHashValues(l, m, probingHashVals[i], ite->hashValues)) {
//                printf("search \n");
//                getchar();
                for (int j = 0; j < l; ++j) {
                    free(probingHashVals[i][j]);
                }
                free(probingHashVals[i]);
                distance = search(dim, ite, query, distance, result);
                break;
            }
            ite = ite->next;
        }
    }
    printf("Closest distance: %f", distance);

    printDataSet(dim, 1, result);

//    getchar();

    return result;
}

double *
LSH_probing(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query, int NUM_VECTORS,
            double *centroid, double *distanceB4Probing) {
    double *result;
    int **hashVal = calculateHashValues(dim, l, m, w, centroid, hashTables, query);
    double distance = MAXDOUBLE;

    int **perturVectors = probing(NUM_VECTORS, dim, l, m, w, query, hashTables[0]);

    result = LSH_search(dim, l, m, w, hashTables, buckets, query, perturVectors, NUM_VECTORS, centroid,
                        distanceB4Probing);

    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }

//    printf("Closest distance: %f", distanceOfTwoPoints(dim, result, query));
    free(hashVal);
//    getchar();
    return result;
}
