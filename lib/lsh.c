//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include <values.h>
#include "lsh.h"
#include "utils.h"
#include "data_structure.h"

HashBucket *hashBuckets = NULL;

int insert(int dim, int l, int m, double w, double ***hashTables, double *ele, double *centroid) {
    int **hashValues = calculateHashValues(dim, l, m, w ,centroid, hashTables, ele);

    if (hashBuckets == NULL) {
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

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets, double *centroid) {
    hashBuckets = buckets;
    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        insert(dim, l, m, w, hashTables, ele, centroid);
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

//??correct formula
double calculateScore(const int *a0, int length, struct pairZ zs[]) {
    double score = 0;
    for (int i = 0; i < length; ++i) {
//        score += zs[a0[i]].x;
        score += zs[a0[i]].x * zs[a0[i]].x;
    }
    return score;
}

struct HeapEle *minHeap(struct HeapEle **heap) {
    struct HeapEle *ite = heap[0];
    struct HeapEle *res = (struct HeapEle *) malloc(sizeof(struct HeapEle));
    struct HeapEle *min = NULL;
    double minScore = MAXDOUBLE;
    while (ite != NULL) {
        if (minScore > ite->score) {
            minScore = ite->score;
            res->length = ite->length;
            res->next = NULL;
            res->prev = NULL;
            res->data = (int *) malloc(res->length * sizeof(int));
            for (int i = 0; i < res->length; ++i) {
                res->data[i] = ite->data[i];
            }
            res->score = ite->score;
            min = ite;
        }
        ite = ite->next;
    }

    if (min->next != NULL && min->prev != NULL) {
        min->prev->next = min->next;
        min->next->prev = min->prev;
        min->next = NULL;
        min->prev = NULL;
    } else if (min->prev == NULL) {
        if (heap[0]->next == NULL) {
            free(heap[0]->data);
            *heap = NULL;
        }
    } else {
        min->prev->next = NULL;
        min->prev = NULL;
    }

    return res;
}

void insertEleHeap(struct HeapEle **heap, struct HeapEle *ele) {
    if (*heap == NULL) {
        *heap = ele;
    } else {
        ele->next = *heap;
        heap[0]->prev = ele;
        *heap = ele;
    }
}

bool isValid(struct HeapEle *ele, int twoM) {
    for (int i = 0; i < ele->length; ++i) {
        if (ele->data[i] > twoM)
            return false;

        int negJ = twoM - 1 - ele->data[i];
        for (int j = 0; j < ele->length; ++j) {
            if (ele->data[j] == negJ)
                return false;
        }
    }
    return true;
}


struct HeapEle *shiftHeap(struct HeapEle *ele, struct pairZ *zs) {
    int *data = (int *) malloc(ele->length * sizeof(int));

    for (int i = 0; i < ele->length - 1; ++i) {
        data[i] = ele->data[i];
    }

    data[ele->length - 1] = ele->data[ele->length - 1] + 1;

    struct HeapEle *shifted = (struct HeapEle *) malloc(sizeof(struct HeapEle));
    shifted->data = data;
    shifted->score = calculateScore(data, ele->length, zs);
    shifted->length = ele->length;
    shifted->next = NULL;
    shifted->prev = NULL;

    return shifted;
}

struct HeapEle *expandHeap(struct HeapEle *ele, struct pairZ *zs) {
    struct HeapEle *expanded = (struct HeapEle *) malloc(sizeof(struct HeapEle));
    expanded->data = (int *) malloc(ele->length + 1 * sizeof(int));
    for (int i = 0; i < ele->length; ++i) {
        expanded->data[i] = ele->data[i];
    }
    expanded->data[ele->length] = ele->data[ele->length - 1] + 1;
    expanded->length = ele->length + 1;
    expanded->score = calculateScore(expanded->data, expanded->length, zs);
    expanded->next = NULL;
    expanded->prev = NULL;
    return expanded;
}

int **probing(int numOfVectors, int dim, int l, int m, double w,
                double *query, double **hashFuncs) {
    int **perturbationVectors = (int **) malloc(numOfVectors * sizeof(int *));
    //find 2M array
    struct pairZ twoM[2 * m];

    int counter = 0;
    for (int i = 0; i < m; ++i) {
        struct pairZ z1 = {.x = distanceToBoundary(dim, w, query, hashFuncs[i], -1), .i = i, .r = -1};
        twoM[counter++] = z1;
        struct pairZ z2 = {.x = distanceToBoundary(dim, w, query, hashFuncs[i], 1), .i = i, .r = 1};
        twoM[counter++] = z2;
    }

    for (int i = 0; i < 2 * m; ++i) {
        printf("%f - %d \n", twoM[i].x, twoM[i].i);
    }

    for (int i = 0; i < 2 * m; ++i) {
        for (int j = i + 1; j < 2 * m - 1; ++j) {
            struct pairZ temp = twoM[i];
            if (twoM[j].x < temp.x) {
                twoM[i] = twoM[j];
                twoM[j] = temp;
            }
        }
    }

    for (int i = 0; i < 2 * m; ++i) {
        printf("%f - %d \n", twoM[i].x, twoM[i].i);
    }

    //generate Heap
    int *a0 = (int *) malloc(sizeof(int));
    a0[0] = 0;
    struct HeapEle *heap = &(struct HeapEle) {.data = a0, .length= 1, .score = calculateScore(a0, 1,
                                                                                              twoM), .next = NULL, .prev = NULL};

    //generate perturbation vectors
    for (int i = 0; i < numOfVectors; ++i) {
        struct HeapEle *minA;
        do {
            minA = minHeap(&heap);

            struct HeapEle *shifted = shiftHeap(minA, twoM);
            insertEleHeap(&heap, shifted);
            struct HeapEle *expanded = expandHeap(minA, twoM);
            insertEleHeap(&heap, expanded);
        } while (!isValid(minA, 2 * m));
        perturbationVectors[i] = (int *) malloc(m * sizeof(int));
        for (int j = 0; j < minA->length; ++j) {
            perturbationVectors[i][twoM[minA->data[j]].i] = twoM[minA->data[j]].r;
        }
        free(minA->data);
        free(minA);
    }

    printf("perturbation vectors: \n");
    for (int i = 0; i < numOfVectors; ++i) {
        printf("%d -- \n", i);
        for (int j = 0; j < m; ++j) {
            printf("%d \n", perturbationVectors[i][j]);
        }
    }

    while (heap->next != NULL) {
        if (heap->prev != NULL)
            free(heap->prev->data);
        free(heap->prev);
        heap = heap->next;
    }
    free(heap->data);
    free(heap);

//    getchar();

    return perturbationVectors;
}

double *LSH_search(int dim, int l, int m, double w, double ***hashTables,
                    HashBucket *buckets, double *query,
                   int **perturVectors, int num_vectors, double *centroid) {
    double distance = MAXDOUBLE;
    double *result = (double *) malloc(dim * sizeof(double));

    int **hashVal = calculateHashValues(dim, l, m, w, centroid ,hashTables, query);

    //convert perturbation vectors
    int ***probingHashVals = (int ***) malloc(num_vectors * sizeof(int **));
    for (int i = 0; i < num_vectors; ++i) {
        printf("pertur vector %d \n", i);
        probingHashVals[i] = (int **) malloc(l * sizeof(int *));
        for (int j = 0; j < l; ++j) {
            probingHashVals[i][j] = (int *) malloc(m * sizeof(int));
            for (int k = 0; k < m; ++k) {
                probingHashVals[i][j][k] = perturVectors[i][k] + hashVal[j][k];
                printf("%d \n", probingHashVals[i][j][k]);
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
    printf("Closest distance: %f",distance);

    return result;
}

double *LSH_probing(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query, int NUM_VECTORS, double *centroid) {
    double *result;
    int **hashVal = calculateHashValues(dim, l, m, w, centroid, hashTables, query);
    double distance = MAXDOUBLE;

    int **perturVectors = probing(NUM_VECTORS, dim, l, m, w, query, hashTables[0]);

    result = LSH_search(dim, l, m, w, hashTables, buckets, query, perturVectors, NUM_VECTORS, centroid);

    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }

//    printf("Closest distance: %f", distanceOfTwoPoints(dim, result, query));
    free(hashVal);
    getchar();
    return result;
}
