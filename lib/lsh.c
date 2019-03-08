//
// Created by huyen on 1/17/19.
//

#include <malloc.h>
#include <values.h>
#include "lsh.h"
#include "utils.h"
#include "data_structure.h"

HashBucket *hashBuckets = NULL;

int insert(int dim, int l, int m, double w, double ***hashTables, double *ele) {
    int **hashValues = calculateHashValues(dim, l, m, w, hashTables, ele);

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

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets) {
    hashBuckets = buckets;
    for (int i = 0; i < n_data; ++i) {
        double *ele = getElementAtIndex(i, dim, n_data, data);
        insert(dim, l, m, w, hashTables, ele);
    }

    return hashBuckets;
}

double search(int dim, HashBucket *bucket, double *query, double *result) {
    double distance = MAXDOUBLE;

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
    printf("Closest distance: %f", distance);

    return distance;
}

//??correct formula
double calculateScore(const int *a0, int length, struct pairZ zs[]) {
    double score = 0;
    for (int i = 0; i < length; ++i) {
        score += zs[a0[i]].x;
//        score += zs[a0[i]].x * zs[a0[i]].x;
    }
    return score;
}

struct HeapEle *minHeap(struct HeapEle *heap) {
    struct HeapEle *ite = heap;
    struct HeapEle *min = (struct HeapEle*) malloc(sizeof(struct HeapEle));
    double minScore = MAXDOUBLE;
    while(ite != NULL) {
        if (minScore > ite->score) {
            minScore = ite->score;
            *min = *ite;
        }
        ite = ite->next;
    }

    //TODO: remove min element

    return min;
}

bool isValid(struct HeapEle *ele, int twoM) {
    for (int i = 0; i < ele->length; ++i) {
        if (ele->data[i] > twoM)
            return false;

        int negJ = twoM - 1 - ele->data[i];
        for (int j = 0; j < ele->length; ++j) {
            if(ele->data[j] == negJ)
                return false;
        }
    }
    return true;
}



struct HeapEle *shiftHeap(struct HeapEle *ele, struct pairZ *zs) {
    int *data = (int *)malloc(ele->length*sizeof(int));

    for (int i = 0; i < ele->length-1; ++i) {
        data[i] = ele->data[i];
    }

    data[ele->length - 1] = ele->data[ele->length - 1] +1;

    struct HeapEle *shifted = (struct HeapEle *)malloc(sizeof(struct HeapEle));
    shifted->data = data;
    shifted->score = calculateScore(data, ele->length, zs);
    shifted->length = ele->length;
    shifted->next = NULL;
    shifted->prev = NULL;

    return shifted;
}

struct HeapEle *expandHeap(struct HeapEle *ele, struct pairZ *zs) {
    struct HeapEle *expanded = (struct HeapEle*)malloc(sizeof(struct HeapEle));
    expanded->data = (int*)malloc(ele->length+1 *sizeof(int));
    for (int i = 0; i < ele->length; ++i) {
        expanded->data[i] = ele->data[i];
    }
    expanded->data[ele->length] = ele->data[ele->length-1]+1;
    expanded->length = ele->length +1;
    expanded->score = calculateScore(expanded->data, expanded->length, zs);
    expanded->next = NULL;
    expanded->prev = NULL;
    return expanded;
}
int **probing(int numOfVectors, int dim, int l, int m, double w, double *query, double **hashFuncs) {
    int **perturbationVectors = (int **)malloc(numOfVectors * sizeof(int*));
    //find 2M array
    struct pairZ twoM[2*m];

    int counter = 0;
    for (int i = 0; i < m; ++i) {
        struct pairZ z1 = { .x = distanceToBoundary(dim, w, query, hashFuncs[i], -1), .i = i, .r = -1};
        twoM[counter++] = z1;
        struct pairZ z2 = { .x = distanceToBoundary(dim, w, query, hashFuncs[i], 1), .i = i, .r = 1};
        twoM[counter++] = z2;
    }

    for (int i = 0; i < 2 * m; ++i) {
        printf("%f - %d \n", twoM[i].x, twoM[i].i);
    }

    for (int i = 0; i < 2 * m; ++i) {
        for (int j = i+1; j < 2 * m - 1; ++j) {
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
    int *a0 = (int*) malloc(sizeof(int));
    a0[0] = 0;
    struct HeapEle *heap = &(struct HeapEle){.data = a0, .length= 1, .score = calculateScore(a0, 1, twoM), .next = NULL, .prev = NULL};

    //generate perturbation vectors
    for (int i = 0; i < numOfVectors; ++i) {
        struct HeapEle *minA;
        do {
            minA = minHeap(heap);
            struct HeapEle *shifted = shiftHeap(minA, twoM);

            shifted->next = heap;
            heap->prev = shifted;
            heap = shifted;

            struct HeapEle *expanded = expandHeap(minA, twoM);

            expanded->next = heap;
            heap->prev = expanded;
            heap = expanded;
        } while (!isValid(minA, 2 * m));
        perturbationVectors[i] = (int*) malloc(m* sizeof(int));
        for (int j = 0; j < minA->length; ++j) {
            perturbationVectors[i][twoM[minA->data[j]].i] = twoM[minA->data[j]].r;
        }
    }

    return perturbationVectors;
}

double *LSH_search(int dim, int l, int m, double w, double ***hashTables, HashBucket *buckets, double *query) {
    double *result = (double *) malloc(sizeof(double) * dim);
    int **hashVal = calculateHashValues(dim, l, m, w, hashTables, query);
    double distance = MAXDOUBLE;
    const int NUM_VECTORS = 50;

    HashBucket *ite = buckets;

    while (ite->next) {
        if (compareHashValues(l, m, hashVal, ite->hashValues)) {
            //do probing here

            int **perturVectors = probing(NUM_VECTORS, dim, l, m, w, query, hashTables[0]);

            double localDistance = search(dim, ite, query, result);
            distance = distance > localDistance? localDistance: distance;

            //add perturVectors to hashVal => go through bucket to find the bucket again
            break;
        }
        ite = ite->next;
    }

    for (int i = 0; i < l; ++i) {
        free(hashVal[i]);
    }

    printf("Closest distance: %f", distance);

    free(hashVal);
    return result;
}



