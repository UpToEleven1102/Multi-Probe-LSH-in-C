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
    int **hashValues = calculateHashValues(dim, l, m, w, centroid, hashTables, ele);

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

HashBucket *LSH(int dim, int n_data, int l, int m, double w, double ***hashTables, double *data, HashBucket *buckets,
                double *centroid) {
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

double calculateScore(const int *a0, int length, struct pairZ zs[]) {
    printf("debug \n");
    for (int j = 0; j < length; ++j) {
        printf("%d   ", a0[j]);
    }

    printf("\n");

    double score = 0;
    for (int i = 0; i < length; ++i) {
//        score += zs[a0[i]].x;
        score += zs[a0[i]].x * zs[a0[i]].x;
    }
    return score;
}

struct HeapTreeNode *minHeap(struct HeapTreeNode **heap) {
    if (*heap == NULL) {
        return NULL;
    }

    if ((*heap)->left == NULL) {
        struct HeapTreeNode *min = (struct HeapTreeNode *) malloc(sizeof(struct HeapTreeNode));
        min->score = (*heap)->score;
        min->length = (*heap)->length;
        min->data = (int *) malloc(min->length * sizeof(int));
        for (int i = 0; i < min->length; i++) {
            min->data[i] = (*heap)->data[i];
        }
        min->left = NULL;
        min->right = NULL;

        //remove heap
        if ((*heap)->right == NULL) {
            //no right child
            if ((*heap)->parent !=NULL) {
                (*heap)->parent->left = NULL;
            }
//            freeHeap(*heap);
            *heap = NULL;
        } else {
            if ((*heap)->parent != NULL) {
                (*heap)->right->parent = (*heap)->parent;
                (*heap)->parent->left = (*heap)->right;
                free((*heap)->data);
                free((*heap));
            } else {
//                freeHeap(*heap);
                (*heap)->right->parent = NULL;
                free((*heap)->data);
                *heap = (*heap)->right;
                (*heap)->parent = NULL;
            }
        }
        return min;
    } else {
        return minHeap(&(*heap)->left);
    }
}

//error might happen here
struct HeapTreeNode *insertEleHeap(struct HeapTreeNode **heap, struct HeapTreeNode *parent, struct HeapTreeNode **ele) {
    //pass double pointer to modify the heap outside the scope of the function
    if (*heap == NULL) {

        //ading parent to element
        (*ele)->parent = parent;
        *heap = *ele;
        return *ele;
    } else {
        if ((*ele)->score <= (*heap)->score) {
            (*heap)->parent = parent;
            (*heap)->left = insertEleHeap(&(*heap)->left, *heap, ele);
        } else {
            (*heap)->parent = parent;
            (*heap)->right = insertEleHeap(&(*heap)->right, *heap, ele);
        }
        return *heap;
    }
}

bool isValid(struct HeapTreeNode *ele, int twoM) {
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


struct HeapTreeNode *shiftHeap(struct HeapTreeNode *ele, struct pairZ *zs) {
    int *data = (int *) malloc(ele->length * sizeof(int));

    for (int i = 0; i < ele->length - 1; ++i) {
        data[i] = ele->data[i];
    }

    data[ele->length - 1] = ele->data[ele->length - 1] + 1;

    struct HeapTreeNode *shifted = (struct HeapTreeNode *) malloc(sizeof(struct HeapTreeNode));
    shifted->data = data;
    shifted->score = calculateScore(data, ele->length, zs);
    shifted->length = ele->length;
    shifted->left = NULL;
    shifted->right = NULL;
    shifted->parent = NULL;

    return shifted;
}

struct HeapTreeNode *expandHeap(struct HeapTreeNode *ele, struct pairZ *zs) {
    struct HeapTreeNode *expanded = (struct HeapTreeNode *) malloc(sizeof(struct HeapTreeNode));
    expanded->data = (int *) malloc(ele->length + 1 * sizeof(int));
    for (int i = 0; i < ele->length; ++i) {
        expanded->data[i] = ele->data[i];
    }
    expanded->data[ele->length] = ele->data[ele->length - 1] + 1;
    expanded->length = ele->length + 1;
    expanded->score = calculateScore(expanded->data, expanded->length, zs);
    expanded->left = NULL;
    expanded->right = NULL;
    expanded->parent = NULL;
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

//    for (int i = 0; i < 2 * m; ++i) {
//        printf("%f - %d \n", twoM[i].x, twoM[i].i);
//    }

    for (int i = 0; i < 2 * m; ++i) {
        for (int j = i + 1; j < 2 * m; ++j) {
            struct pairZ temp = twoM[i];
            if (twoM[j].x < temp.x) {
                twoM[i] = twoM[j];
                twoM[j] = temp;
            }
        }
    }


//    for (int i = 0; i < 2 * m; ++i) {
//        printf("%f - %d \n", twoM[i].x, twoM[i].i);
//    }

    //generate Heap
    int *a0 = (int *) malloc(sizeof(int));
    a0[0] = 0;

    struct HeapTreeNode *heap = &(struct HeapTreeNode) {.data = a0, .length=1, .score = calculateScore(a0, 1,
                                                                                                       twoM), .left=NULL, .right=NULL, .parent=NULL};
    //generate perturbation vectors
    for (int i = 0; i < numOfVectors; ++i) {
        struct HeapTreeNode *minA = NULL;
        do {
            //extract minHeap
            minA = minHeap(&heap);

            struct HeapTreeNode *shifted = shiftHeap(minA, twoM);
            insertEleHeap(&heap, NULL, &shifted);
            struct HeapTreeNode *expanded = expandHeap(minA, twoM);
            insertEleHeap(&heap, NULL, &expanded);
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
