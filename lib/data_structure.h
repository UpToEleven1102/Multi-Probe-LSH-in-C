//
// Created by huyen on 1/23/19.
//

#ifndef LSH_PROBING_DATA_STRUCTURE_H
#define LSH_PROBING_DATA_STRUCTURE_H

#include "utils.h"

//struct for value z in first 2-m array
struct pairZ {
    double x;
    int i;
    int r;
};


//linked list struct for perturbation sets
struct PerturVector {
    int *data;
    double score;
    int length;
    struct PerturVector *next;
    struct PerturVector *prev;

    void (*calculateScore)(struct PerturVector*, struct pairZ*);
    bool (*isValid)(struct PerturVector*, int);
    struct PerturVector* (*shift)(struct PerturVector*);
    struct PerturVector* (*expand)(struct PerturVector*);
};
void calculateScoreA(struct PerturVector*, struct pairZ*);

struct PerturVector *expandA(struct PerturVector *_this);

struct PerturVector *shiftA(struct PerturVector *_this);

bool isValidA(struct PerturVector *_this, int twoM);

struct Heap {
    struct PerturVector *head;
    void (*add)(struct Heap*, struct PerturVector*);
    void (*remove)(struct Heap*, struct PerturVector*);
    struct PerturVector* (*extractMin)(struct Heap*);
};

void addHeapLinkedList(struct Heap*, struct PerturVector*);

void removeHeapLinkedList(struct Heap*, struct PerturVector*);

struct PerturVector *extractMinHeapLinkedList(struct Heap*);

#endif //LSH_PROBING_DATA_STRUCTURE_H
