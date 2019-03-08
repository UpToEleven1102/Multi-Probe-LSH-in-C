//
// Created by huyen on 1/23/19.
//

#ifndef LSH_PROBING_DATA_STRUCTURE_H
#define LSH_PROBING_DATA_STRUCTURE_H

#include "utils.h"

//struct for value z in first 2-m array
struct Z {
    double x;
    int i;
    int r;
};


//linked list struct for perturbation sets
struct A {
    int *data;
    double score;
    int length;
    struct A *next;

    void (*calculateScore)(struct A*, struct Z*);
    bool (*isValid)(struct A*, int);
    struct A* (*shift)(struct A*);
    struct A* (*expand)(struct A*);
};
void calculateScoreA(struct A*, struct Z*);

struct A *expandA(struct A *_this);

struct A *shiftA(struct A *_this);

bool isValidA(struct A *_this, int twoM);

struct Heap {
    struct A *head;
    void (*add)(struct Heap*, struct A*);
    void (*remove)(struct Heap*, struct A*);
    struct A* (*extractMin)(struct Heap*);
};

void addHeapLinkedList(struct Heap*, struct A*);

void removeHeapLinkedList(struct Heap*, struct A*);

struct A * extractMinHeapLinkedList(struct Heap*);

#endif //LSH_PROBING_DATA_STRUCTURE_H
