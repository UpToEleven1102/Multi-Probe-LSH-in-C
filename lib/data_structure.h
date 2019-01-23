//
// Created by huyen on 1/23/19.
//

#ifndef LSH_PROBING_DATA_STRUCTURE_H
#define LSH_PROBING_DATA_STRUCTURE_H

//struct for value z in first 2-m array
struct Z {
    double x;
    int i;
    int r;
};


//linked list struct for perturbation sets
struct A {
    double *data;
    int length;
    struct A *next;
};

struct Heap {
    struct A *head;
    void (*add)(struct Heap*, struct A*);
    void (*remove)(struct Heap*, struct A*);
};

void addHeapLinkedList(struct Heap*, struct A*);

void removeHeapLinkedList(struct Heap*, struct A*);

#endif //LSH_PROBING_DATA_STRUCTURE_H
