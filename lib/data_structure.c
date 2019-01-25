//
// Created by huyen on 1/23/19.
//

#include <stddef.h>
#include <malloc.h>
#include "data_structure.h"
#include "utils.h"

bool isEqualStructA(struct A *_A1, struct A *_A2) {
    if (_A1->length != _A2->length)
        return false;
    for (int i = 0; i < _A1->length; ++i) {
        if (_A1->data[i] != _A2->data[i]) {
            return false
        }
    }
    return true;
}

bool isValidA(struct A *_this, int twoM) {
    for (int i = 0; i < _this->length; ++i) {
        if (_this->data[i] > twoM)
            return false;

        int negJ = twoM +1 -_this->data[i];
        for (int j = 0; j < _this->length; ++j) {
            if(_this->data[j] == negJ)
                return false;
        }
    }
    return true;
}


struct A *expandA(struct A *_this){
    struct A *expanded = (struct A*)malloc(sizeof(struct A));
    expanded->score = 0;
    expanded->length = _this->length +1;
    expanded->next = NULL;
    expanded->calculateScore = calculateScoreA;
    expanded->isValid = isValidA;
    expanded->shift = shiftA;
    expanded->expand = expandA;

    expanded->data = (int*)malloc(_this->length+1 *sizeof(int));
    for (int i = 0; i < _this->length; ++i) {
        expanded->data[i] = _this->data[i];
    }
    expanded->data[_this->length] = _this->data[_this->length-1]+1;

    return expanded;

}

struct A *shiftA(struct A *_this){
    int *data = (int *)malloc(_this->length*sizeof(int));

    for (int i = 0; i < _this->length-1; ++i) {
        data[i] = _this->data[i];
    }

    data[_this->length - 1] = _this->data[_this->length - 1] +1;

    struct A *shifted = (struct A *)malloc(sizeof(struct A));
    shifted->data = data;
    shifted->score = 0;
    shifted->length = _this->length;
    shifted->next = NULL;
    shifted->calculateScore = calculateScoreA;
    shifted->isValid = isValidA;
    shifted->shift = shiftA;
    shifted->expand = expandA;

    return shifted;
}


void calculateScoreA(struct A *_this, struct Z *zs) {
    _this->score = 0;
    for (int i = 0; i < _this->length; ++i) {
        _this->score += zs[_this->data[i]].x * zs[_this->data[i]].x;
    }
}


void addHeapLinkedList(struct Heap *_this, struct A *ele) {
    if (_this->head == NULL) {
        _this->head = ele;
    } else {
        ele->next = _this->head;
        _this->head = ele;
    }
}

void removeHeapLinkedList(struct Heap *_this, struct A *ele) {
    struct A *ite = _this->head;

    if (isEqualStructA(ite, ele)) {
        _this->head = _this->head->next;

        free(ite->data);
        free(ite);
        return;
    }

    while (ite->next != NULL) {
        if (isEqualStructA(ite->next, ele)) {
            struct A *prev = ite->next;
            ite->next = ite->next->next;

            free(prev->data);
            free(prev);

            return;
        }
        ite = ite->next;
    }
}


//note: get the lowest score of the heap or get the top element
struct A *extractMinHeapLinkedList(struct Heap *_this) {
//    struct A *minPtr, *ite;
//
//    ite = _this->head;
//    minPtr = ite;
//
//    while (ite != NULL) {
//        if(ite->score < minPtr->score) {
//            minPtr = ite;
//        }
//        ite = ite->next;
//    }
//
//    return minPtr;
    if(_this->head == NULL)
        return NULL;

    struct A* topNode = _this->head;

    _this->head = _this->head->next;

    return topNode;
}
