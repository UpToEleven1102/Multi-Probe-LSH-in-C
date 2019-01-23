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

struct A *extractMinHeapLinkedList(struct Heap *_this) {
    struct A *minPtr, *ite;

    ite = _this->head;
    minPtr = ite;

    while (ite != NULL) {
        if(ite->score < minPtr->score) {
            minPtr = ite;
        }
        ite = ite->next;
    }

    return minPtr;
}
