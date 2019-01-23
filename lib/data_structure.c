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

    if(isEqualStructA(ite, ele)) {
        _this->head = _this->head->next;

        free(ite->data);
        free(ite);
        return;
    }

    while(ite->next != NULL){
        if(isEqualStructA(ite->next, ele)) {
            struct A *prev = ite->next;
            ite->next = ite->next->next;

            free(prev->data);
            free(prev);

            return;
        }
        ite = ite->next;
    }
}