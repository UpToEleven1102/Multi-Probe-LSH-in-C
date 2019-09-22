#!/bin/bash
gcc main.c lib/utils.h lib/utils.c -o lsh -O3 -std=gnu99 -lm
./lsh