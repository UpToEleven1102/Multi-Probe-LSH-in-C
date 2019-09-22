#!/bin/bash
gcc main.c lib/utils.h lib/utils.c -O3 lsh_hetero -std=gnu99 -lm
./lsh_hetero