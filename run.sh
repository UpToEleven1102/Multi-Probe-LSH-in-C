#!/bin/bash
gcc main.c lib/utils.h lib/utils.c -o lsh_hetero -std=gnu99 -lm
./lsh_hetero