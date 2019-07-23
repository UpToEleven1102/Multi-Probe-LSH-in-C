#!/bin/bash
gcc main.c lib/utils.h lib/utils.c -o lsh -std=gnu99 -lm
./lsh