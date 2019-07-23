#!/bin/bash
gcc main.c lib/utils.h lib/utils.c -o lsh_tlc -std=gnu99 -lm
./lsh