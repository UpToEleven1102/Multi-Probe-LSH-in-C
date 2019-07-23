start:
	./run.sh

build:
	gcc main.c lib/utils.h lib/utils.c -o lsh_hepmass -std=gnu99 -lm

run:
	./lsh
