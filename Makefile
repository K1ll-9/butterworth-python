all:
	gcc -Wall -c filter.c
	gcc -Wall -lm filter.o main.c -o filter-test
