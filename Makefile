all: math
	gcc -Wall -c filter.c
	gcc -Wall -lm filter.o main.c -o filter-test
	gcc -Wall -lm cheby.c -o cheby
	gcc -Wall -lm math.c -o math

clean:
	\rm filter.o filter-test cheby
