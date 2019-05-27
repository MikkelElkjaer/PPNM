CFLAGS = -Wall -std=gnu99 -O1
LDLIBS = -lm 
CC = gcc

.PHONEY: default
default: out.txt
	cat $^

out.txt: exercise-epsilon
	./$< > $@

exercise-epsilon: exercise-epsilon.o eps_part2.o eps_part3.o

clean:
	rm -f *.o out*
	find -type f -executable -delete
