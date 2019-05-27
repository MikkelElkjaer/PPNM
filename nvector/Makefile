CFLAGS := -Wall -std=gnu99
LDLIBS := -lm
CC = gcc

all: out.txt
	cat out.txt

out.txt: checker
	./$< > $@

checker: nvector.o double_equal.o

clean:
	$(RM) *.o out.txt
