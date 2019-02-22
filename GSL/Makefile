CFLAGS := -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS := -lm `gsl-config --libs`

all: out.txt
	cat $<

out.txt: main input.matrix
	cat input.matrix | ./main > $@

main: main.o


clean:
	$(RM) main out.txt *.o
