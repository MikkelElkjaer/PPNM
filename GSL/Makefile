CFLAGS := -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS := -lm `gsl-config --libs`

all: out.txt airy.svg
	cat $<

out.txt airy.txt: main input.matrix
	cat input.matrix | ./main 1> $@ 2> airy.txt

airy.svg: airy.gpi airy.txt
	gnuplot $< >$@

main: main.o airy.o


clean:
	$(RM) main out.txt *.o airy.svg airy.txt
