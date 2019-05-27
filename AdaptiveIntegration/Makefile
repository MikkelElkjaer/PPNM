CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: out.txt

out.txt: main
	./$< > $@

main: Integrator.o main.o

clean:
	$(RM) *.o out.txt main
