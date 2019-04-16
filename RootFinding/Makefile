CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: out.txt

out.txt: main main-gslrootfinding
	./main > $@
	./main-gslrootfinding

main: main.o Newton_with_Jacobian.o QR_GS_decomp.o QR_GS_solve.o functions.o Newton.o Newton_refined_linesearch.o


clean:
	$(RM) out* main *.o main-gslrootfinding
