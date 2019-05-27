CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: out.txt

out.txt: main
	./main > $@

main: Newton.o QR_GS_decomp.o QR_GS_solve.o functions.o Quasi_Newton.o main.o


clean:
	$(RM) *.o out* main
