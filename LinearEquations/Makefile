CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: out_QR_decomp.txt

out_QR_decomp.txt: main
	./$< > $@

main: main.o QR_GS_decomp.o QR_GS_solve.o QR_GS_inverse.o

clean:
	$(RM) *.o main out*
