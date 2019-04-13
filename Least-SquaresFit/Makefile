CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: plot.svg

plot.svg: plot.gpi out_lsfit.txt out_fit_data.txt
	gnuplot $< > $@

out_lsfit.txt out_fit_data.txt: main
	./main $$(wc --lines data.txt) < data.txt 1> out_lsfit.txt 2> out_fit_data.txt

main: main.o lsfit.o QR_GS_decomp.o QR_GS_solve.o QR_GS_inverse.o


clean:
	$(RM) out* main *.o *.svg
