CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: plot_energy.svg

main: main.o ln_over_sqrt.o norm.o H_HO.o

plot_energy.svg: plot_energy.gpi out.txt
	gnuplot $< > $@

out.txt: main
	./main 2> $@


clean:
	$(RM) *.o main out.txt plot_energy.svg
