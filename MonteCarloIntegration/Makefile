CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: out.txt plot_error_Ndep.svg

out.txt: main
	./$< 1> $@ 2> out_error_Ndep.txt

plot_error_Ndep.svg: plot_error_Ndep.gpi out_error_Ndep.txt
	gnuplot $< > $@

main: PlainMCIntegrator.o main.o

clean:
	$(RM) *.o main *.svg out*
