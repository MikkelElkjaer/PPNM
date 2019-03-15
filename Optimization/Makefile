CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: out_Rosenbrock.txt plot_lsf.svg

out_Rosenbrock.txt: main-Ros
	./$< 1> $@

plot_lsf.svg: plot_lsf_data.txt Makefile
	gnuplot -e "\
	set terminal svg background 'white';\
	plot '$<' using 1:2:3 with errorbars title 'Experimental Data',\
	'$<' using 1:4 with lines smooth acsplines title 'Fit';\
	" > $@

plot_lsf_data.txt: main-lsf
	./main-lsf $$(wc --lines decay.data) < decay.data 2>plot_lsf_data.txt


clean:
	$(RM) *.o main-Ros out_Rosenbrock.txt out_minimization.txt plot_lsf_data.txt main-lsf plot_lsf.svg
