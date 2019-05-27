CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc
all: plot.svg

plot.svg: out.txt Makefile
	echo '\
	set terminal svg background "white";\
	set out "$@";\
	set title "Curve Fitting by Neural Network";\
	set xlabel "x";\
	set ylabel "y";\
	plot\
	"$<" index 0 with points pointtype 2 title "data"\
	,"$<" index 1 with lines linestyle 2 title "ANN Fit"\
	'|gnuplot

out.txt: main
	./$< 1> $@ 2> log.txt


main: main.o ann.o

clean:
	$(RM) *.o out.txt main *.svg log.txt
