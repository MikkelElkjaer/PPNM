CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: project.pdf

project.pdf: project.tex plot-cairo.tex
	pdflatex $<
	pdflatex $<

plot-cairo.tex: out.txt Makefile
	gnuplot -e "\
	set term cairolatex pdf;\
	set out '$@';\
	set xlabel 'x';\
	set ylabel 'y';\
	plot \
		'out.txt' u 1:2 with lines title 'My Exp',\
                'out.txt' u 1:3 with lines dashtype (10,20) linewidth 3 title 'Exp(x) from math.h';\
        "

out.txt: main
	./$< > $@

main: main.o myexp.o

clean:
	$(RM) main *.o out.txt project.pdf plot-cairo.tex project.log project.aux plot-cairo.pdf
