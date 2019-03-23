CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: report.pdf

report.pdf: report.tex plot-cairo.tex
	pdflatex $<
	pdflatex $<

plot-cairo.tex: out.txt Makefile
	gnuplot -e "\
	set term cairolatex pdf;\
	set out '$@';\
	set xlabel 'x';\
	set ylabel 'erf(x)';\
	plot \
		'out.txt' u 1:2 with lines title 'My Error Function',\
	"

out.txt: main
	./$< -3 3 0.05 > $@

main: main.o my_errfct.o

clean:
	$(RM) *.o main out.txt report.pdf plot-cairo.tex report.log report.aux plot-cairo.pdf
