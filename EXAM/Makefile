CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: Exam.pdf

Exam.pdf: Exam.tex plot_real.tex plot_imag.tex
	pdflatex $<
	pdflatex $<

plot_real.tex plot_imag.tex: out.txt Makefile
	gnuplot -e "\
        set term cairolatex pdf;\
        set out 'plot_real.tex';\
        set xlabel 'x';\
        set ylabel 'y';\
	set key left top;\
        plot \
                'out.txt' i 0 u 1:3 with lines title 'Real exp from ODE',\
                'out.txt' i 0 u 1:5 with lines dashtype (10,20) linewidth 3 title 'Exp(x) from math.h';\
	"
	gnuplot -e "\
	set term cairolatex pdf;\
	set out 'plot_imag.tex';\
	set xlabel 'x';\
	set ylabel 'y';\
	set key font ',9';\
	plot \
		'out.txt' i 1 u 2:3 with lines title 'Re\{exp(ix)\} along imag-axis from ODE',\
		'out.txt' i 1 u 2:5 with lines dashtype (10,20) linewidth 3 title 'Exp(-x) from math.h',\
		'out.txt' i 1 u 2:4 with lines title 'Im\{exp(ix)\} along imag-axis from ODE',\
		'out.txt' i 2 u 1:3 with lines title 'Re\{exp(ix)\} along real axis from ODE',\
		'out.txt' i 2 u 1:4 with lines title 'Im\{exp(ix)\} along real axis from ODE',\
		'out.txt' i 2 u 1:5 with lines dashtype (10,20) linewidth 3 title 'Cos(x) from math.h',\
		'out.txt' i 2 u 1:6 with lines dashtype (10,20) linewidth 3 title 'Sin(x) from math.h';\
        "

out.txt: main
	./main 2> out.txt

main: main.o rkstep23.o ODE_driver.o

clean:
	$(RM) *.o main plot_real.tex plot_imag.tex *.pdf *.log *.aux out.txt
