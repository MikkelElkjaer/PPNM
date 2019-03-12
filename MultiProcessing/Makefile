CFLAGS = -fopenmp -pthread -Wall -std=gnu99
LDLIBS = -lm -pthread

all: plot.svg

plot.svg: out.txt Makefile
	gnuplot -e "\
	set term svg backg 'white';\
	set logscale x;\
	plot \
		'out.txt' index 0 u 1:3 title 'Residuals using pthreads',\
		'out.txt' index 0 u 1:(0) with lines dashtype (5,10) notitle,\
		'out.txt' index 1 u 1:3 title 'Residuals using OpenMP';\
	" > $@

out.txt: main
	./$< > $@

clean:
	$(RM) *.o main plot.svg out.txt
