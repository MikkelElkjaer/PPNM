CFLAGS := -Wall -std=gnu99
LDLIBS := -lm
plot.svg: plot.gpi data.txt #scalable vector graphics  - other option is png... but use svg! "gpi" er mere eller mindre accepteret standard for gnuplots...
	gnuplot plot.gpi

data.txt: main
	./main > $@

