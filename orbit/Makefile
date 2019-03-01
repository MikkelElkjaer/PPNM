CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = -lm $$(gsl-config --libs)

all: plot_part1.svg plot_part2.svg

plot_part1.svg: plot_part1.gpi plot_part1.txt
	gnuplot $< >$@

plot_part2.svg: plot_part2.gpi plot_part2.txt
	gnuplot $< >$@

plot_part1.txt plot_part2.txt: main
	./main 1> plot_part1.txt 2> plot_part2.txt

main: part1.o part2.o

clean:
	$(RM) *.o *.svg main plot_part1.txt plot_part2.txt
