CFLAGS = -Wall -std=gnu99
CC = gcc
LDLIBS = -lm

all: plot.svg

plot.svg: plot.gpi out.io.txt out.cmd.txt
	gnuplot $< > $@

out.io.txt: main-stdio input.data.txt
	./main-stdio < input.data.txt > $@

out.cmd.txt: main-cmdline input.data.txt
	./main-cmdline `cat input.data.txt` > $@

#exploit built-in defaults for main-stdio and main-cmdline?
main-stdio: main-stdio.o

main-cmdline: main-cmdline.o


clean:
	$(RM) *.o main-stdio main-cmdline test.io.out.txt test.cmd.out.txt plot.svg out.io.txt out.cmd.txt
