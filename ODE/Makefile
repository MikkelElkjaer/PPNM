CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: plot_cubic_root.svg plot_Lidocaine.svg

plot_Lidocaine.svg: plot_Lidocaine.gpi out_Lidocaine.txt
	gnuplot $< > $@

plot_cubic_root.svg: plot_cubic_root.gpi out.txt
	gnuplot $< > $@

out.txt: main
	./main 2> out.txt 1> out_Lidocaine.txt

main: main.o rkstep23.o ODE_driver.o

clean:
	$(RM) *.o main out.txt *.svg out_Lidocaine.txt
