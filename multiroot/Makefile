CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: plot_wavefunction.svg plot_energyconvergence.svg out_Rosenbrock.txt

plot_wavefunction.svg: plot_wavefunction.gpi hydrogen_data.txt
	gnuplot $< > $@

plot_energyconvergence.svg: hydrogen_energies.txt Makefile
	gnuplot -e "\
	set terminal svg background 'white';\
	set yrange [-0.55: -0.4] ;\
	plot 'hydrogen_energies.txt' index 0 using 1:2 title 'Energy as fct of rmax' \
	, 'hydrogen_energies.txt' index 0 using 1:(-0.5) with lines title 'Ground state Energy';\
	" > $@

hydrogen_data.txt: main-hydrogen
	./$< 2> hydrogen_data.txt

main-hydrogen: main-hydrogen.o schroe_sol.o

out_Rosenbrock.txt: main-ros
	./$< 1> $@

clean:
	$(RM) *.o main-hydrogen main-ros hydrogen_data.txt plot_wavefunction.svg plot_energyconvergence.svg hydrogen_energies.txt out_Rosenbrock.txt
