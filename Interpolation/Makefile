CFLAGS = -Wall -std=gnu99
LDLIBS = -lm

all: plot_linear.svg plot_linear_integ.svg plot_quad.svg out_quad_parameters.txt plot_quad_deriv_int.svg

plot_linear.svg: plot_linear.gpi out_linear.txt
	gnuplot $< > $@

plot_linear_integ.svg: plot_linear_integ.gpi out_linear.txt
	gnuplot $< > $@

out_linear.txt: main-linterp data_linear.txt
	cat data_linear.txt | ./main-linterp $$(wc --lines data_linear.txt) < data_linear.txt 1> out_linear.txt

main-linterp: main-linterp.o linterp.o linterp_integ.o

plot_quad.svg: plot_quad.gpi out_quad.txt
	gnuplot $< > $@

plot_quad_deriv_int.svg: plot_quad_deriv_int.gpi out_quad.txt
	gnuplot $< > $@

out_quad.txt out_quad_parameters.txt: main-quad
	cat data_quad.txt | ./main-quad $$(wc --lines data_quad.txt) < data_quad.txt 1> out_quad.txt 2> out_quad_parameters.txt

main-quad: main-quad.o qspline.o


clean:
	$(RM) *.o out* main-linterp *.svg main-quad
