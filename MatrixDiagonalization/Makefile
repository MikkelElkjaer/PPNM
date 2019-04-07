CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: out_sweep.txt time_sweep.txt out_k_eigvals.txt

out_sweep.txt time_sweep.txt: main_sweep Makefile
	./main_sweep 5 1> out_sweep.txt 2> /dev/null ;\
	cat /dev/null > time_sweep.txt
	for n in `seq 100 20 200`; do \
	\time --format "$$n %U" --append --output time_sweep.txt ./main_sweep $$n 1> /dev/null ;\
	done

out_k_eigvals.txt: main_k_eigvals
	./main_k_eigvals 5 3 1 1> out_k_smallest_eigvals.txt 2> /dev/null
	./main_k_eigvals 5 3 -1 1> out_k_largest_eigvals.txt 2> /dev/null
	./main_k_eigvals 100 1 1 1> /dev/null 2> time_compare_lowest_eigval.txt
	./main_k_eigvals 100 100 1 1> /dev/null 2> time_compare_full_diag.txt

main_sweep: main_sweep.o Jacobi_sweep.o

main_k_eigvals: main_k_eigvals.o Jacobi_k_eigvals.o Jacobi_sweep.o

clean:
	$(RM) *.o main_sweep out* time* main_k_eigvals
