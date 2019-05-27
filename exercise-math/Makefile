Exercise2Result.txt: exercise2
	./exercise2 > $@

exercise2: exercise2.o
	gcc exercise2.o -o exercise2 -lm

exercise2.o: exercise2.c
	gcc -c exercise2.c -o exercise2.o 

clean:
	rm -f *.o exercise2
