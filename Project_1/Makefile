# This target will compile all the file.
all: pdeSolver

pdeSolver.o: main.o partialDifferential.o utils.o
	gcc -c main.o partialDifferential.o utils.o -lm -o pdeSolver

main.o: ./src/main.c
	gcc -c ./src/main.c

partialDifferential.o: ./src/partialDifferential.c
	gcc -c ./src/partialDifferential.c

utils.o: ./src/utils.c
	gcc -c ./src/utils.c

clean:
	rm -rf *o pdeSolver
