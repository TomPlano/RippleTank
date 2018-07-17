all:
	gcc ripple.c Utils.c -fopenmp -lm -o ripple.x

debug:
	gcc ripple.c Utils.c -fopenmp -lm -o life.debug -g
clean:
	rm *.x

