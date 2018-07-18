all:
	gcc pointsource.c Utils.c -fopenmp -lm -o ripple.x

debug:
	gcc pointsource.c Utils.c -fopenmp -lm -o ripple.debug -g
clean:
	rm *.x

