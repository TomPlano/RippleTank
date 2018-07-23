all:
	gcc tank.c Utils.c -fopenmp -lm -o ripple.x

debug:
	gcc tank.c Utils.c -fopenmp -lm -o ripple.debug -g
clean:
	rm *.x

