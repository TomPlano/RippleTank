all:
	gcc betterripple.c Utils.c -fopenmp -lm -o ripple.x

debug:
	gcc betterripple.c Utils.c -fopenmp -lm -o ripple.debug -g
clean:
	rm *.x

