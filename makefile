all:
	gcc tank.c Utils.c -fopenmp -lm -o ripple.x
	tau_cc.sh htank.c Utils.c -fopenmp -lm -o hripple.x

debug:
	gcc tank.c Utils.c -fopenmp -lm -o ripple.debug -g
clean:
	rm *.x

