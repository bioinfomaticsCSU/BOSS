CC=g++

CPPFLAGS = -g -Wall -O3


boss:	main.o scaffoldgraph.o scaffolding.o
	$(CC) -o $@ $^ ./lp/liblpsolve55.a -lm -ldl -I include/ -L lib/ -lbamtools

main.o: main.cpp scaffoldgraph.h scaffolding.h
	$(CC) -c main.cpp -lm -ldl -I include/ -L lib/ -lbamtools
	
scaffoldgraph.o: scaffoldgraph.cpp scaffoldgraph.h
	$(CC) -c scaffoldgraph.cpp -lm -ldl -I include/ -L lib/ -lbamtools
	
scaffolding.o: scaffolding.cpp scaffoldgraph.h scaffolding.h
	$(CC) -c scaffolding.cpp -lm -ldl -I include/ -L lib/ -lbamtools
	

all: boss

clean:
	rm -f *.o
	rm boss
