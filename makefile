CC=g++

CPPFLAGS = -g -Wall -O3


boss: main.o scaffoldgraph.o scaffolding.o
	$(CC) -o $@ $^ ./lp/liblpsolve55.a -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lm -ldl -lz

main.o: main.cpp scaffoldgraph.h scaffolding.h
	$(CC) -c main.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
scaffoldgraph.o: scaffoldgraph.cpp scaffoldgraph.h
	$(CC) -c scaffoldgraph.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
scaffolding.o: scaffolding.cpp scaffoldgraph.h scaffolding.h
	$(CC) -c scaffolding.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	

all: boss
	rm -f *.o

clean:
	rm -f *.o
	rm boss
