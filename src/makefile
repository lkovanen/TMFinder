CC = g++
CFLAGS = -O4 -Wall

all: tmf

tmf: main.o events.o edges.o tsubgraph.o subnets.o binner.h motif.o progress_counter.o bin_limits.o
	mkdir -p ../bin
	${CC} ${CFLAGS} -o  ../bin/tmf main.o tsubgraph.o subnets.o events.o edges.o motif.o progress_counter.o bin_limits.o -lstdc++ -lbliss

main.o: events.o tsubgraph.o main.cc subnets.o
	${CC} ${CFLAGS} -c main.cc

tsubgraph.o: tsubgraph.h tsubgraph.cc 
	${CC} ${CFLAGS} -c tsubgraph.cc

events.o: events.h events.cc
	${CC} ${CFLAGS} -c events.cc

edges.o: edges.h edges.cc
	${CC} ${CFLAGS} -c edges.cc

subnets.o: subnets.h subnets.cc
	${CC} ${CFLAGS} -c subnets.cc

bin_limits.o: bin_limits.h bin_limits.cc
	${CC} ${CFLAGS} -c bin_limits.cc

motif.o: motif.h motif.cc
	${CC} ${CFLAGS} -c motif.cc

progress_counter.o: progress_counter.h progress_counter.cc
	${CC} ${CFLAGS} -c progress_counter.cc

clean:
	rm -f ../bin/tmf main.o events.o edges.o tsubgraph.o subnets.o motif.o progress_counter.o bin_limits.o
