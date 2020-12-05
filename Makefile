#Simple Makefile used to compile and run every implementation we had to.

#Make (all) command: compiles and produces executables for all of the implementations
#Make run command:	runs the executable of every simple implementation, with the proper variables
#Make clean command: removes all of the executable files created by the Make command.

CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang        
CFLAGS=-O3
CC=gcc



default:all

V3cilk:
	$(CILKCC) $(CFLAGS) -o V3cilk.out V3cilk.c -fcilkplus

V4cilk:
	$(CILKCC) $(CFLAGS) -o V4cilk.out V4cilk.c -fcilkplus	

V3omp:
	$(CC) $(CFLAGS) -o V3omp.out V3omp.c -fopenmp	

V4omp:
	$(CC) $(CFLAGS) -o V4omp.out V4omp.c -fopenmp	

V4pthreads:
	$(CC) $(CFLAGS) -o V4pthreads.out V4pthreads.c	-lpthread	

V3serial:
	$(CC) $(CFLAGS) -o V3serial.out V3serial.c

V4serial:
	$(CC) $(CFLAGS) -o V4serial.out V4serial.c

all: V3cilk V4cilk V3omp V4omp V4pthreads V4serial V3serial


#To use the run command, the user must give the values to two variables. These are the .mtx filename and the number of threads employed.

#e.g. make run fileName=com-Youtube.mtx threadNum=3
#This command runs the programms, they read the com-Youtube.mtx file and they employ 3 threads on the parallelization parts
run:
	CILK_NWORKERS=$(threadNum)
	./V3serial.out		$(fileName)
	./V3omp.out 		$(fileName) $(threadNum)
	./V3cilk.out		$(fileName)
	./V4serial.out		$(fileName)
	./V4pthreads.out 	$(fileName) $(threadNum)
	./V4omp.out 		$(fileName) $(threadNum)
	./V4cilk.out		$(fileName)
.PHONY: clean

clean:
	rm -f V3cilk.out V4cilk.out V3omp.out V4omp.out V4pthreads.out V4serial.out V3serial.out
