# High Performance Triangle Counting 
Given an adjacency matrix A in matrix market format, our code calculates in how many triangles is each node adjacent to.
This is done using OMP, OpenCilk and pthreads shared memory parallelisation.

For V3:
1)Serial implementation: V3serial.c
2)Parallel with openMP: V3omp.c
3)Parallel with openCilk: V3cilk.c

For V4:
1)Serial implementation: V4serial.c
2)Parallel with pthreads: V4pthreads.c
3)Parallel with openMP: V4omp.c
4)Parallel with openCilk: V4cilk.c

Finally, a Makefile to make the executables is included. A report (in greek) for more comments and insights regarding the implementations and the code performance is also included.

Co-written with [George Tsoumplekas](https://github.com/GeorgeTsoumplekas)


