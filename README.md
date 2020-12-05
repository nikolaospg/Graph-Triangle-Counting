# Parallel-Systems-Exercise-1

This directory contains all the code written for the first exercise for the course 'Parallel and Distributed Systems', semester 7 at ECE, AUTh.

Given an adjacency matrix A in matrix market format, our code calculates in how many triangles is each node adjacent to.
This is done using two different methods, namely V3 using a triple for loop and V4 that takes advantage of A to the square matrix thus reducing complexity.

Both of these methods are implemented both serially and in parallel.

For V3:
1)Serial implementation: V3serial.c
2)Parallel with openMP: V3omp.c
3)Parallel with openCilk: V3cilk.c

For V4:
1)Serial implementation: V4serial.c
2)Parallel with pthreads: V4pthreads.c
3)Parallel with openMP: V4omp.c
4)Parallel with openCilk: V4cilk.c

Finally, a Makefile to make the executables is included.

Authors: 
1)Nikolaos Papageorgiou, nikolaospg@ece.auth.gr
2)Georgios Tsoumplekas, gktsoump@ece.auth.gr

