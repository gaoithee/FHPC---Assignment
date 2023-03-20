## Exercise 1

This directory contains all the files needed to compile and run an executable file:

  - 2 makefiles, one to produce a executable which will run the code using a program that uses `MPI_Scatterv` to share the workload with all the MPI processes (this can be produced by using `make -f Makefile_scat`), and another that produces a faster version of the program which doesn't use `MPI_Scatterv` (the default one, can be produced by simply using `make`);

  - some folders: 
	
	  -`batch` contains the batch files used to submit jobs to Orfeo; 
	  
	  -`data` contains the collected data; 
	  
	  -`obj` collects the object files obtained from compilation; 
	  
	  -`snap` is the folder where the output images will be stored;
	  
	  -`src` contains all the source codes.

Going more deeply, `src` contains:

  - `w_init.c`, the file used to initialize a matrix, it initializes a matrix with random values (in serial or in parallel according to the value of `np`);

  - `grow_ordered.c`, updates a matrix using "ordered" method, it is serial with no parallel version;

  - `grow_static_scat.c`, updates a matrix using "static" method, it has both a serial and a parallel method, the latter using `Scatterv`;

  - `grow_static.c`, updates a matrix using "static" method, it has both a serial and a parallel method, the latter without using `Scatterv`;

  - `read_write_pgm.c`, contains the functions used to read and write pgm images.

To use the program, you must initialize a matrix (by calling `mpirun -np 1 main.x -i -k 100 -f file`, it will produce a *pgm* image of a matrix with size `100` called `file`), then you can call `mpirun -np 1 main.x -r -e 0 -s 5 -n 50 -f file` to evolve the initial status of `file` for `50` iterations, with `ordered` method, producing a snapshot each `5` iterations.
