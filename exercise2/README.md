## Exercise 2

This folder contains 2 directories: 

- `part1` contains all the files used to measure scalability over the size of the matrix (from 2000 to 20000 with step 1000) at fixed number of cores (64). It is divided in `single-precision` and `double-precision`, each of them contains:
	
	- the `gemm.c` file;
	
	- the folders `mkl` and `oblas`, each of them containing the batch files and the output files for that library;

	- the `Makefile`, set to `DUSE_FLOAT` or `DUSE_DOUBLE` for each case, which can be called with `make` and will produce the two executables in the proper folder;

- `part2` contains all the files used to measure scalability over the number of cores (from 1 to 64 with step 1) at fixed size (12000). The structure is the same of the previous folder.

The batch and output files are named so that you can easily understand what benchmark it corresponds to, following the order of the directories: for example, `p1-sp-mkl-cl.sh` is the batch file of part1, single precision, MKL library, with close binding policy.
